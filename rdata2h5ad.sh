#!/usr/bin/env bash
# =============================================================================
# rdata2h5ad.sh — Convert .Rdata / .RData files to h5ad format
# =============================================================================
# Converts single-cell data stored in legacy Rdata (workspace dump) format to
# h5ad (AnnData). Unlike .rds (single object), .Rdata can contain multiple
# named objects — the script auto-detects the count matrix and optional
# metadata (PCA, t-SNE) using heuristics, with CLI overrides.
#
# Usage:
#   bash rdata2h5ad.sh <input.Rdata> [output.h5ad] [options]
#
# Options:
#   --list                 List all objects in the .Rdata and exit
#   --count-object NAME    Name of the count matrix object (auto-detected)
#   --filter-object NAME   Object whose rownames define which cells to keep
#   --pca-object NAME      PCA coordinates object name
#   --tsne-object NAME     t-SNE coordinates object name
#   --umap-object NAME     UMAP coordinates object name
#   --method mtx           Force MTX extraction (default; the only option)
#   --force                Overwrite existing output file
#   --verbose              Show detailed debug logs
#   --skip-deps-check      Skip the pre-flight dependency check
#
# Environment variables:
#   R2H5AD_CONDA_ENV       Name of conda env to use (auto-detected otherwise)
#
# Examples:
#   bash rdata2h5ad.sh my_data.Rdata
#   bash rdata2h5ad.sh my_data.Rdata output.h5ad --force
#   bash rdata2h5ad.sh my_data.Rdata --count-object exprs --filter-object pca
#   bash rdata2h5ad.sh my_data.Rdata --list
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_DIR="${SCRIPT_DIR}/R"

# ── Defaults ────────────────────────────────────────────────────────────────
FORCE=false
VERBOSE=false
SKIP_DEPS_CHECK=false
LIST_ONLY=false
INPUT_FILE=""
OUTPUT_FILE=""
PASSTHROUGH_ARGS=()

# ── Functions ───────────────────────────────────────────────────────────────

log() {
    local level="$1"; shift
    local msg="$*"
    case "$level" in
        INFO)  echo "[INFO] $msg" ;;
        WARN)  echo -e "\033[33m[WARN]\033[0m $msg" >&2 ;;
        ERROR) echo -e "\033[31m[ERROR]\033[0m $msg" >&2 ;;
        DEBUG) $VERBOSE && echo "[DEBUG] $msg" || true ;;
    esac
}

have_cmd() { command -v "$1" &>/dev/null; }

# Convert Windows path to WSL path
win_to_wsl() {
    local path="$1"
    if [[ "$path" =~ ^([A-Za-z]): ]]; then
        local drive="${BASH_REMATCH[1],,}"
        local rest="${path#?:}"
        rest="${rest//\\//}"
        echo "/mnt/${drive}${rest}"
    else
        echo "$path"
    fi
}

wsl_to_win() {
    local path="$1"
    if [[ "$path" =~ ^/mnt/([a-z]) ]]; then
        local drive="${BASH_REMATCH[1]^^}"
        local rest="${path#/mnt/?}"
        echo "${drive}:${rest}"
    else
        echo "$path"
    fi
}

# Resolve output path
resolve_output() {
    local input="$1"
    local specified="$2"
    if [ -n "$specified" ]; then
        echo "$specified"
    else
        local dir="$(dirname "$input")"
        local base="$(basename "$input")"
        local stem="${base%.*}"
        stem="${stem%.Rdata}"
        stem="${stem%.rdata}"
        stem="${stem%.RData}"
        stem="${stem%.rda}"
        [ "$dir" = "." ] && echo "${stem}.h5ad" || echo "${dir}/${stem}.h5ad"
    fi
}

# Conda activation (reuses r2h5ad.sh logic)
try_activate_conda() {
    local conda_sh_paths=(
        "${HOME}/miniforge3/etc/profile.d/conda.sh"
        "${HOME}/miniconda3/etc/profile.d/conda.sh"
        "${HOME}/anaconda3/etc/profile.d/conda.sh"
        "/opt/miniforge3/etc/profile.d/conda.sh"
        "/opt/conda/etc/profile.d/conda.sh"
    )
    for csh in "${conda_sh_paths[@]}"; do
        if [ -f "$csh" ]; then
            log DEBUG "Found conda at: $csh"
            # shellcheck disable=SC1090
            source "$csh" 2>/dev/null || continue
            local env_name="${R2H5AD_CONDA_ENV:-}"
            if [ -n "$env_name" ]; then
                conda activate "$env_name" 2>/dev/null && return 0
            fi
            for try_env in r2h5ad atacseq-archr singlecell base; do
                if conda activate "$try_env" 2>/dev/null; then
                    log DEBUG "Activated conda env: $try_env"
                    return 0
                fi
            done
            conda activate base 2>/dev/null || true
            return 0
        fi
    done
    log DEBUG "No conda found — using system R and Python"
    return 1
}

# Pre-flight checks
check_prerequisites() {
    # Try conda activation
    try_activate_conda || true

    if ! have_cmd Rscript; then
        log ERROR "Rscript not found on PATH"
        exit 1
    fi
    if ! have_cmd python3; then
        log ERROR "python3 not found on PATH"
        exit 1
    fi
}

# Detect conda env for R/Python tools
CONDA_RUN_R=
CONDA_RUN_PY=
if command -v conda &>/dev/null; then
  C="conda run -n ${R2H5AD_CONDA_ENV:-r2h5ad} --"
  $C Rscript -e 'cat("ok")' 2>/dev/null 1>&2 && CONDA_RUN_R="$C Rscript"
  $C python3 -c 'print(1)' 2>/dev/null 1>&2 && CONDA_RUN_PY="$C python3"
fi

# Validate h5ad
validate_h5ad() {
    local h5ad_path="$1"
    ${CONDA_RUN_PY:-python3} -c "import anndata; a=anndata.read_h5ad('${h5ad_path}'); print('Shape:',a.shape); print('obsm:',list(a.obsm.keys()) if len(a.obsm)>0 else 'none')" 2>&1 || true
}

# ── Argument Parsing ────────────────────────────────────────────────────────

while [ $# -gt 0 ]; do
    case "$1" in
        --list)            LIST_ONLY=true ;;
        --count-object)    shift; PASSTHROUGH_ARGS+=(count-object="$1") ;;
        --filter-object)   shift; PASSTHROUGH_ARGS+=(filter-object="$1") ;;
        --pca-object)      shift; PASSTHROUGH_ARGS+=(pca-object="$1") ;;
        --tsne-object)     shift; PASSTHROUGH_ARGS+=(tsne-object="$1") ;;
        --umap-object)     shift; PASSTHROUGH_ARGS+=(umap-object="$1") ;;
        --method)          shift ;;
        --force)           FORCE=true ;;
        --verbose)         VERBOSE=true ;;
        --skip-deps-check) SKIP_DEPS_CHECK=true ;;
        --help|-h)
            head -40 "$0" | grep -E '^#' | sed 's/^# \?//'
            exit 0 ;;
        -*)
            log ERROR "Unknown option: $1"
            exit 1 ;;
        *)
            if [ -z "$INPUT_FILE" ]; then INPUT_FILE="$1"
            elif [ -z "$OUTPUT_FILE" ]; then OUTPUT_FILE="$1"
            else log ERROR "Unexpected argument: $1"; exit 1
            fi
            ;;
    esac
    shift
done

[ -z "$INPUT_FILE" ] && { log ERROR "No input file specified"; usage; exit 1; }

INPUT_FILE=$(win_to_wsl "$INPUT_FILE")
OUTPUT_FILE=$(resolve_output "$INPUT_FILE" "$(win_to_wsl "${OUTPUT_FILE:-}")")

log INFO "======== rdata2h5ad: .Rdata → h5ad Conversion ========"
log INFO "Input:  $(wsl_to_win "$INPUT_FILE")"
log INFO "Output: $(wsl_to_win "$OUTPUT_FILE")"

[ -f "$INPUT_FILE" ] || { log ERROR "Input file not found: $INPUT_FILE"; exit 1; }

# ── List mode ───────────────────────────────────────────────────────────────
if [ "$LIST_ONLY" = true ]; then
    log INFO "Listing objects in .Rdata..."
    check_prerequisites
    ${CONDA_RUN_R:-Rscript} -e "
load('${INPUT_FILE}')
for (nm in ls()) {
    obj <- get(nm)
    cls <- class(obj)[1]
    if (is.data.frame(obj) || is.matrix(obj) || inherits(obj, 'Matrix')) {
        cat(sprintf('  %-25s  %-20s  %d x %d\n', nm, cls, nrow(obj), ncol(obj)))
    } else {
        sz <- round(object.size(obj) / 1e6, 1)
        cat(sprintf('  %-25s  %-20s  %.1f MB\n', nm, cls, sz))
    }
}
" 2>&1
    exit 0
fi

# ── Pre-flight ──────────────────────────────────────────────────────────────
check_prerequisites

if [ -f "$OUTPUT_FILE" ] && [ "$FORCE" != true ]; then
    log ERROR "Output file already exists: $OUTPUT_FILE"
    log INFO  "Use --force to overwrite"
    exit 1
fi

[ -f "$OUTPUT_FILE" ] && [ "$FORCE" = true ] && {
    log WARN "Overwriting existing file: $OUTPUT_FILE"
    rm -f "$OUTPUT_FILE"
}

# ── Create working directory ────────────────────────────────────────────────
WORK_DIR=$(mktemp -d -t rdata2h5ad.XXXXXXXXXX)
trap 'rm -rf "$WORK_DIR"' EXIT
log DEBUG "Work dir: $WORK_DIR"

# ── Step 1: Extract MTX + metadata via R ───────────────────────────────────
log INFO "[1/2] Extracting data from .Rdata..."

R_SCRIPT="${R_DIR}/convert_rdata.R"
[ -f "$R_SCRIPT" ] || { log ERROR "Missing R script: $R_SCRIPT"; exit 1; }

export R2H5AD_VERBOSE=$VERBOSE
R_ARGS=("$INPUT_FILE" "$WORK_DIR")
for arg in "${PASSTHROUGH_ARGS[@]}"; do
    R_ARGS+=("--$arg")
done
[ "$VERBOSE" = true ] && R_ARGS+=("--verbose=true")

R_OUTPUT=$(${CONDA_RUN_R:-Rscript} "$R_SCRIPT" "${R_ARGS[@]}" 2>&1) || {
    log ERROR "R conversion failed. Stderr/out:"
    echo "$R_OUTPUT" >&2
    exit 1
}

# Parse JSON summary (last line starting with '{')
R_SUMMARY=$(echo "$R_OUTPUT" | grep -E '^\s*\{' | tail -1)

if echo "$R_SUMMARY" | ${CONDA_RUN_PY:-python3} -c "import sys,json; d=json.load(sys.stdin); sys.exit(0 if d.get('status')=='success' else 1)" 2>/dev/null; then
    COUNT_OBJ=$(echo "$R_SUMMARY" | ${CONDA_RUN_PY:-python3} -c "import sys,json; print(json.load(sys.stdin).get('count_object','?'))")
    N_GENES=$(echo "$R_SUMMARY" | ${CONDA_RUN_PY:-python3} -c "import sys,json; print(json.load(sys.stdin).get('n_genes','?'))")
    N_CELLS=$(echo "$R_SUMMARY" | ${CONDA_RUN_PY:-python3} -c "import sys,json; print(json.load(sys.stdin).get('n_cells','?'))")
    log INFO "Count matrix: '$COUNT_OBJ' — $N_GENES genes x $N_CELLS cells"
else
    log WARN "Could not parse R summary. Raw output:"
    echo "$R_OUTPUT" | grep -E '^\s*\{' | head -1
fi

# ── Step 2: Assemble h5ad via Python ────────────────────────────────────────
log INFO "[2/2] Assembling h5ad..."

PY_SCRIPT="${R_DIR}/assemble_h5ad.py"
[ -f "$PY_SCRIPT" ] || { log ERROR "Missing Python script: $PY_SCRIPT"; exit 1; }

${CONDA_RUN_PY:-python3} "$PY_SCRIPT" "$WORK_DIR" "$OUTPUT_FILE" 2>&1 || {
    log ERROR "h5ad assembly failed"
    exit 1
}

# ── Validate ────────────────────────────────────────────────────────────────
log INFO "Validating output h5ad..."
VALIDATION=$(validate_h5ad "$OUTPUT_FILE" 2>&1) && VALIDATE_EXIT=0 || VALIDATE_EXIT=1
if [ $VALIDATE_EXIT -eq 0 ]; then
    log INFO "Validation passed:"
    echo "$VALIDATION" | while IFS= read -r line; do log INFO "  $line"; done
else
    log WARN "Validation returned warnings (file may still be usable):"
    echo "$VALIDATION" >&2
fi

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo 0)
OUTPUT_MB=$(python3 -c "print(round($OUTPUT_SIZE/1e6,1))")
log INFO "Output size: $(printf "%'d" "$OUTPUT_SIZE") bytes (${OUTPUT_MB} MB)"

log INFO "======== Conversion complete ========"
echo "$OUTPUT_FILE"
