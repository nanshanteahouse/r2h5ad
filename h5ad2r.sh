#!/usr/bin/env bash
# =============================================================================
# h5ad2r.sh — Convert h5ad (AnnData) files back to R-native formats
# =============================================================================
#
# Converts AnnData .h5ad files to .rds / .qs (Seurat) or .Rdata (workspace dump).
# Reverse of r2h5ad.sh.
#
# Output formats (auto-detected by extension):
#   .rds                  — Single Seurat object (via saveRDS)
#   .qs                   — Single Seurat object (via qs::qsave)
#   .Rdata / .RData / .rda — Workspace dump (counts, cell_meta, reductions...)
#
# Three conversion paths:
#   SeuratDisk reverse  — h5ad → .h5Seurat → LoadH5Seurat → Seurat → .rds/.qs
#   MTX intermediate    — h5ad → extract_h5ad.py → MTX+CSV → assemble_seurat.R → .rds/.qs
#   Rdata export        — h5ad → extract_h5ad.py → MTX+CSV → assemble_rdata.R → .Rdata
#
# Usage:
#   bash h5ad2r.sh <input.h5ad> [output] [options]
#
# Options:
#   --method seuratdisk|mtx   Force conversion method
#   --force                   Overwrite existing output file
#   --verbose                 Show detailed debug logs
#   --processed               Use adata.X instead of adata.raw.X (default: prefer raw.X)
#   --assay NAME              Assay name (default: RNA)
#   --no-cleanup              Keep temporary files on failure
#   --skip-deps-check         Skip the pre-flight dependency check
#
# Examples:
#   bash h5ad2r.sh data/processed.h5ad
#   bash h5ad2r.sh data/obj.h5ad output.qs --force
#   bash h5ad2r.sh data/obj.h5ad output.Rdata
#   bash h5ad2r.sh data/obj.h5ad --method mtx --processed --verbose
#
# From Windows host (via WSL):
#   wsl bash D:/Projects/r2h5ad/h5ad2r.sh D:/data/obj.h5ad
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_DIR="${SCRIPT_DIR}/R"

# ── Defaults ──────────────────────────────────────────────────────────────────
FORCE=false
VERBOSE=false
METHOD="auto"
ASSAY="RNA"
NO_CLEANUP=false
SKIP_DEPS_CHECK=false
PROCESSED=false
INPUT_FILE=""
OUTPUT_FILE=""

# ── Functions ─────────────────────────────────────────────────────────────────

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

# Convert WSL path to Windows path
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

# Resolve output path from input + optional specified output
resolve_output() {
    local input="$1"
    local specified="$2"
    if [ -n "$specified" ]; then
        echo "$specified"
    else
        local dir="$(dirname "$input")"
        local base="$(basename "$input")"
        local stem="${base%.h5ad}"
        stem="${stem%.h5ad}"
        [ "$dir" = "." ] && echo "${stem}.rds" || echo "${dir}/${stem}.rds"
    fi
}

# Detect output file type by extension
detect_output_type() {
    local fname="$1"
    local lower
    lower="$(echo "$fname" | tr '[:upper:]' '[:lower:]')"
    case "$lower" in
        *.rdata|*.rda) echo "rdata" ;;
        *.qs)           echo "qs" ;;
        *)              echo "rds" ;;
    esac
}

# Validate input is .h5ad
check_input_h5ad() {
    local fname="$1"
    local lower
    lower="$(echo "$fname" | tr '[:upper:]' '[:lower:]')"
    if [[ ! "$lower" =~ \.h5ad$ ]]; then
        log ERROR "Input file must be .h5ad, got: $fname"
        exit 1
    fi
}

# Check an R package is loadable
r_has_pkg() {
    Rscript -e "cat(requireNamespace('${1}', quietly=TRUE))" 2>/dev/null | grep -q TRUE
}

# Check a Python package is importable
py_has_pkg() {
    python3 -c "import ${1}" 2>/dev/null
}

print_deps_help() {
    cat << 'DEPS_HELP'

==================== Dependency Setup ====================

h5ad2r requires R and Python with these packages:

  R:  Seurat, SeuratDisk, qs, jsonlite, Matrix
  Python: anndata, scanpy

Quick setup with conda (recommended):
  conda create -n r2h5ad -c conda-forge \
      r-seurat r-qs r-jsonlite r-matrix \
      python=3.10 anndata scanpy
  conda activate r2h5ad
  Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

Set env var R2H5AD_CONDA_ENV to point to an existing env:
  export R2H5AD_CONDA_ENV=my_env

===========================================================

DEPS_HELP
}

print_usage() {
    cat << 'EOF'
h5ad2r.sh — Convert AnnData .h5ad back to R-native formats

Usage: bash h5ad2r.sh <input.h5ad> [output_file] [options]

Output formats (auto-detected by extension):
  .rds                  Single Seurat object (via saveRDS)
  .qs                   Single Seurat object (via qs::qsave)
  .Rdata / .RData / .rda Workspace dump with multiple R objects

Options:
  --method seuratdisk|mtx   Force conversion method
  --force                   Overwrite existing output file
  --verbose                 Show detailed logs
  --processed               Use adata.X instead of adata.raw.X (default: prefer raw.X)
  --assay NAME              Assay name (default: RNA)
  --no-cleanup              Keep temporary files on failure
  --skip-deps-check         Skip the pre-flight dependency check

Examples:
  bash h5ad2r.sh data/processed.h5ad
  bash h5ad2r.sh data/obj.h5ad output.qs --force
  bash h5ad2r.sh data/obj.h5ad output.Rdata --processed
  bash h5ad2r.sh data/obj.h5ad --method mtx --verbose
EOF
}

# ── Conda activation ──────────────────────────────────────────────────────────
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

# ── Pre-flight checks ─────────────────────────────────────────────────────────
check_prerequisites() {
    try_activate_conda || true

    if ! have_cmd Rscript; then
        log ERROR "Rscript not found on PATH"
        log INFO  "Install R (>=4.0) or activate a conda environment that includes R."
        exit 1
    fi

    if ! have_cmd python3; then
        log ERROR "python3 not found on PATH"
        exit 1
    fi

    local missing_scripts=()
    for script in utils.R convert_h5ad_seuratdisk.R extract_h5ad.py assemble_seurat.R assemble_rdata.R; do
        [ -f "${R_DIR}/${script}" ] || missing_scripts+=("$script")
    done
    if [ ${#missing_scripts[@]} -gt 0 ]; then
        log ERROR "Missing helper scripts: ${missing_scripts[*]}"
        log ERROR "Expected in: ${R_DIR}/"
        exit 1
    fi
}

# ── Dependency check ──────────────────────────────────────────────────────────
check_deps() {
    if [ "$SKIP_DEPS_CHECK" = true ]; then
        log DEBUG "Skipping dependency check (--skip-deps-check)"
        return 0
    fi

    local missing_r=""
    local missing_py=""

    for pkg in Seurat SeuratDisk qs jsonlite Matrix; do
        r_has_pkg "$pkg" || missing_r="${missing_r} ${pkg}"
    done

    for pkg in anndata scanpy; do
        py_has_pkg "$pkg" || missing_py="${missing_py} ${pkg}"
    done

    if [ -n "$missing_r" ] || [ -n "$missing_py" ]; then
        log ERROR "Missing dependencies detected:"
        [ -n "$missing_r" ] && log ERROR "  R packages:${missing_r}"
        [ -n "$missing_py" ] && log ERROR "  Python packages:${missing_py}"
        print_deps_help
        log INFO  "Or use --skip-deps-check to proceed anyway (at your own risk)."
        exit 1
    fi

    log DEBUG "All dependencies present."
}

# ── Validate RDS/QS output ────────────────────────────────────────────────────
validate_r_output() {
    local r_path="$1"
    local r_ext="$2"

    if [ "$r_ext" = "rds" ]; then
        Rscript -e "
obj <- readRDS('${r_path}')
cls <- class(obj)[1]
cat(sprintf('Class: %s\n', cls))
if (inherits(obj, 'Seurat')) {
cat(sprintf('Cells: %d\n', ncol(obj)))
cat(sprintf('Features: %d\n', nrow(obj)))
cat(sprintf('Assays: %s\n', paste(names(obj@assays), collapse=', ')))
}" 2>&1
    elif [ "$r_ext" = "qs" ]; then
        Rscript -e "
obj <- qs::qread('${r_path}')
cls <- class(obj)[1]
cat(sprintf('Class: %s\n', cls))
if (inherits(obj, 'Seurat')) {
cat(sprintf('Cells: %d\n', ncol(obj)))
cat(sprintf('Features: %d\n', nrow(obj)))
cat(sprintf('Assays: %s\n', paste(names(obj@assays), collapse=', ')))
}" 2>&1
    fi
}

# ── Validate Rdata output ─────────────────────────────────────────────────────
validate_rdata_output() {
    local r_path="$1"
    Rscript -e "
load('${r_path}')
objs <- ls()
cat(sprintf('Objects loaded: %s\n', paste(objs, collapse=', ')))
for (nm in objs) {
    obj <- get(nm)
    if (is.matrix(obj) || inherits(obj, 'Matrix')) {
        cat(sprintf('  %s: %d x %d\n', nm, nrow(obj), ncol(obj)))
    } else if (is.data.frame(obj)) {
        cat(sprintf('  %s: %d rows x %d cols\n', nm, nrow(obj), ncol(obj)))
    } else {
        cat(sprintf('  %s: %s\n', nm, class(obj)[1]))
    }
}" 2>&1
}

# ── Argument Parsing ──────────────────────────────────────────────────────────

while [ $# -gt 0 ]; do
    case "$1" in
        --method)           shift; METHOD="$1" ;;
        --force)            FORCE=true ;;
        --verbose)          VERBOSE=true ;;
        --processed)       PROCESSED=true ;;
        --assay)            shift; ASSAY="$1" ;;
        --no-cleanup)       NO_CLEANUP=true ;;
        --skip-deps-check)  SKIP_DEPS_CHECK=true ;;
        --help|-h)          print_usage; exit 0 ;;
        -*)
            log ERROR "Unknown option: $1"
            print_usage
            exit 1 ;;
        *)
            if [ -z "$INPUT_FILE" ]; then INPUT_FILE="$1"
            elif [ -z "$OUTPUT_FILE" ]; then OUTPUT_FILE="$1"
            else log ERROR "Unexpected argument: $1"; print_usage; exit 1
            fi
            ;;
    esac
    shift
done

# ── Main ──────────────────────────────────────────────────────────────────────

[ -z "$INPUT_FILE" ] && { log ERROR "No input file specified"; print_usage; exit 1; }

# Convert Windows paths to WSL
INPUT_FILE=$(win_to_wsl "$INPUT_FILE")
OUTPUT_FILE=$(resolve_output "$INPUT_FILE" "$(win_to_wsl "${OUTPUT_FILE:-}")")

# Validate input
check_input_h5ad "$INPUT_FILE"

# Detect output type
OUTPUT_TYPE=$(detect_output_type "$OUTPUT_FILE")

# Validate method + output type combos
if [ "$OUTPUT_TYPE" = "rdata" ] && [ "$METHOD" != "auto" ] && [ "$METHOD" != "mtx" ]; then
    log ERROR "--method seuratdisk is not valid for .Rdata output; use --method mtx or omit"
    exit 1
fi

log INFO "======== h5ad2r: $(wsl_to_win "$INPUT_FILE") → R ========"
log DEBUG "Input file (WSL):   $INPUT_FILE"
log DEBUG "Output file (WSL):  $OUTPUT_FILE"
log DEBUG "Output type:        $OUTPUT_TYPE"
log DEBUG "Method:             $METHOD"
log DEBUG "Processed mode:      $PROCESSED"

# Validate input file exists
[ -f "$INPUT_FILE" ] || { log ERROR "Input file not found: $INPUT_FILE"; exit 1; }

# Pre-flight checks
check_prerequisites

# Dependency check
check_deps

# Output file check
if [ -f "$OUTPUT_FILE" ] && [ "$FORCE" != true ]; then
    log ERROR "Output file already exists: $OUTPUT_FILE"
    log INFO  "Use --force to overwrite"
    exit 1
fi
[ -f "$OUTPUT_FILE" ] && [ "$FORCE" = true ] && { log WARN "Overwriting existing file: $OUTPUT_FILE"; rm -f "$OUTPUT_FILE"; }

# Create temp directory
TMPDIR=$(mktemp -d -t h5ad2r.XXXXXXXXXX)
trap 'if [ "$NO_CLEANUP" != true ]; then rm -rf "$TMPDIR"; else log WARN "Temp files kept at: $TMPDIR"; fi' EXIT
log DEBUG "Temp dir: $TMPDIR"

# ── Conversion ────────────────────────────────────────────────────────────────

CONVERT_EXIT=0

if [ "$OUTPUT_TYPE" = "rdata" ] || [ "$METHOD" = "mtx" ]; then
    # ── MTX intermediate path (for .Rdata or forced --method mtx) ──────────────
    log INFO "[1/2] Extracting data from h5ad..."

    PY_SCRIPT="${R_DIR}/extract_h5ad.py"
    [ -f "$PY_SCRIPT" ] || { log ERROR "Missing Python script: $PY_SCRIPT"; exit 1; }

    PY_ARGS=("$INPUT_FILE" "$TMPDIR" "--assay" "$ASSAY")
    [ "$PROCESSED" = true ] && PY_ARGS+=("--processed")

    python3 "$PY_SCRIPT" "${PY_ARGS[@]}" 2>&1 || {
        log ERROR "h5ad extraction failed"
        exit 1
    }

    # Parse summary
    SUMMARY_PATH="${TMPDIR}/summary.json"
    if [ -f "$SUMMARY_PATH" ]; then
        N_CELLS=$(python3 -c "import json; d=json.load(open('$SUMMARY_PATH')); print(d['shape'][0])" 2>/dev/null || echo "?")
        N_GENES=$(python3 -c "import json; d=json.load(open('$SUMMARY_PATH')); print(d['shape'][1])" 2>/dev/null || echo "?")
        log INFO "Extracted: $N_CELLS cells × $N_GENES genes"
    fi

    # Route to assembler
    if [ "$OUTPUT_TYPE" = "rdata" ]; then
        log INFO "[2/2] Assembling .Rdata..."
        R_SCRIPT="${R_DIR}/assemble_rdata.R"
        [ -f "$R_SCRIPT" ] || { log ERROR "Missing R script: $R_SCRIPT"; exit 1; }
        export R2H5AD_VERBOSE=$VERBOSE
        R_OUTPUT=$(Rscript "$R_SCRIPT" "$TMPDIR" "$OUTPUT_FILE" 2>&1) || {
            log ERROR "Rdata assembly failed:"
            echo "$R_OUTPUT" >&2
            exit 1
        }
    else
        log INFO "[2/2] Assembling Seurat object..."
        R_SCRIPT="${R_DIR}/assemble_seurat.R"
        [ -f "$R_SCRIPT" ] || { log ERROR "Missing R script: $R_SCRIPT"; exit 1; }
        export R2H5AD_VERBOSE=$VERBOSE
        R_OUTPUT=$(Rscript "$R_SCRIPT" "$TMPDIR" "$OUTPUT_FILE" 2>&1) || {
            log ERROR "Seurat assembly failed:"
            echo "$R_OUTPUT" >&2
            exit 1
        }
    fi

    # Print R summary
    R_SUMMARY=$(echo "$R_OUTPUT" | grep -E '^\s*\{' | tail -1)
    if echo "$R_SUMMARY" | python3 -c "import sys,json; d=json.load(sys.stdin); sys.exit(0 if d.get('status')=='success' else 1)" 2>/dev/null; then
        N_CELLS=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('cells','?'))")
        N_FEAT=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('features','?'))")
        log INFO "Created: $N_CELLS cells × $N_FEAT features"
    fi

else
    # ── SeuratDisk reverse path (.rds / .qs) ───────────────────────────────────
    log INFO "Converting h5ad → Seurat via SeuratDisk..."

    R_SCRIPT="${R_DIR}/convert_h5ad_seuratdisk.R"
    [ -f "$R_SCRIPT" ] || { log ERROR "Missing R script: $R_SCRIPT"; exit 1; }

    export R2H5AD_VERBOSE=$VERBOSE
    R_OUTPUT=$(Rscript "$R_SCRIPT" "$INPUT_FILE" "$OUTPUT_FILE" 2>&1) || CONVERT_EXIT=$?

    if [ $CONVERT_EXIT -ne 0 ]; then
        log WARN "SeuratDisk method failed, falling back to MTX path..."
        log DEBUG "SeuratDisk error: $(echo "$R_OUTPUT" | tail -5)"
        rm -f "$OUTPUT_FILE"

        # MTX fallback
        log INFO "[1/2] Extracting data from h5ad (MTX fallback)..."
        PY_SCRIPT="${R_DIR}/extract_h5ad.py"
        PY_ARGS=("$INPUT_FILE" "$TMPDIR" "--assay" "$ASSAY")
        [ "$PROCESSED" = true ] && PY_ARGS+=("--processed")
        python3 "$PY_SCRIPT" "${PY_ARGS[@]}" 2>&1 || {
            log ERROR "h5ad extraction failed"
            exit 1
        }

        log INFO "[2/2] Assembling Seurat object (MTX fallback)..."
        R_SCRIPT="${R_DIR}/assemble_seurat.R"
        export R2H5AD_VERBOSE=$VERBOSE
        R_OUTPUT=$(Rscript "$R_SCRIPT" "$TMPDIR" "$OUTPUT_FILE" 2>&1) || {
            log ERROR "Seurat assembly failed:"
            echo "$R_OUTPUT" >&2
            exit 1
        }
        CONVERT_EXIT=0
    fi

    if [ $CONVERT_EXIT -ne 0 ]; then
        log ERROR "Conversion failed (exit code: $CONVERT_EXIT)"
        exit $CONVERT_EXIT
    fi

    # Print R summary
    R_SUMMARY=$(echo "$R_OUTPUT" | grep -E '^\s*\{' | tail -1)
    if echo "$R_SUMMARY" | python3 -c "import sys,json; d=json.load(sys.stdin); sys.exit(0 if d.get('status')=='success' else 1)" 2>/dev/null; then
        METHOD_USED=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('method','?'))")
        N_CELLS=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('cells','?'))")
        N_FEAT=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('features','?'))")
        log INFO "Method: $METHOD_USED, $N_CELLS cells × $N_FEAT features"
    fi
fi

[ -f "$OUTPUT_FILE" ] || { log ERROR "Output file was not created: $OUTPUT_FILE"; exit 1; }

# ── Validate ──────────────────────────────────────────────────────────────────

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo 0)
[ "$OUTPUT_SIZE" -gt 0 ] || { log ERROR "Output file is empty: $OUTPUT_FILE"; exit 1; }
log INFO "Output size: $(printf "%'d" "$OUTPUT_SIZE") bytes ($(python3 -c "print(round($OUTPUT_SIZE/1e6,1))") MB)"

log INFO "Validating output..."
if [ "$OUTPUT_TYPE" = "rdata" ]; then
    VALIDATION=$(validate_rdata_output "$OUTPUT_FILE" 2>&1) && VALIDATE_EXIT=0 || VALIDATE_EXIT=1
else
    VALIDATION=$(validate_r_output "$OUTPUT_FILE" "$OUTPUT_TYPE" 2>&1) && VALIDATE_EXIT=0 || VALIDATE_EXIT=1
fi

if [ $VALIDATE_EXIT -eq 0 ]; then
    log INFO "Validation passed:"
    echo "$VALIDATION" | while IFS= read -r line; do log INFO "  $line"; done
else
    log WARN "Validation returned warnings (file may still be usable):"
    echo "$VALIDATION" >&2
fi

log INFO "======== Conversion complete ========"
log INFO "Output: $(wsl_to_win "$OUTPUT_FILE")"
echo "$OUTPUT_FILE"
