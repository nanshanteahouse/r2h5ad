#!/usr/bin/env bash
# =============================================================================
# r2h5ad.sh — Convert R single-cell files (.rds, .qs, .Rdata) to h5ad
# =============================================================================
# Converts single-cell data from R-native formats to AnnData (.h5ad).
#
# Supported input formats:
#   .rds / .qs             — Single serialized R objects (Seurat, SCE, etc.)
#                            Uses detect_format.R → SeuratDisk or MTX export
#   .rds.gz / .qs.gz       — Gzip-compressed serialized objects
#   .Rdata / .RData / .rda — Legacy workspace dumps (multiple named objects)
#                            Uses convert_rdata.R → assemble_h5ad.py
#
# Usage:
#   bash r2h5ad.sh <input> [output] [options]
#
# Options:
#   --method seuratdisk|mtx   Force conversion method (serialized objects only)
#   --force                   Overwrite existing output file
#   --verbose                 Show detailed debug logs
#   --assay NAME              Specify assay to use (default: RNA, serialized objects)
#   --no-cleanup              Keep temporary files on failure
#   --skip-deps-check         Skip the pre-flight dependency check
#   --list                    List objects in .Rdata and exit (Rdata only)
#   --count-object NAME       Name of count matrix object (Rdata only, auto-detected)
#   --filter-object NAME      Object whose rownames define cells to keep (Rdata only)
#   --pca-object NAME         PCA coordinates object name (Rdata only)
#   --tsne-object NAME        t-SNE coordinates object name (Rdata only)
#   --umap-object NAME        UMAP coordinates object name (Rdata only)
#
# Environment variables:
#   R2H5AD_CONDA_ENV          Name of conda env to use (auto-detected otherwise)
#
# Examples:
#   bash r2h5ad.sh seurat_obj.rds
#   bash r2h5ad.sh obj.Rdata output.h5ad --force
#   bash r2h5ad.sh obj.Rdata --list
#   bash r2h5ad.sh sce.rds --method mtx
#   bash r2h5ad.sh obj.qs --verbose
#
# From Windows host (via WSL):
#   wsl bash D:/Projects/r2h5ad/r2h5ad.sh D:/data/obj.rds
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
LIST_ONLY=false
INPUT_FILE=""
OUTPUT_FILE=""
PASSTHROUGH_ARGS=()

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

# Convert Windows path to WSL path (e.g. D:\Projects\data.rds -> /mnt/d/Projects/data.rds)
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

# Convert WSL path to Windows path (e.g. /mnt/d/Projects/data.rds -> D:/Projects/data.rds)
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
        local stem="${base%.*}"
        stem="${stem%.rds}"
        stem="${stem%.qs}"
        stem="${stem%.Rdata}"
        stem="${stem%.rdata}"
        stem="${stem%.RData}"
        stem="${stem%.rda}"
        [ "$dir" = "." ] && echo "${stem}.h5ad" || echo "${dir}/${stem}.h5ad"
    fi
}

# Detect file type by extension
detect_file_type() {
    local fname="$1"
    local lower
    lower="$(echo "$fname" | tr '[:upper:]' '[:lower:]')"
    case "$lower" in
        *.rdata|*.rda) echo "rdata" ;;
        *.rds|*.qs|*.rds.gz|*.qs.gz) echo "serialized" ;;
        *) echo "unknown" ;;
    esac
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

r2h5ad requires R and Python with these packages:

  R:  Seurat, SeuratDisk, qs, jsonlite, Matrix, hdf5r
  Python: anndata, scanpy

Quick setup with conda (recommended):
  conda create -n r2h5ad -c conda-forge \
      r-seurat r-qs r-hdf5r r-jsonlite r-matrix \
      python=3.10 anndata scanpy
  conda activate r2h5ad
  Rscript -e 'remotes::install_github("mojaveazure/seurat-disk")'

Alternatively, set env var R2H5AD_CONDA_ENV to point to an existing env:
  export R2H5AD_CONDA_ENV=my_env

Without conda (manual install):
  R:   install.packages(c("Seurat","qs","jsonlite","Matrix"))
       remotes::install_github("mojaveazure/seurat-disk")
  Python: pip install anndata scanpy

===========================================================

DEPS_HELP
}

print_usage() {
    cat << 'EOF'
r2h5ad.sh — Convert R single-cell files (.rds, .qs, .Rdata) to h5ad

Usage: bash r2h5ad.sh <input_file> [output_file] [options]

Options:
  --method seuratdisk|mtx   Force conversion method (serialized objects only)
  --force                   Overwrite existing output file
  --verbose                 Show detailed logs
  --assay NAME              Specify assay to use (serialized objects, default: RNA)
  --no-cleanup              Keep temporary files on failure
  --skip-deps-check         Skip the pre-flight dependency check
  --list                    List objects in .Rdata and exit (Rdata only)
  --count-object NAME       Name of count matrix (Rdata only, auto-detected)
  --filter-object NAME      Object whose rownames define cells to keep (Rdata only)
  --pca-object NAME         PCA coordinates object name (Rdata only)
  --tsne-object NAME        t-SNE coordinates object name (Rdata only)
  --umap-object NAME        UMAP coordinates object name (Rdata only)

Input formats:
  .rds / .qs / .rds.gz / .qs.gz   Serialized R objects (detected automatically)
  .Rdata / .RData / .rda          Legacy workspace dumps (detected automatically)

Examples:
  bash r2h5ad.sh seurat_obj.rds
  bash r2h5ad.sh obj.Rdata output.h5ad --force
  bash r2h5ad.sh obj.Rdata --list
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

# ── Pre-flight checks ─────────────────────────────────────────────────────────
check_prerequisites() {
    local file_type="$1"

    # Conda activation (best-effort)
    try_activate_conda || true

    # Rscript
    if ! have_cmd Rscript; then
        log ERROR "Rscript not found on PATH"
        log INFO  "Install R (>=4.0) or activate a conda environment that includes R."
        exit 1
    fi

    # Python3
    if ! have_cmd python3; then
        log ERROR "python3 not found on PATH"
        exit 1
    fi

    # Helper scripts for the detected file type
    local missing_scripts=()
    if [ "$file_type" = "rdata" ]; then
        for script in convert_rdata.R assemble_h5ad.py; do
            [ -f "${R_DIR}/${script}" ] || missing_scripts+=("$script")
        done
    elif [ "$file_type" = "serialized" ]; then
        for script in utils.R detect_format.R convert_seuratdisk.R convert_mtx.R; do
            [ -f "${R_DIR}/${script}" ] || missing_scripts+=("$script")
        done
    fi
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

    for pkg in Seurat SeuratDisk qs jsonlite Matrix hdf5r; do
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

# ── Validate h5ad output ──────────────────────────────────────────────────────
validate_h5ad() {
    local h5ad_path="$1"
    python3 -c "
import anndata
adata = anndata.read_h5ad('${h5ad_path}')
print('Shape:', adata.shape)
obs = list(adata.obs.columns) if adata.obs.shape[1] > 0 else 'none'
var = list(adata.var.columns) if adata.var.shape[1] > 0 else 'none'
print('Obs keys:', obs)
print('Var keys:', var)
" 2>&1
}

# ── List Rdata objects ────────────────────────────────────────────────────────
list_rdata_objects() {
    local input="$1"
    log INFO "Listing objects in .Rdata..."
    Rscript -e "
load('${input}')
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
}

# ── Argument Parsing ──────────────────────────────────────────────────────────

while [ $# -gt 0 ]; do
    case "$1" in
        --method)           shift; METHOD="$1" ;;
        --force)            FORCE=true ;;
        --verbose)          VERBOSE=true ;;
        --assay)            shift; ASSAY="$1" ;;
        --no-cleanup)       NO_CLEANUP=true ;;
        --skip-deps-check)  SKIP_DEPS_CHECK=true ;;
        --list)             LIST_ONLY=true ;;
        --count-object)     shift; PASSTHROUGH_ARGS+=(count-object="$1") ;;
        --filter-object)    shift; PASSTHROUGH_ARGS+=(filter-object="$1") ;;
        --pca-object)       shift; PASSTHROUGH_ARGS+=(pca-object="$1") ;;
        --tsne-object)      shift; PASSTHROUGH_ARGS+=(tsne-object="$1") ;;
        --umap-object)      shift; PASSTHROUGH_ARGS+=(umap-object="$1") ;;
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

# Detect file type
FILE_TYPE=$(detect_file_type "$INPUT_FILE")
[ "$FILE_TYPE" = "unknown" ] && { log ERROR "Unsupported file type: $INPUT_FILE"; print_usage; exit 1; }

# Validate option compatibility with file type
if [ "$FILE_TYPE" != "rdata" ] && [ "$LIST_ONLY" = true ]; then
    log ERROR "--list is only valid for .Rdata files"
    exit 1
fi
if [ "$FILE_TYPE" != "rdata" ] && [ ${#PASSTHROUGH_ARGS[@]} -gt 0 ]; then
    log ERROR "--count-object, --filter-object, --pca-object, --tsne-object, --umap-object are only valid for .Rdata files"
    exit 1
fi
if [ "$FILE_TYPE" = "rdata" ]; then
    if [ "$METHOD" != "auto" ]; then
        log ERROR "--method is not supported for .Rdata files (only MTX export is available)"
        exit 1
    fi
    if [ "$ASSAY" != "RNA" ]; then
        log ERROR "--assay is not supported for .Rdata files"
        exit 1
    fi
fi

log INFO "======== r2h5ad: $(wsl_to_win "$INPUT_FILE") → h5ad ========"
log DEBUG "Input file (WSL):  $INPUT_FILE"
log DEBUG "Output file (WSL): $OUTPUT_FILE"
log DEBUG "File type: $FILE_TYPE"

# Validate input file exists
[ -f "$INPUT_FILE" ] || { log ERROR "Input file not found: $INPUT_FILE"; exit 1; }

# Handle --list (lightweight — does not require full dependency check)
if [ "$LIST_ONLY" = true ]; then
    try_activate_conda || true
    have_cmd Rscript || { log ERROR "Rscript not found"; exit 1; }
    list_rdata_objects "$INPUT_FILE"
    exit 0
fi

# Pre-flight checks
check_prerequisites "$FILE_TYPE"

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
TMPDIR=$(mktemp -d -t r2h5ad.XXXXXXXXXX)
trap 'if [ "$NO_CLEANUP" != true ]; then rm -rf "$TMPDIR"; else log WARN "Temp files kept at: $TMPDIR"; fi' EXIT
log DEBUG "Temp dir: $TMPDIR"

# ── Conversion ────────────────────────────────────────────────────────────────

CONVERT_EXIT=0

if [ "$FILE_TYPE" = "rdata" ]; then
    # ── .Rdata path: convert_rdata.R → assemble_h5ad.py ──────────────────────
    log INFO "[1/2] Extracting data from .Rdata..."

    R_SCRIPT="${R_DIR}/convert_rdata.R"
    [ -f "$R_SCRIPT" ] || { log ERROR "Missing R script: $R_SCRIPT"; exit 1; }

    export R2H5AD_VERBOSE=$VERBOSE
    R_ARGS=("$INPUT_FILE" "$TMPDIR")
    for arg in "${PASSTHROUGH_ARGS[@]}"; do
        R_ARGS+=("--$arg")
    done
    [ "$VERBOSE" = true ] && R_ARGS+=("--verbose=true")

    R_OUTPUT=$(Rscript "$R_SCRIPT" "${R_ARGS[@]}" 2>&1) || {
        log ERROR "R conversion failed:"
        echo "$R_OUTPUT" >&2
        exit 1
    }

    # Parse JSON summary (last line starting with '{')
    R_SUMMARY=$(echo "$R_OUTPUT" | grep -E '^\s*\{' | tail -1)
    if echo "$R_SUMMARY" | python3 -c "import sys,json; d=json.load(sys.stdin); sys.exit(0 if d.get('status')=='success' else 1)" 2>/dev/null; then
        COUNT_OBJ=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('count_object','?'))")
        N_GENES=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('n_genes','?'))")
        N_CELLS=$(echo "$R_SUMMARY" | python3 -c "import sys,json; print(json.load(sys.stdin).get('n_cells','?'))")
        log INFO "Count matrix: '$COUNT_OBJ' — $N_GENES genes x $N_CELLS cells"
    else
        log WARN "Could not parse R summary. Raw output:"
        echo "$R_OUTPUT" | grep -E '^\s*\{' | head -1
    fi

    # Assemble h5ad via Python
    log INFO "[2/2] Assembling h5ad..."
    PY_SCRIPT="${R_DIR}/assemble_h5ad.py"
    [ -f "$PY_SCRIPT" ] || { log ERROR "Missing Python script: $PY_SCRIPT"; exit 1; }
    python3 "$PY_SCRIPT" "$TMPDIR" "$OUTPUT_FILE" 2>&1 || {
        log ERROR "h5ad assembly failed"
        exit 1
    }

else
    # ── RDS/QS path: detect_format.R → SeuratDisk or MTX ────────────────────
    log INFO "Detecting file format..."

    DETECT_RESULT=$(Rscript "${R_DIR}/detect_format.R" "$INPUT_FILE" 2>"$TMPDIR/detect_stderr.txt") || {
        log ERROR "Format detection failed:"
        cat "$TMPDIR/detect_stderr.txt" >&2
        exit 1
    }
    log DEBUG "Detection result: $DETECT_RESULT"

    # Parse JSON from detection
    OBJ_CLASS=$(echo "$DETECT_RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('object_class','unknown'))" 2>/dev/null || echo "unknown")
    OBJ_DIMS=$(echo "$DETECT_RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); dims=d.get('dims',{}); print(f\"{dims.get('cells',0)} cells x {dims.get('features',0)} features\")" 2>/dev/null || echo "?")

    log INFO "Object class: $OBJ_CLASS"
    log INFO "Dimensions:   $OBJ_DIMS"

    # Route to converter
    if [ "$METHOD" = "mtx" ]; then
        log INFO "Using MTX export method (forced)"
        Rscript "${R_DIR}/convert_mtx.R" "$INPUT_FILE" "$OUTPUT_FILE" "$ASSAY" \
            > "$TMPDIR/convert_stdout.txt" 2>"$TMPDIR/convert_stderr.txt" || CONVERT_EXIT=$?

    elif [ "$METHOD" = "seuratdisk" ]; then
        log INFO "Using SeuratDisk method (forced)"
        Rscript "${R_DIR}/convert_seuratdisk.R" "$INPUT_FILE" "$OUTPUT_FILE" \
            > "$TMPDIR/convert_stdout.txt" 2>"$TMPDIR/convert_stderr.txt" || CONVERT_EXIT=$?

    elif [[ "$OBJ_CLASS" =~ (Seurat|seurat) ]]; then
        log INFO "Detected Seurat object, trying SeuratDisk method..."
        Rscript "${R_DIR}/convert_seuratdisk.R" "$INPUT_FILE" "$OUTPUT_FILE" \
            > "$TMPDIR/convert_stdout.txt" 2>"$TMPDIR/convert_stderr.txt" || CONVERT_EXIT=$?

        if [ $CONVERT_EXIT -ne 0 ]; then
            log WARN "SeuratDisk method failed, falling back to MTX export..."
            log DEBUG "SeuratDisk error: $(cat "$TMPDIR/convert_stderr.txt")"
            rm -f "$OUTPUT_FILE"
            Rscript "${R_DIR}/convert_mtx.R" "$INPUT_FILE" "$OUTPUT_FILE" "$ASSAY" \
                > "$TMPDIR/convert_stdout.txt" 2>"$TMPDIR/convert_stderr.txt" || CONVERT_EXIT=$?
        fi
    else
        log INFO "Non-Seurat object, using MTX export method..."
        Rscript "${R_DIR}/convert_mtx.R" "$INPUT_FILE" "$OUTPUT_FILE" "$ASSAY" \
            > "$TMPDIR/convert_stdout.txt" 2>"$TMPDIR/convert_stderr.txt" || CONVERT_EXIT=$?
    fi

    if [ $CONVERT_EXIT -ne 0 ]; then
        log ERROR "Conversion failed (exit code: $CONVERT_EXIT)"
        [ -f "$TMPDIR/convert_stderr.txt" ] && { log ERROR "Stderr:"; cat "$TMPDIR/convert_stderr.txt" >&2; }
        exit $CONVERT_EXIT
    fi

    [ -f "$OUTPUT_FILE" ] || { log ERROR "Output file was not created: $OUTPUT_FILE"; exit 1; }
fi

# ── Validate ──────────────────────────────────────────────────────────────────

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo 0)
[ "$OUTPUT_SIZE" -gt 0 ] || { log ERROR "Output file is empty: $OUTPUT_FILE"; exit 1; }
log INFO "Output size: $(printf "%'d" "$OUTPUT_SIZE") bytes ($(python3 -c "print(round($OUTPUT_SIZE/1e6,1))") MB)"

log INFO "Validating output h5ad..."
VALIDATION=$(validate_h5ad "$OUTPUT_FILE" 2>&1) && VALIDATE_EXIT=0 || VALIDATE_EXIT=1
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
