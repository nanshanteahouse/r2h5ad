#!/usr/bin/env bash
# =============================================================================
# r2h5ad.sh - RDS/QS to h5ad Conversion Tool
# =============================================================================
# Converts single-cell data files in RDS (.rds) or QS (.qs) format to h5ad
# for use in the RNA/ATAC pipeline loaders.
#
# Usage:
#   bash r2h5ad.sh <input_file> [output_file] [options]
#
# Options:
#   --method seuratdisk|mtx   Force specific conversion method
#   --force                   Overwrite existing output file
#   --verbose                 Show detailed logs
#   --assay NAME              Specify assay to use (default: RNA)
#   --no-cleanup              Keep temporary files on failure
#   --skip-deps-check         Skip the pre-flight dependency check
#
# Environment variables:
#   R2H5AD_CONDA_ENV          Name of conda env to use (auto-detected otherwise)
#
# Examples:
#   bash r2h5ad.sh data/seurat_obj.rds
#   bash r2h5ad.sh data/seurat_obj.rds output/processed.h5ad --force
#   bash r2h5ad.sh data/sce_obj.rds --method mtx
#   bash r2h5ad.sh data/obj.qs --verbose
#
# From Windows host:
#   wsl bash D:/Projects/r2h5ad/r2h5ad.sh D:/data/obj.rds
# =============================================================================

set -euo pipefail

# --- Variables ---------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_DIR="${SCRIPT_DIR}/R"

# Defaults
FORCE=false
VERBOSE=false
METHOD="auto"
ASSAY="RNA"
NO_CLEANUP=false
SKIP_DEPS_CHECK=false
INPUT_FILE=""
OUTPUT_FILE=""

# --- Functions ---------------------------------------------------------------

log() {
    local level="$1"; shift
    local msg="$*"
    case "$level" in
        INFO)  echo "[INFO] $msg" ;;
        WARN)  echo -e "\033[33m[WARN]\033[0m $msg" >&2 ;;
        ERROR) echo -e "\033[31m[ERROR]\033[0m $msg" >&2 ;;
        DEBUG) $VERBOSE && echo "[DEBUG] $msg" ;;
    esac
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

usage() {
    cat << 'EOF'
r2h5ad.sh - Convert RDS/QS single-cell files to h5ad format

Usage: bash r2h5ad.sh <input_file> [output_file] [options]

Options:
  --method seuratdisk|mtx   Force specific conversion method
  --force                   Overwrite existing output file
  --verbose                 Show detailed logs
  --assay NAME              Specify assay to use (default: RNA)
  --no-cleanup              Keep temporary files on failure
  --skip-deps-check         Skip the pre-flight dependency check

Examples:
  bash r2h5ad.sh data/seurat_obj.rds
  bash r2h5ad.sh data/seurat_obj.rds output.h5ad --force
  bash r2h5ad.sh data/obj.qs --verbose
EOF
    exit 1
}

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
        if [ "$dir" = "." ]; then echo "${stem}.h5ad"
        else echo "${dir}/${stem}.h5ad"; fi
    fi
}

# Check a binary is on PATH
have_cmd() { command -v "$1" &>/dev/null; }

# Check an R package is loadable
r_has_pkg() {
    Rscript -e "cat(requireNamespace('${1}', quietly=TRUE))" 2>/dev/null | grep -q TRUE
}

# Check a Python package is importable
py_has_pkg() {
    python3 -c "import ${1}" 2>/dev/null
}

# --- Auto-detect conda (try common paths) ------------------------------------
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

# --- Pre-flight checks ------------------------------------------------------
check_prerequisites() {
    # 1. Try conda activation (best-effort)
    try_activate_conda || true

    # 2. Check Rscript
    if ! have_cmd Rscript; then
        log ERROR "Rscript not found on PATH"
        log INFO  "Install R (>=4.0) or activate a conda environment that includes R."
        exit 1
    fi

    # 3. Check python3
    if ! have_cmd python3; then
        log ERROR "python3 not found on PATH"
        exit 1
    fi

    # 4. Check helper scripts exist
    local missing_scripts=()
    for script in utils.R detect_format.R convert_seuratdisk.R convert_mtx.R; do
        [ -f "${R_DIR}/${script}" ] || missing_scripts+=("$script")
    done
    if [ ${#missing_scripts[@]} -gt 0 ]; then
        log ERROR "Missing R scripts: ${missing_scripts[*]}"
        log ERROR "Expected in: ${R_DIR}/"
        exit 1
    fi
}

# --- Dependency check --------------------------------------------------------
check_deps() {
    if [ "$SKIP_DEPS_CHECK" = true ]; then
        log DEBUG "Skipping dependency check (--skip-deps-check)"
        return 0
    fi

    local missing_r=""
    local missing_py=""

    # Required R packages
    for pkg in Seurat SeuratDisk qs jsonlite Matrix hdf5r; do
        r_has_pkg "$pkg" || missing_r="${missing_r} ${pkg}"
    done

    # Required Python packages
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

# --- Validate h5ad output ---------------------------------------------------
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

# --- Argument Parsing --------------------------------------------------------

while [ $# -gt 0 ]; do
    case "$1" in
        --method)       shift; METHOD="$1" ;;
        --force)        FORCE=true ;;
        --verbose)      VERBOSE=true ;;
        --assay)        shift; ASSAY="$1" ;;
        --no-cleanup)   NO_CLEANUP=true ;;
        --skip-deps-check) SKIP_DEPS_CHECK=true ;;
        --help|-h)      usage ;;
        -*)             log ERROR "Unknown option: $1"; usage ;;
        *)
            if [ -z "$INPUT_FILE" ]; then INPUT_FILE="$1"
            elif [ -z "$OUTPUT_FILE" ]; then OUTPUT_FILE="$1"
            else log ERROR "Unexpected argument: $1"; usage
            fi
            ;;
    esac
    shift
done

# --- Main --------------------------------------------------------------------

[ -z "$INPUT_FILE" ] && { log ERROR "No input file specified"; usage; }

log INFO "======== r2h5ad: RDS/QS to h5ad Conversion ========"

# Convert Windows paths to WSL
INPUT_FILE=$(win_to_wsl "$INPUT_FILE")
OUTPUT_FILE=$(resolve_output "$INPUT_FILE" "$(win_to_wsl "${OUTPUT_FILE:-}")")

log DEBUG "Input file (WSL):  $INPUT_FILE"
log DEBUG "Output file (WSL): $OUTPUT_FILE"

# Validate input
[ -f "$INPUT_FILE" ] || { log ERROR "Input file not found: $INPUT_FILE"; exit 1; }

if [ -f "$OUTPUT_FILE" ] && [ "$FORCE" != true ]; then
    log ERROR "Output file already exists: $OUTPUT_FILE"
    log INFO  "Use --force to overwrite"
    exit 1
fi

# Pre-flight (conda, Rscript, python3, helper scripts)
check_prerequisites

# Dependency check (can be skipped)
check_deps

[ -f "$OUTPUT_FILE" ] && [ "$FORCE" = true ] && {
    log WARN "Overwriting existing file: $OUTPUT_FILE"
    rm -f "$OUTPUT_FILE"
}

log INFO "Input:  $(wsl_to_win "$INPUT_FILE")"
log INFO "Output: $(wsl_to_win "$OUTPUT_FILE")"
log INFO "Method: ${METHOD}"

# Create temp directory
TMPDIR=$(mktemp -d -t r2h5ad.XXXXXXXXXX)
trap 'if [ "$NO_CLEANUP" != true ]; then rm -rf "$TMPDIR"; else log WARN "Temp files kept at: $TMPDIR"; fi' EXIT

log DEBUG "Temp dir: $TMPDIR"

# --- Step 1: Detect format ---------------------------------------------------
log INFO "Detecting file format..."

DETECT_RESULT=$(
    Rscript "${R_DIR}/detect_format.R" "$INPUT_FILE" 2>"$TMPDIR/detect_stderr.txt"
) || {
    log ERROR "Format detection failed:"
    cat "$TMPDIR/detect_stderr.txt" >&2
    exit 1
}

log DEBUG "Detection result: $DETECT_RESULT"

# Parse JSON
FILE_TYPE=$(echo "$DETECT_RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('file_type','unknown'))" 2>/dev/null || echo "unknown")
OBJ_CLASS=$(echo "$DETECT_RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); print(d.get('object_class','unknown'))" 2>/dev/null || echo "unknown")
OBJ_DIMS=$(echo "$DETECT_RESULT" | python3 -c "import sys,json; d=json.load(sys.stdin); dims=d.get('dims',{}); print(f\"{dims.get('cells',0)} cells x {dims.get('features',0)} features\")" 2>/dev/null || echo "?")

log INFO "File type:    $FILE_TYPE"
log INFO "Object class: $OBJ_CLASS"
log INFO "Dimensions:   $OBJ_DIMS"

# --- Step 2: Route to converter ----------------------------------------------
CONVERT_EXIT=0

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

# --- Check conversion result ------------------------------------------------
if [ $CONVERT_EXIT -ne 0 ]; then
    log ERROR "Conversion failed (exit code: $CONVERT_EXIT)"
    [ -f "$TMPDIR/convert_stderr.txt" ] && { log ERROR "Stderr:"; cat "$TMPDIR/convert_stderr.txt" >&2; }
    exit $CONVERT_EXIT
fi

[ -f "$OUTPUT_FILE" ] || { log ERROR "Output file was not created: $OUTPUT_FILE"; exit 1; }

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_FILE" 2>/dev/null || echo 0)
[ "$OUTPUT_SIZE" -gt 0 ] || { log ERROR "Output file is empty: $OUTPUT_FILE"; exit 1; }

log INFO "Output size: $(printf "%'d" "$OUTPUT_SIZE") bytes"

# Validate h5ad
log INFO "Validating output h5ad..."
VALIDATION=$(validate_h5ad "$OUTPUT_FILE" 2>&1) && VALIDATE_EXIT=0 || VALIDATE_EXIT=1
if [ $VALIDATE_EXIT -eq 0 ]; then
    log INFO "Validation passed:"
    echo "$VALIDATION" | while IFS= read -r line; do log INFO "  $line"; done
else
    log WARN "Validation returned warnings (file may still be usable):"
    echo "$VALIDATION" >&2
fi

# --- Done --------------------------------------------------------------------
log INFO "======== Conversion complete ========"
log INFO "Output file: $(wsl_to_win "$OUTPUT_FILE")"
echo "$OUTPUT_FILE"
