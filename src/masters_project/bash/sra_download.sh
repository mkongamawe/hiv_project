#!/usr/bin/env bash
 
# Exit immediately if a command exits with non-zero status
set -euo pipefail
 
# --- Configuration ---
# Default configuration parameters (used $HOME instead of ~ to ensure expansion inside quotes)
ACCESSION_LIST=${1:-"$HOME/Downloads/SRR_Acc_List.txt"}
OUTPUT_DIR=${2:-"./fastq_output"}
THREADS=${THREADS:-4}
 
echo "=================================================="
echo "SRA to FASTQ Pipeline Initialized"
echo "Accession List: $ACCESSION_LIST"
echo "Output Directory: $OUTPUT_DIR"
echo "Threads Allocated: $THREADS"
echo "=================================================="
 
# 1. Validate inputs and environment
if [ ! -f "$ACCESSION_LIST" ]; then
    echo "ERROR: Accession list file '$ACCESSION_LIST' not found."
    echo "Usage: $0 [path/to/SraAccList.txt] [output_directory]"
    exit 1
fi
 
# Ensure required commands are available in the active conda environment
MISSING_TOOLS=()
for cmd in prefetch fasterq-dump pigz; do
    if ! command -v "$cmd" &> /dev/null; then
        # Map the command name to its respective Conda package name
        if [ "$cmd" = "prefetch" ] || [ "$cmd" = "fasterq-dump" ]; then
            MISSING_TOOLS+=("bioconda::sra-tools")
        elif [ "$cmd" = "pigz" ]; then
            MISSING_TOOLS+=("conda-forge::pigz")
        fi
    fi
done
 
# If there are missing tools, prompt the user for installation
if [ ${#MISSING_TOOLS[@]} -gt 0 ]; then
    # Remove duplicate package names
    UNIQUE_TOOLS=($(echo "${MISSING_TOOLS[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    echo "WARNING: The following required packages are missing from your environment:"
    printf ' - %s\n' "${UNIQUE_TOOLS[@]}"
    echo "--------------------------------------------------"
    read -p "Would you like to install them now? (y/N): " choice
    echo "--------------------------------------------------"
    case "$choice" in
        [yY][eE][sS]|[yY])
            CONDA_BIN=$(command -v mamba || command -v conda || echo "")
            if [ -z "$CONDA_BIN" ]; then
                echo "ERROR: Neither 'conda' nor 'mamba' was found in your PATH."
                echo "Please install the missing packages manually."
                exit 1
            fi
            echo "Installing missing packages via $CONDA_BIN..."
            "$CONDA_BIN" install -y "${UNIQUE_TOOLS[@]}"
            echo "Installation complete. Resuming script..."
            ;;
        *)
            echo "ERROR: Missing required tools. Script aborted."
            exit 1
            ;;
    esac
fi
 
# Create necessary directories
# Note: Using OUTPUT_DIR as the parent for temp files avoids filling up a small /tmp RAM disk
mkdir -p "$OUTPUT_DIR"
TMP_DIR=$(mktemp -d -p "$OUTPUT_DIR" .sra_scratch_XXXXXX)
 
# Ensure temporary directories are cleaned up on script exit or interruption
cleanup() {
    echo "Cleaning up temporary scratch directory ($TMP_DIR)..."
    rm -rf "$TMP_DIR"
}
trap cleanup EXIT
 
# 2. Process each accession number sequentially
echo -e "\nProcessing SRA accessions..."
while IFS= read -r srr || [[ -n "$srr" ]]; do
    # Skip empty lines or commented lines
    [[ -z "$srr" || "$srr" =~ ^# ]] && continue
    # Trim whitespaces/carriage returns (common in Windows-edited text files)
    srr=$(echo "$srr" | xargs)
    echo "--------------------------------------------------"
    echo "Checking Accession: $srr"
    echo "--------------------------------------------------"
    # RESUME CHECK: Look for existing compressed .fastq.gz files for this SRR
    if compgen -G "$OUTPUT_DIR/${srr}*.fastq.gz" > /dev/null; then
        echo "--> Found existing compressed outputs for $srr. Skipping!"
        continue
    fi
    # Step 2a: Prefetch individual SRA file into the local scratch directory
    echo "Prefetching $srr..."
    if ! prefetch "$srr" \
        --output-directory "$TMP_DIR" \
        --max-size 100G \
        --progress \
        --resume yes \
        --verify yes; then
        echo "WARNING: prefetch failed for $srr. Skipping to next."
        continue
    fi
    # Locate the exact prefetched .sra file so fasterq-dump doesn't re-download it from NCBI
    if [ -f "$TMP_DIR/$srr/$srr.sra" ]; then
        SRA_TARGET="$TMP_DIR/$srr/$srr.sra"
    elif [ -f "$TMP_DIR/$srr/$srr.sralite" ]; then
        SRA_TARGET="$TMP_DIR/$srr/$srr.sralite"
    else
        SRA_TARGET="$TMP_DIR/$srr"
    fi
 
    # Step 2b: Run fasterq-dump directly against the local file path
    echo "Extracting FASTQ files from local archive: $SRA_TARGET..."
    if fasterq-dump "$SRA_TARGET" \
        --outdir "$OUTPUT_DIR" \
        --temp "$TMP_DIR" \
        --split-files \
        --skip-technical \
        --threads "$THREADS" \
        --progress; then
        echo "Successfully extracted FASTQ files for $srr."
        # Parallel compression using pigz for any new fastq files belonging to this SRR
        echo "Compressing files for $srr..."
        pigz -p "$THREADS" "$OUTPUT_DIR"/${srr}*.fastq
        # Clean up the local .sra folder immediately to free up disk space for the next run
        echo "Cleaning up local .sra archive for $srr..."
        rm -rf "$TMP_DIR/${srr}" 2>/dev/null || true
    else
        echo "WARNING: extraction failed for $srr. Cleaning up partial output."
        rm -f "$OUTPUT_DIR"/${srr}*.fastq 2>/dev/null || true
        rm -rf "$TMP_DIR/${srr}" 2>/dev/null || true
    fi
 
done < "$ACCESSION_LIST"
 
echo -e "\n=================================================="
echo "Pipeline execution finished successfully!"
echo "All outputs are located in: $OUTPUT_DIR"
echo "=================================================="