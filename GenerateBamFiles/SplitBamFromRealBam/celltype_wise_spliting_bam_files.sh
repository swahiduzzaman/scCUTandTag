#!/bin/bash

############################################################################################################################################################
##### Splitting bam from real HumanPBMC bam for all histone marks
############################################################################################################################################################

# Define the histone marks to process
HISTONE_MARKS=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")

# Define cell types based on your PBMC data
CELL_TYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Base directory - use absolute path or $HOME
BASE_DIR="$HOME/project_scHMTF/GSE195725_processed_data/splitbam_realbam"

echo "####################################################################################################################"
echo "### STARTING BAM SPLITTING FOR ALL HISTONE MARKS"
echo "####################################################################################################################"
echo "Base directory: $BASE_DIR"

# Function to split BAM for a specific histone mark
split_bam_for_histone() {
    local histone="$1"
    local input_bam="${histone}.bam"
    local output_dir="${BASE_DIR}/${histone}/split_celltype_bams"
    
    echo "=== Processing $histone ==="
    echo "Input BAM: $input_bam"
    echo "Output directory: $output_dir"
    
    # Check if directory exists
    if [[ ! -d "${BASE_DIR}/${histone}" ]]; then
        echo "❌ ERROR: Directory ${BASE_DIR}/${histone} not found!"
        return 1
    fi
    
    # Check if input BAM exists
    if [[ ! -f "${BASE_DIR}/${histone}/$input_bam" ]]; then
        echo "❌ ERROR: Input BAM file $input_bam not found in ${BASE_DIR}/${histone}!"
        return 1
    fi
    
    # Change to histone directory
    cd "${BASE_DIR}/${histone}" || {
        echo "❌ ERROR: Cannot change to directory ${BASE_DIR}/${histone}"
        return 1
    }
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Check if BAM is indexed, if not create index
    if [[ ! -f "${input_bam}.bai" ]]; then
        echo "Creating index for $input_bam..."
        samtools index "$input_bam"
    fi
    
    # Process each cell type
    for celltype in "${CELL_TYPES[@]}"; do
        local barcode_file="${celltype}_barcodes.txt"
        
        # Check if barcode file exists
        if [[ ! -f "$barcode_file" ]]; then
            echo "❌ WARNING: Barcode file $barcode_file not found, skipping $celltype"
            continue
        fi
        
        echo "Splitting $histone for $celltype..."
        
        # Method 1: Using awk (most reliable)
        samtools view -h "$input_bam" | \
        awk -F'\t' -v barcode_file="$barcode_file" '
        BEGIN {
            # Load barcodes into array
            while((getline line < barcode_file) > 0) {
                barcodes[line] = 1
            }
            close(barcode_file)
        }
        /^@/ { 
            # Print all header lines
            print 
            next 
        }
        {
            # Check for CB tag in optional fields
            for(i=12; i<=NF; i++) {
                if($i ~ /^CB:Z:/) {
                    split($i, arr, ":")
                    cb_value = arr[3]
                    if(cb_value in barcodes) {
                        print
                        break
                    }
                }
            }
        }' | samtools view -b - > "${output_dir}/${histone}_${celltype}.bam"
        
        # Index the output BAM
        samtools index "${output_dir}/${histone}_${celltype}.bam"
        
        # Count reads and report
        read_count=$(samtools view -c "${output_dir}/${histone}_${celltype}.bam")
        echo "✓ $histone $celltype: $read_count reads"
        
    done
    
    echo "✓ Completed $histone splitting"
    echo ""
}

# Main function to process all histone marks
process_all_histones() {
    # Check if base directory exists
    if [[ ! -d "$BASE_DIR" ]]; then
        echo "❌ ERROR: Base directory $BASE_DIR does not exist!"
        exit 1
    fi
    
    for histone in "${HISTONE_MARKS[@]}"; do
        split_bam_for_histone "$histone"
    done
}

# Run quality check
quality_check() {
    echo "=== QUALITY CHECK ==="
    for histone in "${HISTONE_MARKS[@]}"; do
        local output_dir="${BASE_DIR}/${histone}/split_celltype_bams"
        if [[ -d "$output_dir" ]]; then
            echo "--- $histone ---"
            for celltype in "${CELL_TYPES[@]}"; do
                local bam_file="${output_dir}/${histone}_${celltype}.bam"
                if [[ -f "$bam_file" ]]; then
                    read_count=$(samtools view -c "$bam_file" 2>/dev/null || echo "0")
                    echo "  $celltype: $read_count reads"
                else
                    echo "  $celltype: NOT FOUND"
                fi
            done
        fi
        echo ""
    done
}

# Main execution
main() {
    echo "Starting BAM splitting process..."
    echo "Histone marks: ${HISTONE_MARKS[*]}"
    echo "Cell types: ${CELL_TYPES[*]}"
    echo ""
    
    # Process all histone marks
    process_all_histones
    
    # Quality check
    quality_check
    
    echo "####################################################################################################################"
    echo "### BAM SPLITTING COMPLETED!"
    echo "####################################################################################################################"
}

# Execute main function
main


############################################################################################################################################################
##### Splitting bam from real MouseBrain bam for all histone marks
############################################################################################################################################################

#!/bin/bash

############################################################################################################################################################
##### Splitting bam from real MouseBrain bam for all histone marks
############################################################################################################################################################

# Define the histone marks to process
HISTONE_MARKS=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")

# Define cell types for each histone mark
declare -A CELLTYPES_MAP=(
    ["H3K27ac"]="mOL Astrocytes VLMC OPC OEC"
    ["H3K27me3"]="mOL Astrocytes VLMC OPC OEC Microglia Neurons_1 Neurons_3"
    ["H3K36me3"]="mOL Astrocytes OPC OEC"
    ["H3K4me3"]="mOL Astrocytes VLMC OPC OEC Microglia Neurons_1 Neurons_2 Neurons_3"
    ["Olig2"]="mOL Astrocytes OEC Unknown"
    ["Rad21"]="mOL Astrocytes OEC Unknown"
)

# Base directory - use absolute path or $HOME
BASE_DIR="$HOME/project_scHMTF/GSE157637_processed_data/splitbam_realbam"  # Fixed for MouseBrain data

echo "####################################################################################################################"
echo "### STARTING BAM SPLITTING FOR ALL HISTONE MARKS - MOUSE BRAIN"
echo "####################################################################################################################"
echo "Base directory: $BASE_DIR"

# Function to split BAM for a specific histone mark
split_bam_for_histone() {
    local histone="$1"
    local input_bam="${histone}.bam"
    local output_dir="${BASE_DIR}/${histone}/split_celltype_bams"
    
    # Get cell types for this specific histone mark
    local cell_types=(${CELLTYPES_MAP[$histone]})
    
    if [[ ${#cell_types[@]} -eq 0 ]]; then
        echo "❌ ERROR: No cell types defined for $histone"
        return 1
    fi
    
    echo "=== Processing $histone ==="
    echo "Input BAM: $input_bam"
    echo "Output directory: $output_dir"
    echo "Cell types: ${cell_types[*]}"
    
    # Check if directory exists
    if [[ ! -d "${BASE_DIR}/${histone}" ]]; then
        echo "❌ ERROR: Directory ${BASE_DIR}/${histone} not found!"
        return 1
    fi
    
    # Check if input BAM exists
    if [[ ! -f "${BASE_DIR}/${histone}/$input_bam" ]]; then
        echo "❌ ERROR: Input BAM file $input_bam not found in ${BASE_DIR}/${histone}!"
        return 1
    fi
    
    # Change to histone directory
    cd "${BASE_DIR}/${histone}" || {
        echo "❌ ERROR: Cannot change to directory ${BASE_DIR}/${histone}"
        return 1
    }
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Check if BAM is indexed, if not create index
    if [[ ! -f "${input_bam}.bai" ]]; then
        echo "Creating index for $input_bam..."
        samtools index "$input_bam"
    fi
    
    # Process each cell type
    for celltype in "${cell_types[@]}"; do
        local barcode_file="${celltype}_barcodes.txt"
        
        # Check if barcode file exists
        if [[ ! -f "$barcode_file" ]]; then
            echo "❌ WARNING: Barcode file $barcode_file not found, skipping $celltype"
            continue
        fi
        
        echo "Splitting $histone for $celltype..."
        
        # Method 1: Using awk (most reliable)
        samtools view -h "$input_bam" | \
        awk -F'\t' -v barcode_file="$barcode_file" '
        BEGIN {
            # Load barcodes into array
            while((getline line < barcode_file) > 0) {
                barcodes[line] = 1
            }
            close(barcode_file)
        }
        /^@/ { 
            # Print all header lines
            print 
            next 
        }
        {
            # Check for CB tag in optional fields
            for(i=12; i<=NF; i++) {
                if($i ~ /^CB:Z:/) {
                    split($i, arr, ":")
                    cb_value = arr[3]
                    if(cb_value in barcodes) {
                        print
                        break
                    }
                }
            }
        }' | samtools view -b - > "${output_dir}/${histone}_${celltype}.bam"
        
        # Index the output BAM
        samtools index "${output_dir}/${histone}_${celltype}.bam"
        
        # Count reads and report
        read_count=$(samtools view -c "${output_dir}/${histone}_${celltype}.bam")
        echo "✓ $histone $celltype: $read_count reads"
        
    done
    
    echo "✓ Completed $histone splitting"
    echo ""
}

# Main function to process all histone marks
process_all_histones() {
    # Check if base directory exists
    if [[ ! -d "$BASE_DIR" ]]; then
        echo "❌ ERROR: Base directory $BASE_DIR does not exist!"
        exit 1
    fi
    
    for histone in "${HISTONE_MARKS[@]}"; do
        split_bam_for_histone "$histone"
    done
}

# Run quality check
quality_check() {
    echo "=== QUALITY CHECK ==="
    for histone in "${HISTONE_MARKS[@]}"; do
        local output_dir="${BASE_DIR}/${histone}/split_celltype_bams"
        local cell_types=(${CELLTYPES_MAP[$histone]})
        
        if [[ -d "$output_dir" ]]; then
            echo "--- $histone ---"
            for celltype in "${cell_types[@]}"; do
                local bam_file="${output_dir}/${histone}_${celltype}.bam"
                if [[ -f "$bam_file" ]]; then
                    read_count=$(samtools view -c "$bam_file" 2>/dev/null || echo "0")
                    echo "  $celltype: $read_count reads"
                else
                    echo "  $celltype: NOT FOUND"
                fi
            done
        else
            echo "--- $histone: OUTPUT DIRECTORY NOT FOUND ---"
        fi
        echo ""
    done
}

# Pre-run check to verify files exist
pre_run_check() {
    echo "=== PRE-RUN CHECK ==="
    for histone in "${HISTONE_MARKS[@]}"; do
        local cell_types=(${CELLTYPES_MAP[$histone]})
        echo "--- Checking $histone ---"
        
        # Check BAM file
        if [[ -f "${BASE_DIR}/${histone}/${histone}.bam" ]]; then
            echo "  ✓ BAM file exists"
        else
            echo "  ❌ BAM file missing: ${BASE_DIR}/${histone}/${histone}.bam"
        fi
        
        # Check barcode files
        for celltype in "${cell_types[@]}"; do
            if [[ -f "${BASE_DIR}/${histone}/${celltype}_barcodes.txt" ]]; then
                echo "  ✓ $celltype barcodes exist"
            else
                echo "  ❌ $celltype barcodes missing"
            fi
        done
        echo ""
    done
    echo "=== PRE-RUN CHECK COMPLETE ==="
    echo ""
}

# Main execution
main() {
    echo "Starting BAM splitting process for MouseBrain..."
    echo "Histone marks: ${HISTONE_MARKS[*]}"
    echo ""
    
    # Pre-run check
    pre_run_check
    
    # Ask for confirmation
    read -p "Continue with BAM splitting? (y/n): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
    
    # Process all histone marks
    process_all_histones
    
    # Quality check
    quality_check
    
    echo "####################################################################################################################"
    echo "### MOUSE BRAIN BAM SPLITTING COMPLETED!"
    echo "####################################################################################################################"
}

# Execute main function
main
