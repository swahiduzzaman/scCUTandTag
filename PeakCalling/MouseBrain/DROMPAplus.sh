#!/bin/bash

####################################################################################################################
### Application DROMPAplus on MouseBrain data
### Flexible function for multiple histones
####################################################################################################################

run_drompa_for_histone() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Processing $histone"
    echo "####################################################################################################################"
    
    # Create and change to output directory
    mkdir -p ~/project_scHMTF/mm_processed_data/splitbam_realbam/MouseBrain_peakbed/DROMPAplus_peakbed/${histone}/
    cd ~/project_scHMTF/mm_processed_data/splitbam_realbam/MouseBrain_peakbed/DROMPAplus_peakbed/${histone}/
    
    # Mouse brain cell types based on your data
    local cell_types=("Astrocytes" "mOL" "OPC" "OEC" "VLMC" "Microglia" "Neurons1" "Neurons2" "Neurons3")
    
    # Process treatment BAMs
    for cell in "${cell_types[@]}"; do
        if [ -f ~/project_scHMTF/mm_processed_data/splitbam_realbam/${histone}/split_celltype_bams/${histone}_${cell}.bam ]; then
            echo "Processing ${histone} treatment - $cell"
            parse2wig+ -i ~/project_scHMTF/mm_processed_data/splitbam_realbam/${histone}/split_celltype_bams/${histone}_${cell}.bam \
                       -o ${histone}.${cell} \
                       --pair \
                       --gt ~/project_scHMTF/mm_processed_data/ref/genome_file.txt \
                       -n GR
        else
            echo "Warning: ${histone}_${cell}.bam not found, skipping..."
        fi
    done
    
    # Process input BAMs  
    for cell in "${cell_types[@]}"; do
        if [ -f ~/project_scHMTF/mm_processed_data/splitbam_realbam/${histone}/split_celltype_bams/input_${cell}.bam ]; then
            echo "Processing input - $cell"
            parse2wig+ -i ~/project_scHMTF/mm_processed_data/splitbam_realbam/${histone}/split_celltype_bams/input_${cell}.bam \
                       -o input.${cell} \
                       --pair \
                       --gt ~/project_scHMTF/mm_processed_data/ref/genome_file.txt \
                       -n GR
        else
            echo "Warning: input_${cell}.bam not found, skipping..."
        fi
    done
    
    echo "✓ Completed $histone processing"
    echo ""
}

# Main function to run for all histones
run_drompa_all_histones() {
    local histones=("H3K27ac" "H3K27me3" "H3K36me3" "H3K4me3" "Olig2" "Rad21")
    
    for histone in "${histones[@]}"; do
        run_drompa_for_histone "$histone"
    done
    
    echo "####################################################################################################################"
    echo "### ALL HISTONES PROCESSING COMPLETED!"
    echo "####################################################################################################################"
}

# AUTO-RUN WHEN SCRIPT IS EXECUTED
run_drompa_all_histones





