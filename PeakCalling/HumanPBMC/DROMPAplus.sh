
#!/bin/bash

####################################################################################################################
### Application DROMPAplus on HumanPBMC data
### Flexible function for multiple histones
####################################################################################################################

run_drompa_for_histone() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Processing $histone"
    echo "####################################################################################################################"
    
    # Create and change to output directory
    cd ~/project_scHMTF/hs_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/${histone}/
    
    # Cell types array
    local cell_types=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")
    
    # Process treatment BAMs
    for cell in "${cell_types[@]}"; do
        echo "Processing ${histone} treatment - $cell"
        parse2wig+ -i ~/project_scHMTF/hs_processed_data/splitbam_realbam/${histone}/split_celltype_bams/${histone}_${cell}.bam \
                   -o ${histone}.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/hs_processed_data/ref/genome_file.txt \
                   -n GR
    done
    
    # Process input BAMs  
    for cell in "${cell_types[@]}"; do
        echo "Processing input - $cell"
        parse2wig+ -i ~/project_scHMTF/hs_processed_data/splitbam_realbam/${histone}/split_celltype_bams/input_${cell}.bam \
                   -o input.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/Hs_processed_data/ref/genome_file.txt \
                   -n GR
    done
    
    echo "✓ Completed $histone processing"
    echo ""
}

# Main function to run for all histones
run_drompa_all_histones() {
    local histones=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
    
    for histone in "${histones[@]}"; do
        run_drompa_for_histone "$histone"
    done
    
    echo "####################################################################################################################"
    echo "### ALL HISTONES PROCESSING COMPLETED!"
    echo "####################################################################################################################"
}

# AUTO-RUN WHEN SCRIPT IS EXECUTED
run_drompa_all_histones

