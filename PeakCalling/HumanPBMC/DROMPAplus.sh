
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
    cd ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/${histone}/
    
    # Cell types array
    local cell_types=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")
    
    # Process treatment BAMs
    for cell in "${cell_types[@]}"; do
        echo "Processing ${histone} treatment - $cell"
        parse2wig+ -i ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/${histone}/split_celltype_bams/${histone}_${cell}.bam \
                   -o ${histone}.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
                   -n GR
    done
    
    # Process input BAMs  
    for cell in "${cell_types[@]}"; do
        echo "Processing input - $cell"
        parse2wig+ -i ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/${histone}/split_celltype_bams/input_${cell}.bam \
                   -o input.${cell} \
                   --pair \
                   --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
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



I need to do this for all histone
histones=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
cd ~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27ac/

#HM1 H3K27ac With  INPUT
dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K27ac.B.100.bw, $dir/input.B.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.DC.100.bw, $dir/input.DC.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.Mono.100.bw, $dir/input.Mono.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.NK.100.bw, $dir/input.NK.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.otherT.100.bw, $dir/input.otherT.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.other.100.bw, $dir/input.other.100.bw, H3K27ac,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27ac/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K27ac.B.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.CD4T.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.CD8T.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.DC.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.Mono.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.NK.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.otherT.100.bw, H3K27ac,,,100 \
-i $dir/H3K27ac.other.100.bw, H3K27ac,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27ac/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5


cd /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27me3/
#HM2 H3K27me3 With  INPUT
dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K27me3.B.100.bw, $dir/input.B.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.DC.100.bw, $dir/input.DC.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.Mono.100.bw, $dir/input.Mono.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.NK.100.bw, $dir/input.NK.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.otherT.100.bw, $dir/input.otherT.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.other.100.bw, $dir/input.other.100.bw, H3K27me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27me3/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K27me3.B.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.CD4T.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.CD8T.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.DC.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.Mono.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.NK.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.otherT.100.bw, H3K27me3,,,100 \
-i $dir/H3K27me3.other.100.bw, H3K27me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K27me3/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5


cd /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me1/
#HM3 H3K4me1 With  INPUT
dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me1.B.100.bw, $dir/input.B.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.DC.100.bw, $dir/input.DC.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.Mono.100.bw, $dir/input.Mono.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.NK.100.bw, $dir/input.NK.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.otherT.100.bw, $dir/input.otherT.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.other.100.bw, $dir/input.other.100.bw, H3K4me1,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me1/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me1.B.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.CD4T.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.CD8T.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.DC.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.Mono.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.NK.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.otherT.100.bw, H3K4me1,,,100 \
-i $dir/H3K4me1.other.100.bw, H3K4me1,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me1/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5

cd /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me2/
#HM4 H3K4me2 With  INPUT
dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me2.B.100.bw, $dir/input.B.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.DC.100.bw, $dir/input.DC.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.Mono.100.bw, $dir/input.Mono.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.NK.100.bw, $dir/input.NK.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.otherT.100.bw, $dir/input.otherT.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.other.100.bw, $dir/input.other.100.bw, H3K4me2,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me2/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me2.B.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.CD4T.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.CD8T.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.DC.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.Mono.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.NK.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.otherT.100.bw, H3K4me2,,,100 \
-i $dir/H3K4me2.other.100.bw, H3K4me2,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me2/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5



cd /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me3/
#HM5 H3K4me3 With  INPUT
dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me3.B.100.bw, $dir/input.B.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.DC.100.bw, $dir/input.DC.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.Mono.100.bw, $dir/input.Mono.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.NK.100.bw, $dir/input.NK.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.otherT.100.bw, $dir/input.otherT.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.other.100.bw, $dir/input.other.100.bw, H3K4me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me3/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_SHARP \
-i $dir/H3K4me3.B.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.CD4T.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.CD8T.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.DC.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.Mono.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.NK.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.otherT.100.bw, H3K4me3,,,100 \
-i $dir/H3K4me3.other.100.bw, H3K4me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K4me3/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5


cd /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K9me3/
#HM6 H3K9me3 With  INPUT
dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K9me3.B.100.bw, $dir/input.B.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.CD4T.100.bw, $dir/input.CD4T.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.CD8T.100.bw, $dir/input.CD8T.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.DC.100.bw, $dir/input.DC.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.Mono.100.bw, $dir/input.Mono.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.NK.100.bw, $dir/input.NK.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.otherT.100.bw, $dir/input.otherT.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.other.100.bw, $dir/input.other.100.bw, H3K9me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K9me3/peakbed_with_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3


dir=parse2wigdir+
drompa+ PC_BROAD \
-i $dir/H3K9me3.B.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.CD4T.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.CD8T.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.DC.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.Mono.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.NK.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.otherT.100.bw, H3K9me3,,,100 \
-i $dir/H3K9me3.other.100.bw, H3K9me3,,,100 \
-o /home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed/H3K9me3/peakbed_without_input/DROMPAplus \
--gt /home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt \
-g /home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt \
--lpp 5 --showitag 1 --callpeak --pthre_internal 5



MARKS=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
celltypes=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

# Loop over all marks
for cell in "${celltypes[@]}"; do
for MARK in "${MARKS[@]}"; do
    input_file="DROMPAplus.${MARK}.${cell}.100.bw.peak.bed"
    output_file="DROMPAplus_${MARK}_${cell}.sorted.bed"
    
    if [[ -f "$input_file" ]]; then
        # Sort by chr (1st col) and start (2nd col)
        sort -k1,1 -k2,2n "$input_file" > "$output_file"
        echo "✅ Sorted and saved: $input_file → $output_file"
    else
        echo "⚠️  Warning: $input_file not found!"
    fi
done
done



#!/bin/bash

####################################################################################################################
### DROMPA+ Peak Calling for All Histones
### Simple function for peak calling only
####################################################################################################################

# Configuration
BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/splitbam_realbam"
PEAKBED_DIR="${BASE_DIR}/HumanPBMC_peakbed/DROMPAplus_peakbed"
GENOME_TABLE="/home/wahid/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt"
REF_FLAT="/home/wahid/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt"

# Histones array
HISTONES=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")

# Cell types array
CELL_TYPES=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

####################################################################################################################
### Peak Calling WITH INPUT
####################################################################################################################

run_peak_calling_with_input() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Peak Calling WITH INPUT for $histone"
    echo "####################################################################################################################"
    
    cd "${PEAKBED_DIR}/${histone}/"
    
    local dir="parse2wigdir+"
    local output_dir="peakbed_with_input"
    mkdir -p "$output_dir"
    
    echo "Running DROMPA+ with input control for $histone..."
    
    drompa+ PC_BROAD \
        -i $dir/${histone}.B.100.bw,$dir/input.B.100.bw,${histone},,,100 \
        -i $dir/${histone}.CD4T.100.bw,$dir/input.CD4T.100.bw,${histone},,,100 \
        -i $dir/${histone}.CD8T.100.bw,$dir/input.CD8T.100.bw,${histone},,,100 \
        -i $dir/${histone}.DC.100.bw,$dir/input.DC.100.bw,${histone},,,100 \
        -i $dir/${histone}.Mono.100.bw,$dir/input.Mono.100.bw,${histone},,,100 \
        -i $dir/${histone}.NK.100.bw,$dir/input.NK.100.bw,${histone},,,100 \
        -i $dir/${histone}.otherT.100.bw,$dir/input.otherT.100.bw,${histone},,,100 \
        -i $dir/${histone}.other.100.bw,$dir/input.other.100.bw,${histone},,,100 \
        -o ${output_dir}/DROMPAplus \
        --gt "$GENOME_TABLE" \
        -g "$REF_FLAT" \
        --lpp 5 --showitag 1 --callpeak --pthre_internal 4 --pthre_enrich 3
    
    echo "✓ Completed peak calling with input for $histone"
    echo ""
}

####################################################################################################################
### Peak Calling WITHOUT INPUT
####################################################################################################################

run_peak_calling_without_input() {
    local histone="$1"
    
    echo "####################################################################################################################"
    echo "### Peak Calling WITHOUT INPUT for $histone"
    echo "####################################################################################################################"
    
    cd "${PEAKBED_DIR}/${histone}/"
    
    local dir="parse2wigdir+"
    local output_dir="peakbed_without_input"
    mkdir -p "$output_dir"
    
    echo "Running DROMPA+ without input control for $histone..."
    
    drompa+ PC_BROAD \
        -i $dir/${histone}.B.100.bw,${histone},,,100 \
        -i $dir/${histone}.CD4T.100.bw,${histone},,,100 \
        -i $dir/${histone}.CD8T.100.bw,${histone},,,100 \
        -i $dir/${histone}.DC.100.bw,${histone},,,100 \
        -i $dir/${histone}.Mono.100.bw,${histone},,,100 \
        -i $dir/${histone}.NK.100.bw,${histone},,,100 \
        -i $dir/${histone}.otherT.100.bw,${histone},,,100 \
        -i $dir/${histone}.other.100.bw,${histone},,,100 \
        -o ${output_dir}/DROMPAplus \
        --gt "$GENOME_TABLE" \
        -g "$REF_FLAT" \
        --lpp 5 --showitag 1 --callpeak --pthre_internal 5
    
    echo "✓ Completed peak calling without input for $histone"
    echo ""
}

####################################################################################################################
### Complete Peak Calling for Single Histone
####################################################################################################################

run_peak_calling_for_histone() {
    local histone="$1"
    
    echo "===================================================================="
    echo "STARTING PEAK CALLING FOR: $histone"
    echo "===================================================================="
    
    # Check if parse2wigdir+ exists
    if [[ ! -d "${PEAKBED_DIR}/${histone}/parse2wigdir+" ]]; then
        echo "Error: parse2wigdir+ not found for $histone. Run parse2wig+ first."
        return 1
    fi
    
    # With input
    run_peak_calling_with_input "$histone"
    
    # Without input
    run_peak_calling_without_input "$histone"
    
    echo "===================================================================="
    echo "COMPLETED PEAK CALLING FOR: $histone"
    echo "===================================================================="
    echo ""
}

####################################################################################################################
### Run for All Histones
####################################################################################################################

run_peak_calling_all_histones() {
    echo "####################################################################################################################"
    echo "### DROMPA+ PEAK CALLING FOR ALL HISTONES"
    echo "### Histones: ${HISTONES[*]}"
    echo "####################################################################################################################"
    echo ""
    
    for histone in "${HISTONES[@]}"; do
        run_peak_calling_for_histone "$histone"
    done
    
    echo "####################################################################################################################"
    echo "### ALL PEAK CALLING COMPLETED!"
    echo "### Processed: ${HISTONES[*]}"
    echo "####################################################################################################################"
}

####################################################################################################################
### Individual Functions
####################################################################################################################

run_only_with_input() {
    local histone="$1"
    run_peak_calling_with_input "$histone"
}

run_only_without_input() {
    local histone="$1"
    run_peak_calling_without_input "$histone"
}

####################################################################################################################
### Main Execution
####################################################################################################################

main() {
    local command="${1:-all}"
    local histone="${2:-}"
    
    case "$command" in
        "all")
            run_peak_calling_all_histones
            ;;
        "H3K27ac"|"H3K27me3"|"H3K4me1"|"H3K4me2"|"H3K4me3"|"H3K9me3")
            run_peak_calling_for_histone "$command"
            ;;
        "with_input")
            if [[ -n "$histone" ]]; then
                run_only_with_input "$histone"
            else
                echo "Error: Histone name required"
            fi
            ;;
        "without_input")
            if [[ -n "$histone" ]]; then
                run_only_without_input "$histone"
            else
                echo "Error: Histone name required"
            fi
            ;;
        "help"|"-h"|"--help")
            echo "Usage: $0 [OPTION] [HISTONE]"
            echo ""
            echo "Options:"
            echo "  all                    Run peak calling for all histones"
            echo "  H3K27ac               Run for H3K27ac"
            echo "  with_input [HISTONE]   Run only with input control"
            echo "  without_input [HISTONE] Run only without input control"
            echo ""
            echo "Examples:"
            echo "  $0 all                    # All histones"
            echo "  $0 H3K27ac               # Only H3K27ac"
            echo "  $0 with_input H3K27ac    # Only with input for H3K27ac"
            ;;
        *)
            echo "Error: Unknown command '$command'"
            echo "Use '$0 help' for usage information"
            exit 1
            ;;
    esac
}

# Run if executed directly
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi
