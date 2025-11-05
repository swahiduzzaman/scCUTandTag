#!/bin/bash
#######################################################################
# DROMPAplus Peak Calling: Compact Version
# Broad (PC_BROAD) and Sharp (PC_SHARP) peaks
# Works with or without input BAMs
#######################################################################

#------------------------------------------
# Histones and cell types
#------------------------------------------
histones=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me3")
cells=("B" "CD4T" "CD8T" "DC" "Mono" "NK" "otherT" "other")

#------------------------------------------
# Base directories
#------------------------------------------
BASE_DIR=~/project_scHMTF/GSE195725_processed_data/splitbam_realbam/HumanPBMC_peakbed/DROMPAplus_peakbed
REF_GENOME=~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt
REF_FLAT=~/project_scHMTF/GSE195725_processed_data/ref/refFlat.dupremoved.txt

#------------------------------------------
# Peak calling parameters
#------------------------------------------
LPP=5
SHOWITAG=1
CALLPEAK="--callpeak"
PTHRE_INTERNAL_WITH_INPUT=4
PTHRE_INTERNAL_NO_INPUT=5
PTHRE_ENRICH_WITH_INPUT=3

#------------------------------------------
# Function: run DROMPAplus
#------------------------------------------
run_drompa() {
    local histone=$1
    local peak_type=$2      # PC_BROAD or PC_SHARP
    local with_input=$3     # "yes" or "no"
    
    # cd into histone folder
    local histone_dir="${BASE_DIR}/${histone}"
    cd "$histone_dir" || { echo "⚠️ Cannot cd to $histone_dir"; return; }

    # parse2wigdir+ path inside this folder
    local PARSE_DIR="./parse2wigdir+"

    # output directory
    local out_dir="${histone_dir}/peakbed_${with_input}_${peak_type}/DROMPAplus"
    mkdir -p "$out_dir"

    # internal threshold
    local pthre_internal=$([ "$with_input" == "yes" ] && echo $PTHRE_INTERNAL_WITH_INPUT || echo $PTHRE_INTERNAL_NO_INPUT)
    local pthre_enrich=$([ "$with_input" == "yes" ] && echo $PTHRE_ENRICH_WITH_INPUT || echo 0)

    # build input arguments
    local args=()
    for cell in "${cells[@]}"; do
        if [ "$with_input" == "yes" ]; then
            args+=("-i $PARSE_DIR/${histone}.${cell}.100.bw, $PARSE_DIR/input.${cell}.100.bw, ${histone},,,100")
        else
            args+=("-i $PARSE_DIR/${histone}.${cell}.100.bw, ${histone},,,100")
        fi
    done

    # run DROMPAplus
    drompa+ "$peak_type" \
        "${args[@]}" \
        -o "$out_dir" \
        --gt "$REF_GENOME" \
        -g "$REF_FLAT" \
        --lpp $LPP \
        --showitag $SHOWITAG \
        $CALLPEAK \
        --pthre_internal $pthre_internal \
        $( [ "$with_input" == "yes" ] && echo "--pthre_enrich $pthre_enrich" )

    echo "✅ Finished $histone ($peak_type, with_input=$with_input)"
}

#------------------------------------------
# Main loop over histones
#------------------------------------------
for histone in "${histones[@]}"; do
    # Broad peaks
    run_drompa "$histone" "PC_BROAD" "yes"
    run_drompa "$histone" "PC_BROAD" "no"
    
    # Sharp peaks
    run_drompa "$histone" "PC_SHARP" "yes"
    run_drompa "$histone" "PC_SHARP" "no"
done
