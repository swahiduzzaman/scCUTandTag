# bedGraphDir="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/celltypeswisebam_l1/HM1"



write the code for 

CELLTYPES=("B" "CD4_T" "CD8_T" "DC" "Mono" "NK" "other" "other_T")

b" "cd4t" "cd8t" "dc" "mono" "nk" "other" "othert

conda activate PeakCalling_analysis
#HM1=H3K27ac
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K27ac.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K27ac.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K27ac.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K27ac.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K27ac.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K27ac.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K27ac.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K27ac.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR

#HM2=H3K27me3
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K27me3.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K27me3.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K27me3.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K27me3.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K27me3.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K27me3.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K27me3.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K27me3.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR



#HM3=H3K4me1
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K4me1.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K4me1.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K4me1.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K4me1.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K4me1.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K4me1.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K4me1.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K4me1.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR



#HM4=H3K4me2
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K4me2.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K4me2.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K4me2.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K4me2.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K4me2.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K4me2.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K4me2.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K4me2.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR


#HM5=H3K4me3
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K4me3.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K4me3.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K4me3.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K4me3.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K4me3.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K4me3.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K4me3.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K4me3.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR



#HM6=H3K9me3
parse2wig+ -i control_not_B_sorted.bam  -o input.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_B_sorted.bam  -o H3K9me3.b  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD4_T_sorted.bam  -o input.cd4t  --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD4_T_sorted.bam  -o H3K9me3.cd4t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_CD8_T_sorted.bam  -o input.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_CD8_T_sorted.bam  -o H3K9me3.cd8t --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_DC_sorted.bam  -o input.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_DC_sorted.bam  -o H3K9me3.dc --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_Mono_sorted.bam  -o input.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_Mono_sorted.bam  -o H3K9me3.mono --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_NK_sorted.bam  -o input.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_NK_sorted.bam  -o H3K9me3.nk --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_not_other_sorted.bam  -o input.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_sorted.bam  -o H3K9me3.other --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i control_other_T_sorted.bam  -o input.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
parse2wig+ -i treatment_other_T_sorted.bam  -o H3K9me3.othert --pair --gt ~/project_scHMTF/GSE195725_processed_data/ref/genome_file.txt -n GR
