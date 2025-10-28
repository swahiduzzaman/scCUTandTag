#############################################################################################################################
#### Making BAM file for all HMs  for Human PBMC dataset
#### CellTypes wise
#### Author: Md Wahiduzzaman
#### This Code for all HMs 
#############################################################################################################################

#############################################################################################################################
## Required packages
#############################################################################################################################


library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(data.table)

#############################################################################################################################
## Make Sure the histones name
#############################################################################################################################

histones<-c("H3K27ac", "H3K27me3", "H3K4me1", "H3K4me2", "H3K4me3", "H3K9me3")

#############################################################################################################################
## Setting the "for loop" 
#############################################################################################################################

for(hist in histones){

histone <- readRDS(paste0("/home/wahid/project_scHMTF/hs_processed_data/", hist, ".rds"))
# --- 2. Read and create Fragment object ---
fragments <- CreateFragmentObject(
  path = paste0("/home/wahid/project_scHMTF/hs_processed_data/", hist, "_fragments.tsv.gz"),
  cells = colnames(histone)
)
# --- 3. Create new ChromatinAssay and reassign it to object ---
library(GenomeInfoDb)

# Read your genome file (chromosome sizes)
genome_file <- "/home/wahid/project_scHMTF/hs_processed_data/ref/hg38.genome"
chrom_sizes <- read.table(genome_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("seqnames", "seqlengths")

# Create header
header <- c("@HD\tVN:1.6\tSO:unsorted")
header <- c(header, paste0("@SQ\tSN:", chrom_sizes$seqnames, "\tLN:", chrom_sizes$seqlengths))
# Create Seqinfo manually
hg38_seqinfo <- Seqinfo(seqnames = chrom_sizes$seqnames,
                        seqlengths = chrom_sizes$seqlengths,
                        genome = "hg38")

                        # Now pass it to CreateChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = GetAssayData(histone[["tiles"]], slot = "counts"),
  fragments = fragments,
  sep = c("-", "-"),
  genome = hg38_seqinfo  # ✅ Custom genome
)

histone[["tiles"]] <- chrom_assay



# Path to original fragment file
frag_path <- Fragments(histone)[[1]]@path

# Read compressed fragment file
frags <- fread(
  cmd = paste("zcat", frag_path),
  header = FALSE,
  colClasses = c("character", "integer", "integer", "character", "integer")
)

colnames(frags) <- c("chr", "start", "end", "cell", "score")

# === 2. Output directory ===

## "/home/wahid/project_scHMTF/hs_processed_data/BAM/celltypeswisebam_l1/HM1"
outdir <- paste0("/home/wahid/project_scHMTF/hs_processed_data/BAM/celltypeswisebam_l1/",hist)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# === 4. Paired-end SAM writer ===
# === 5. Loop through clusters to write paired-end SAM files ===
# === 4. Paired-end SAM writer ===
write_pe_sam <- function(df, outfile) {
  sam_lines <- lapply(1:nrow(df), function(i) {
    row <- df[i]
    read_id <- paste0(row$cell, "_", i)
    chr <- row$chr
    start1 <- format(as.integer(row$start) + 1, scientific = FALSE)
    start2 <- format(as.integer(row$end) - 50, scientific = FALSE)
    if (as.integer(start2) < 1) start2 <- 1
    flag1 <- 99
    flag2 <- 147
    r1 <- paste(read_id, flag1, chr, start1, 255, "50M", "=", start2, as.integer(start2) - as.integer(start1), "*", "*", sep = "\t")
    r2 <- paste(read_id, flag2, chr, start2, 255, "50M", "=", start1, as.integer(start1) - as.integer(start2), "*", "*", sep = "\t")
    c(r1, r2)
  })
  writeLines(c(header, unlist(sam_lines)), con = outfile)
}

#Setting Idents instead of celltype.l1 table(H3K27me3$predicted.celltype.l1)
Idents(histone) <- histone$predicted.celltype.l1



# === 5. Iterate over all clusters ===
all_clusters <- unique(Idents(histone))
library(parallel)
# Delete old corrupted SAMs if needed
unlink(list.files(outdir, pattern = "\\.sam$", full.names = TRUE))


# Re-run parallel SAM generation
library(parallel)
mclapply(all_clusters, function(cluster) {
  message("Processing cluster: ", cluster)
  treatment_cells <- colnames(histone)[Idents(histone) == cluster]
  control_cells   <- colnames(histone)[Idents(histone) != cluster]

  frag_treat <- frags[cell %in% treatment_cells]
  frag_ctrl  <- frags[cell %in% control_cells]

  sam_treat <- file.path(outdir, paste0("treatment_", cluster, ".sam"))
  sam_ctrl  <- file.path(outdir, paste0("control_not_", cluster, ".sam"))

  write_pe_sam(frag_treat, sam_treat)
  write_pe_sam(frag_ctrl, sam_ctrl)
}, mc.cores = 64)
}


#############################################################################################################################
#### Making BAM file for all HMs and TFs for Mouse Brain dataset
#### CellTypes wise
#### Author: Md Wahiduzzaman
#### This Code for HMs and TFs are separately Applied[ 1st Part for HMs and 3nd Part for TFs]
#############################################################################################################################

####################
## 1st PART [HMs]
####################
#############################################################################################################################
## Required packages
#############################################################################################################################
library(Signac)
library(Seurat)
library(GenomicAlignments)
library(GenomeInfoDb)


histones<-c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3")
for(hist in histones){
# --- 1. Read in the Seurat object ---
histone <- readRDS(paste0("/home/wahid/project_scHMTF/mm_processed_data/", hist,"_seurat_object.Rds"))
histone <- UpdateSeuratObject(histone)

# --- 2. Read and create Fragment object ---
fragments <- CreateFragmentObject(
  path = paste0("/home/wahid/project_scHMTF/mm_processed_data/", hist, "_fragments.tsv.gz"),
  cells = colnames(histone)
)

# --- 3. Create new ChromatinAssay and reassign it to object ---


# Read your genome file (chromosome sizes)
genome_file <- "/home/wahid/project_scHMTF/mm_processed_data/ref/mm10.genome"
chrom_sizes <- read.table(genome_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("seqnames", "seqlengths")

# Create Seqinfo manually
mm10_seqinfo <- Seqinfo(seqnames = chrom_sizes$seqnames,
                        seqlengths = chrom_sizes$seqlengths,
                        genome = "mm10")

# Now pass it to CreateChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = GetAssayData(histone[["peaks"]], slot = "counts"),
  fragments = fragments,
  sep = c("-", "-"),
  genome = mm10_seqinfo  # ✅ Custom genome
)


histone[["peaks"]] <- chrom_assay
DefaultAssay(histone) <- "peaks"

library(data.table)

# Path to original fragment file
frag_path <- Fragments(histone)[[1]]@path

# Read compressed fragment file
frags <- fread(
  cmd = paste("zcat", frag_path),
  header = FALSE,
  colClasses = c("character", "integer", "integer", "character", "integer")
)

colnames(frags) <- c("chr", "start", "end", "cell", "score")

# === 2. Output directory ===
### NOTE for my work /home/wahid/project_scHMTF/mm_processed_data/result/sccuttag_cellcluster/HM3
### I did HM1=H3K27ac, HM2=H3K27me3, H3K36me3=HM3, HM4=H3K4me3

outdir <- paste0("/home/wahid/project_scHMTF/mm_processed_data/result/sccuttag_cellcluster/", hist)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


header <- c(
  "@HD\tVN:1.6\tSO:unsorted",
  "@SQ\tSN:chr1\tLN:195471971",
  "@SQ\tSN:chr2\tLN:182113224",
  "@SQ\tSN:chr3\tLN:160039680",
  "@SQ\tSN:chr4\tLN:156508116",
  "@SQ\tSN:chr5\tLN:151834684",
  "@SQ\tSN:chr6\tLN:149736546",
  "@SQ\tSN:chr7\tLN:145441459",
  "@SQ\tSN:chr8\tLN:129401213",
  "@SQ\tSN:chr9\tLN:124595110",
  "@SQ\tSN:chr10\tLN:130694993",
  "@SQ\tSN:chr11\tLN:122082543",
  "@SQ\tSN:chr12\tLN:120129022",
  "@SQ\tSN:chr13\tLN:120421639",
  "@SQ\tSN:chr14\tLN:124902244",
  "@SQ\tSN:chr15\tLN:104043685",
  "@SQ\tSN:chr16\tLN:98207768",
  "@SQ\tSN:chr17\tLN:94987271",
  "@SQ\tSN:chr18\tLN:90702639",
  "@SQ\tSN:chr19\tLN:61431566",
  "@SQ\tSN:chrX\tLN:171031299",
  "@SQ\tSN:chrY\tLN:91744698",
  "@SQ\tSN:chrM\tLN:16299"
)


# === 4. Paired-end SAM writer ===
# === 5. Loop through clusters to write paired-end SAM files ===
# === 4. Paired-end SAM writer ===
write_pe_sam <- function(df, outfile) {
  sam_lines <- lapply(1:nrow(df), function(i) {
    row <- df[i]
    read_id <- paste0(row$cell, "_", i)
    chr <- row$chr
    start1 <- format(as.integer(row$start) + 1, scientific = FALSE)
    start2 <- format(as.integer(row$end) - 50, scientific = FALSE)
    if (as.integer(start2) < 1) start2 <- 1
    flag1 <- 99
    flag2 <- 147
    r1 <- paste(read_id, flag1, chr, start1, 255, "50M", "=", start2, as.integer(start2) - as.integer(start1), "*", "*", sep = "\t")
    r2 <- paste(read_id, flag2, chr, start2, 255, "50M", "=", start1, as.integer(start1) - as.integer(start2), "*", "*", sep = "\t")
    c(r1, r2)
  })
  writeLines(c(header, unlist(sam_lines)), con = outfile)
}


# === 5. Iterate over all clusters ===
all_clusters <- unique(Idents(histone))
library(parallel)
# Delete old corrupted SAMs if needed
unlink(list.files(outdir, pattern = "\\.sam$", full.names = TRUE))

# Re-run parallel SAM generation
library(parallel)
mclapply(all_clusters, function(cluster) {
  message("Processing cluster: ", cluster)
  treatment_cells <- colnames(histone)[Idents(histone) == cluster]
  control_cells   <- colnames(histone)[Idents(histone) != cluster]

  frag_treat <- frags[cell %in% treatment_cells]
  frag_ctrl  <- frags[cell %in% control_cells]

  sam_treat <- file.path(outdir, paste0("treatment_", cluster, ".sam"))
  sam_ctrl  <- file.path(outdir, paste0("control_not_", cluster, ".sam"))

  write_pe_sam(frag_treat, sam_treat)
  write_pe_sam(frag_ctrl, sam_ctrl)
}, mc.cores = 64)
}


####################
## 2nd PART [TFs]
####################
#############################################################################################################################
## Required packages
#############################################################################################################################


library(Signac)
library(Seurat)
library(GenomicAlignments)
library(GenomeInfoDb)
library(data.table)

#############################################################################################################################
# specificed the TFs
#############################################################################################################################

transcriptionalfactors<-c("Olig2", "Rad21")

#############################################################################################################################
# Setting the for loop 
#############################################################################################################################

for(TF in transcriptionalfactors){
transcription <- readRDS(paste0("/home/wahid/project_scHMTF/mm_processed_data/", TF, "_seurat_object.Rds"))
transcription <- UpdateSeuratObject(transcription)

# --- 2. Read and create Fragment object ---
fragments <- CreateFragmentObject(
  path = paste0("/home/wahid/project_scHMTF/mm_processed_data/", TF, "_fragments.tsv.gz"),
  cells = colnames(transcription)
)

# --- 3. Create new ChromatinAssay and reassign it to object ---
# Read your genome file (chromosome sizes)
genome_file <- "/home/wahid/project_scHMTF/mm_processed_data/ref/mm10.genome"
chrom_sizes <- read.table(genome_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(chrom_sizes) <- c("seqnames", "seqlengths")

# Create Seqinfo manually
mm10_seqinfo <- Seqinfo(seqnames = chrom_sizes$seqnames,
                        seqlengths = chrom_sizes$seqlengths,
                        genome = "mm10")

# Now pass it to CreateChromatinAssay
chrom_assay <- CreateChromatinAssay(
  counts = GetAssayData(transcription[["peaks"]], slot = "counts"),
  fragments = fragments,
  sep = c("-", "-"),
  genome = mm10_seqinfo  # ✅ Custom genome
)

transcription[["peaks"]] <- chrom_assay
DefaultAssay(transcription) <- "peaks"

######################################
# Path to original fragment file
######################################

frag_path <- Fragments(transcription)[[1]]@path
######################################
# Read compressed fragment file
######################################
frags <- fread(
  cmd = paste("zcat", frag_path),
  header = FALSE,
  colClasses = c("character", "integer", "integer", "character", "integer")
)

colnames(frags) <- c("chr", "start", "end", "cell", "score")

######################################
# === 2. Output directory ===
######################################
outdir <- paste0("/home/wahid/project_scHMTF/mm_processed_data/result/sccuttag_cellcluster/", TF)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


header <- c(
  "@HD\tVN:1.6\tSO:unsorted",
  "@SQ\tSN:chr1\tLN:195471971",
  "@SQ\tSN:chr2\tLN:182113224",
  "@SQ\tSN:chr3\tLN:160039680",
  "@SQ\tSN:chr4\tLN:156508116",
  "@SQ\tSN:chr5\tLN:151834684",
  "@SQ\tSN:chr6\tLN:149736546",
  "@SQ\tSN:chr7\tLN:145441459",
  "@SQ\tSN:chr8\tLN:129401213",
  "@SQ\tSN:chr9\tLN:124595110",
  "@SQ\tSN:chr10\tLN:130694993",
  "@SQ\tSN:chr11\tLN:122082543",
  "@SQ\tSN:chr12\tLN:120129022",
  "@SQ\tSN:chr13\tLN:120421639",
  "@SQ\tSN:chr14\tLN:124902244",
  "@SQ\tSN:chr15\tLN:104043685",
  "@SQ\tSN:chr16\tLN:98207768",
  "@SQ\tSN:chr17\tLN:94987271",
  "@SQ\tSN:chr18\tLN:90702639",
  "@SQ\tSN:chr19\tLN:61431566",
  "@SQ\tSN:chrX\tLN:171031299",
  "@SQ\tSN:chrY\tLN:91744698",
  "@SQ\tSN:chrM\tLN:16299"
)


# === 4. Paired-end SAM writer ===
# === 5. Loop through clusters to write paired-end SAM files ===
# === 4. Paired-end SAM writer ===
write_pe_sam <- function(df, outfile) {
  sam_lines <- lapply(1:nrow(df), function(i) {
    row <- df[i]
    read_id <- paste0(row$cell, "_", i)
    chr <- row$chr
    start1 <- format(as.integer(row$start) + 1, scientific = FALSE)
    start2 <- format(as.integer(row$end) - 50, scientific = FALSE)
    if (as.integer(start2) < 1) start2 <- 1
    flag1 <- 99
    flag2 <- 147
    r1 <- paste(read_id, flag1, chr, start1, 255, "50M", "=", start2, as.integer(start2) - as.integer(start1), "*", "*", sep = "\t")
    r2 <- paste(read_id, flag2, chr, start2, 255, "50M", "=", start1, as.integer(start1) - as.integer(start2), "*", "*", sep = "\t")
    c(r1, r2)
  })
  writeLines(c(header, unlist(sam_lines)), con = outfile)
}


# === 5. Iterate over all clusters ===
# === 5. Iterate over all clusters ===
library(dplyr)

# Step 1: Assign initial celltype based on highest marker
transcription$celltype <- apply(
  transcription@meta.data[, grep("^marker_", colnames(transcription@meta.data))],
  1,
  function(row) {
    if (sum(row, na.rm = TRUE) == 0) {
      return("Unknown")
    } else {
      return(gsub("marker_", "", names(row)[which.max(row)]))  # highest marker
    }
  }
)

# Step 2: Merge into biologically valid groups
map_to_final_celltype <- function(ct) {
  if (ct == "Astrocytes") {
    return("Astrocytes")
  } else if (ct %in% c("OEC", "Pericytes", "VEC")) {
    return("OEC")
  } else if (ct %in% c("mOL", "OPC", "Oligodendrocytes", "COP.NFOL")) {
    return("mOL")
  } else if (ct == "Unknown") {
    return("Unknown")
  } else {
    return("Unknown")  # fallback
  }
}

transcription$merged_celltype <- sapply(transcription$celltype, map_to_final_celltype)

# Step 3: Ensure Astrocytes have ≥ 500 cells
astro_count <- sum(transcription$merged_celltype == "Astrocytes")

if (astro_count < 500) {
  extra_needed <- 500 - astro_count

  # Marker column for Astrocytes
  astro_marker_col <- "marker_Astrocytes"
  meta <- transcription@meta.data

  # Find cells from COP.NFOL and Pericytes with marker_Astrocytes > 0
  candidates <- which(
    transcription$celltype %in% c("COP.NFOL", "Pericytes") &
    meta[[astro_marker_col]] > 0
  )

  # Pick top cells with highest marker_Astrocytes
  top_astro_like <- head(
    candidates[order(-meta[candidates, astro_marker_col])],
    extra_needed
  )

  # Reassign these to Astrocytes
  transcription$merged_celltype[top_astro_like] <- "Astrocytes"
}

# Step 4: Final assignment
transcription$celltype <- transcription$merged_celltype

# Optional: View final counts
#table(Olig2_obj$celltype)
Idents(transcription)<- transcription$celltype
all_clusters <- unique(Idents(transcription))

library(parallel)
# Delete old corrupted SAMs if needed
unlink(list.files(outdir, pattern = "\\.sam$", full.names = TRUE))

# Re-run parallel SAM generation
library(parallel)
mclapply(all_clusters, function(cluster) {
  message("Processing cluster: ", cluster)
  treatment_cells <- colnames(transcription)[Idents(transcription) == cluster]
  control_cells   <- colnames(transcription)[Idents(transcription) != cluster]

  frag_treat <- frags[cell %in% treatment_cells]
  frag_ctrl  <- frags[cell %in% control_cells]

  sam_treat <- file.path(outdir, paste0("treatment_", cluster, ".sam"))
  sam_ctrl  <- file.path(outdir, paste0("control_not_", cluster, ".sam"))

  write_pe_sam(frag_treat, sam_treat)
  write_pe_sam(frag_ctrl, sam_ctrl)
}, mc.cores = 64)
}



