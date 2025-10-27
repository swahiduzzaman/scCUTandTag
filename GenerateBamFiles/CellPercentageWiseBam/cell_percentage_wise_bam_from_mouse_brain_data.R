###############################################################################
#### FINAL: Mouse Brain Data Random Subsets for All HMs and TF
#### Author: Md. Wahiduzzaman | Date: 2025-10-11
#### Mouse Data: GSE157637 - mm10 genome 
#### H3K27ac , H3K27me3, H3K36me3, H3K4me3, Olig2 and Rad21 has # cells
#### 10414 , 13932, 4350, 13739, 4544, and 4375 respectively
###############################################################################

library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(data.table)


histones<-c("H3K27ac", "H3K27me3", "H3K36me3", "H3K4me3", "Olig2", "Rad21")
for(hist in histones){
# --- Configuration for MOUSE data ---
CONFIG <- list(
  seurat_rds = paste0("/home/wahid/project_scHMTF/GSE157637_processed_data/", hist, "_seurat_object.Rds"),
  fragment_path = paste0("/home/wahid/project_scHMTF/GSE157637_processed_data/", hist, "_fragments.tsv.gz"),
  genome_file = "/home/wahid/project_scHMTF/GSE157637_processed_data/ref/mm10.genome",
  output_dir = "/home/wahid/project_scHMTF/GSE157637_processed_data/BAM/random_subsets_final_all_histones",  # ✅ FIXED PATH
  seed = 2025,
  percent_steps = seq(10, 100, by = 10),
  min_fragment_length = 50,
  max_fragment_length = 1000
)

# --- Input Validation ---
validate_inputs <- function(config) {
  stopifnot("Seurat RDS not found" = file.exists(config$seurat_rds))
  stopifnot("Fragment file not found" = file.exists(config$fragment_path))
  stopifnot("Genome file not found" = file.exists(config$genome_file))
  
  dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)
  message(" All MOUSE input files validated")
}

# --- Load Data ---
load_data <- function(config) {
  message(" Loading Mouse Seurat object...")
  histone <- readRDS(config$seurat_rds)
  histone <- UpdateSeuratObject(histone)  # Update if needed
  
  # Get actual cells
  actual_cells <- colnames(histone)
  message(" Mouse Seurat object loaded: ", length(actual_cells), " cells")
  
  # Read fragments
  message(" Reading fragment file...")
  frags <- fread(
    cmd = paste("zcat", config$fragment_path),
    header = FALSE,
    col.names = c("chr", "start", "end", "cell", "score"),
    showProgress = TRUE
  )
  
  # Read chromosome sizes (mm10)
  chrom_sizes <- fread(config$genome_file, col.names = c("seqnames", "seqlengths"))
  
  # Filter for mouse chromosomes only
  mouse_chroms <- c(paste0("chr", 1:19), "chrX", "chrY", "chrM")
  chrom_sizes <- chrom_sizes[seqnames %in% mouse_chroms]
  
  return(list(
    actual_cells = actual_cells,
    fragments = frags,
    chrom_sizes = chrom_sizes
  ))
}

# --- Filter Fragments for MOUSE ---
filter_fragments <- function(frags, actual_cells, chrom_sizes, config) {
  original_count <- nrow(frags)
  message(" Initial fragments: ", original_count, " from ", length(unique(frags$cell)), " cells")
  
  # Mouse-specific filtering
  mouse_chroms <- chrom_sizes$seqnames
  
  frags_filtered <- frags[
    cell %in% actual_cells &
    chr %in% mouse_chroms &  # Mouse chromosomes only
    start >= 0 &
    end > start &
    (end - start) >= config$min_fragment_length &
    (end - start) <= config$max_fragment_length &
    score >= 0
  ]
  
  filtered_count <- original_count - nrow(frags_filtered)
  message(" Filtered fragments: ", nrow(frags_filtered), " (removed ", filtered_count, " invalid)")
  message(" Unique cells after filtering: ", length(unique(frags_filtered$cell)))
  message(" Mouse chromosomes used: ", paste(mouse_chroms, collapse = ", "))
  
  return(frags_filtered)
}

# --- Generate SAM Header for MOUSE ---
generate_sam_header <- function(chrom_sizes) {
  c(
    "@HD\tVN:1.6\tSO:unsorted",
    paste0("@SQ\tSN:", chrom_sizes$seqnames, "\tLN:", chrom_sizes$seqlengths),
    "@PG\tID:mouse_random_subsets\tPN:R_script\tVN:1.0"
  )
}

# --- Perfect SAM Writer (SAME) ---
write_perfect_sam <- function(fragments, output_file, header) {
  message(" Writing SAM: ", basename(output_file))
  
  total_frags <- nrow(fragments)
  if (total_frags == 0) {
    warning(" No fragments to write for ", output_file)
    return(0)
  }
  
  con <- file(output_file, "w")
  on.exit(close(con))
  
  writeLines(header, con)
  
  written <- 0
  chunk_size <- 50000
  
  for (start_idx in seq(1, total_frags, by = chunk_size)) {
    end_idx <- min(start_idx + chunk_size - 1, total_frags)
    chunk <- fragments[start_idx:end_idx, ]
    
    sam_lines <- character()
    
    for (i in 1:nrow(chunk)) {
      row <- chunk[i, ]
      global_idx <- start_idx + i - 1
      
      # Calculate coordinates
      start1 <- as.integer(row$start) + 1L
      start2 <- as.integer(row$end) - 49L  # 50bp reads
      
      # Validate coordinates
      if (start1 < 1L || start2 < 1L || start1 >= start2) next
      
      tlen <- as.integer(row$end - row$start)
      
      # Safe read name
      safe_cell <- gsub("[^A-Za-z0-9_-]", "_", row$cell)
      read_id <- paste0("frag_", safe_cell, "_", global_idx)
      
      # Sequences
      seq_50bp <- strrep("N", 50)
      qual_50bp <- strrep("I", 50)
      
      # Read 1
      r1 <- paste(
        read_id, "99", row$chr, start1, "255", "50M", 
        "=", start2, tlen, seq_50bp, qual_50bp,
        sep = "\t"
      )
      
      # Read 2
      r2 <- paste(
        read_id, "147", row$chr, start2, "255", "50M", 
        "=", start1, -tlen, seq_50bp, qual_50bp,
        sep = "\t"
      )
      
      sam_lines <- c(sam_lines, r1, r2)
      written <- written + 1
    }
    
    if (length(sam_lines) > 0) {
      writeLines(sam_lines, con)
    }
    
    if (end_idx %% 100000 == 0) {
      message("    Progress: ", end_idx, "/", total_frags, " fragments")
    }
  }
  
  message(" Finished: ", written, " fragments written")
  return(written)
}

# --- Validate and Convert SAM ---
validate_and_convert <- function(sam_file) {
  message(" Validating SAM file...")
  
  validation_result <- system2("samtools", c("quickcheck", sam_file), 
                              stdout = NULL, stderr = NULL)
  
  if (validation_result != 0) {
    message("  ❌ SAM validation failed")
    return(FALSE)
  }
  
  message(" SAM validation passed")
  
  # Convert to BAM
  bam_file <- sub("\\.sam$", ".bam", sam_file)
  message(" Converting to BAM...")
  
  convert_result <- system2("samtools", c("view", "-Sb", sam_file, "-o", bam_file))
  
  if (convert_result == 0 && file.exists(bam_file)) {
    read_count <- system2("samtools", c("view", "-c", bam_file), stdout = TRUE)
    message("   BAM created: ", trimws(read_count), " reads")
    return(TRUE)
  } else {
    message("  ❌ BAM conversion failed")
    return(FALSE)
  }
}

# --- Main Execution ---
main <- function() {
  message(paste0(" Starting","",hist, "","random subset generation..."))
  
  # Validate inputs
  validate_inputs(CONFIG)
  
  # Load data
  data <- load_data(CONFIG)
  
  # Filter fragments
  fragments <- filter_fragments(data$fragments, data$actual_cells, data$chrom_sizes, CONFIG)
  
  # Generate header
  sam_header <- generate_sam_header(data$chrom_sizes)
  
  # Setup random sampling
  set.seed(CONFIG$seed)
  all_cells <- data$actual_cells
  total_cells <- length(all_cells)
  
  message(" Total MOUSE cells for sampling: ", total_cells)
  
  # Show subset sizes
  message(" Subset sizes:")
  for (p in CONFIG$percent_steps) {
    n_cells <- round((p / 100) * total_cells)
    message(sprintf("  %3d%% → %5d cells", p, n_cells))
  }
  
  # Generate subsets
  results <- list()
  
  for (percent in CONFIG$percent_steps) {
    message("\n", strrep("=", 50))
    message("Processing MOUSE ", percent, "% subset")
    message(strrep("=", 50))
    
    # Calculate cell count
    n_cells <- round((percent / 100) * total_cells)
    sampled_cells <- sample(all_cells, n_cells)
    
    # Subset fragments
    subset_frags <- fragments[cell %in% sampled_cells]
    subset_frags <- subset_frags[order(chr, start, end)]
    
    message("  Cells: ", n_cells, " | Fragments: ", nrow(subset_frags))
    
    # Generate output file
    sam_file <- file.path(CONFIG$output_dir, sprintf(paste0(hist, "_%dcell.sam"), percent))
    
    # Write SAM file
    fragments_written <- write_perfect_sam(subset_frags, sam_file, sam_header)
    
    # Validate and convert
    validation_success <- validate_and_convert(sam_file)
    
    # Store results
    results[[as.character(percent)]] <- list(
      percent = percent,
      cells = n_cells,
      fragments_input = nrow(subset_frags),
      fragments_written = fragments_written,
      sam_file = sam_file,
      valid = validation_success
    )
  }
  
  # Print final summary
  print_summary(results)
  invisible(results)
}

# --- Summary Report ---
print_summary <- function(results) {
  message("\n", strrep("=", 70))
  message("FINAL SUMMARY - MOUSE H3K27ac")
  message(strrep("=", 70))
  
  for (percent in names(results)) {
    res <- results[[percent]]
    status <- ifelse(res$valid, " SUCCESS", "❌ FAILED")
    
    message(sprintf("  %3s%%: %5d cells | %7d fragments | %s",
                    percent, res$cells, res$fragments_written, status))
  }
  
  total_written <- sum(sapply(results, function(x) x$fragments_written))
  message("\n Total fragments written: ", format(total_written, big.mark = ","))
  message(" Output directory: ", CONFIG$output_dir)
  message(" MOUSE random subsets generated successfully!")
}

# --- Execute ---
if (!interactive()) {
  results <- main()
}
}
