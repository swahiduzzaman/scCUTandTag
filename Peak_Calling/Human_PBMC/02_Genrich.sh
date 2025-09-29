run_genrich_all_histone_marks() {
  # Configuration
  local BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/BAM/celltypeswisebam_l1"
  local OUT_BASE_DIR="/home/wahid/project_scHMTF/GSE195725_processed_data/result/l1_withcontrol"
  local CELLTYPES=("B" "CD4_T" "CD8_T" "DC" "Mono" "NK" "other" "other_T")
  local HISTONE_MARKS=("H3K27ac" "H3K27me3" "H3K4me1" "H3K4me2" "H3K4me3" "H3K9me2")
  local THREADS=32

  # Process each histone mark
  for MARK in "${HISTONE_MARKS[@]}"; do
    echo -e "\n🔬 ========== PROCESSING HISTONE MARK: $MARK =========="
    
    # Set mark-specific directories
    local MARK_DIR="${BASE_DIR}/${MARK}"
    local SORTED_DIR="${MARK_DIR}/qname_sorted"
    local OUT_DIR="${OUT_BASE_DIR}/${MARK}_peaks/Genrich"
    
    # Create directories
    mkdir -p "$OUT_DIR" "$SORTED_DIR" || {
      echo "❌ Failed to create directories for $MARK"
      continue
    }

    # Process each cell type for this histone mark
    for CELL in "${CELLTYPES[@]}"; do
      echo -e "\n📊 Processing $CELL for $MARK..."
      
      # Define file paths
      local RAW_TREATMENT="${MARK_DIR}/treatment_${CELL}.bam"
      local RAW_CONTROL="${MARK_DIR}/control_not_${CELL}.bam"
      local TREATMENT="${SORTED_DIR}/treatment_${CELL}_qnamesorted.bam"
      local CONTROL="${SORTED_DIR}/control_not_${CELL}_qnamesorted.bam"
      local OUTPUT_FILE="${OUT_DIR}/${CELL}_${MARK}_peaks.narrowPeak"
      local LOG_FILE="${OUT_DIR}/${CELL}_genrich.log"

      # Check input files
      if [[ ! -f "$RAW_TREATMENT" ]]; then
        echo "❌ Treatment BAM not found: $RAW_TREATMENT"
        continue
      fi
      if [[ ! -f "$RAW_CONTROL" ]]; then
        echo "❌ Control BAM not found: $RAW_CONTROL"
        continue
      fi

      # Sort by queryname (if not already done)
      if [[ ! -f "$TREATMENT" ]]; then
        echo "🔄 Sorting treatment BAM by queryname..."
        if ! samtools sort -n -o "$TREATMENT" -@ $THREADS "$RAW_TREATMENT" >> "$LOG_FILE" 2>&1; then
          echo "❌ Failed to sort treatment BAM for $CELL ($MARK)"
          continue
        fi
      fi

      if [[ ! -f "$CONTROL" ]]; then
        echo "🔄 Sorting control BAM by queryname..."
        if ! samtools sort -n -o "$CONTROL" -@ $THREADS "$RAW_CONTROL" >> "$LOG_FILE" 2>&1; then
          echo "❌ Failed to sort control BAM for $CELL ($MARK)"
          continue
        fi
      fi

      # Run Genrich with mark-specific parameters
      echo "🎯 Running Genrich peak calling..."
      if ! Genrich \
        -t "$TREATMENT" \
        -c "$CONTROL" \
        -o "$OUTPUT_FILE" \
        -j -y -r -v \
        -a 200 \
        -l 100 \
        -g 100 \
        -p 0.01 \
        -f BAM \
        >> "$LOG_FILE" 2>&1; then
        echo "❌ Genrich failed for $CELL ($MARK)"
        continue
      fi

      echo "✅ Successfully processed $CELL for $MARK"
      echo "📁 Peaks saved to: $OUTPUT_FILE"
    done

    # Generate summary report for this histone mark
    echo -e "\n📈 ========== SUMMARY FOR $MARK =========="
    local total_peaks=0
    local processed_cells=0
    
    for peakfile in "${OUT_DIR}"/*.narrowPeak; do
      if [[ -f "$peakfile" ]]; then
        count=$(wc -l < "$peakfile")
        echo "  $(basename "$peakfile"): $count peaks"
        total_peaks=$((total_peaks + count))
        processed_cells=$((processed_cells + 1))
      fi
    done
    
    echo "  Total: $processed_cells cell types processed, $total_peaks peaks called"
  done
}

# Run the function
run_genrich_all_histone_marks
