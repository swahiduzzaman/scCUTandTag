<img width="833" height="811" alt="image" src="https://github.com/user-attachments/assets/af2565e3-bcc2-4adf-8f5d-3fddda8f51ca" />





# scCTPEAK: a benchmark suite for scCUT&Tag peak-calling across diverse epigenomic contexts

**scCTPEAK** is a fully reproducible benchmarking framework designed to systematically evaluate seven widely used peak-calling algorithms on single-cell CUT&Tag (scCUT&Tag) datasets.
The pipeline benchmarks the following tools:

**DROMPAplus**

**Genrich**

**GoPeaks**

**HOMER**

**MACS2**

**SEACR**

**SICER2**

By applying these seven tools to multiple histone modifications and transcriptional factors (e.g., H3K27ac, H3K27me3, H3K36me3, H3K9me3, H3K4me1, H3K4me2, H3K4me3, Olig2, Rad21), the pipeline captures performance differences across broad, and sharp (narrow) epigenomic contexts.


## Key Features
1. Unified peak-calling interface
For each cell type, we constructed a pseudo-input by pooling fragments from all remaining cell types.
In this design, the input control for a given cell type is derived from the complementary cell populations. This cross–cell-type input strategy provides a balanced background signal, reduces cell-type–specific bias, and improves the robustness of peak calling in scCUT&Tag data.
