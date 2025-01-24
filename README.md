# Scripts for bulk RNA sequencing

- `bulk/trimmomaticTrim.sh` - trim fastq files
- `bulk/starAlign.sh` - align fastq files to produce a bam file
- `bulk/salmonIndex.sh` - create salmon-specific index for mapping
- `bulk/salmonMap.sh` - create gene quantification file from fastq files
- `bulk/importSalmon.R` - import salmon output to R
- `bulk/makeCm.R` - create count matrices with different normalisation schemes
- `bulk/dgeGsea.R` - calculate DGEs and GSEAs

# Scripts for single nuclei RNA & ATAC sequencing

- `sn/makeCellrangerDfs.R` - create csv files with sample names and fastq locations for cellranger
- `sn/cellrangerCount.sh` - create count matrices with cellranger-arc from fastq files
- `sn/functionsSeuratSignac.R` - generic Seurat and Signac functions in use for proccessing single nuclei RNA & ATAC data
- `sn/snPipeline.R` - apply function above, import cellranger matrices and procces then via the Seurat Signac framework

# Scripts for methylation EPIC850K illumina assay

- `meth/dmp.R` - import EPIC850K assay via minfi, perform quality control and calculate differentially methylated probes
- `meth/elmer.R` - ELMER pipeline for finding associations leveraging both bulk RNA and methylation data

# Scripts for multiplex immunohistochemistry proccessing

- `mIHC/qupathStitch.groovy` - stitch multiple inmixed tiff files as outputed by InForm
- `mIHC/stitch.sh` - submit stitching job to CECI with slurm
- `mIHC/qpathPipeline.groovy` - detect and classify cells in QuPath
- `mIHC/importQupath.R` - import QuPath results to R
- `mIHC/importInform.R` - import InForm results to R