# Chitinophaga paper (ideARG - Corso, Iburg et al. 2026)

This repository contains the scripts used to perform the computational analyses described in:

> **[Full citation of the paper once published or submitted]**  
> [Authors]. *[Title of the paper]*. [Journal or preprint], [Year].  
> DOI: [to be added]

---

## 📁 Repository Overview

The repository includes scripts for both **read-based** and **assembly/MAG-based** analyses of environmental metagenomic samples, focusing on the detection and integration of antimicrobial resistance genes (ARGs), heavy metal resistance genes (HMRs), and mobile genetic elements (MGEs).

### 1. Read-Based Analysis (`args_oap`)
- SLURM job scripts (`AMR_soap`, `MRG_soap`) for large-scale screening of short reads using the **args_oap** tool.  
- Performs gene detection on raw reads against (respectively):
  - AMRFinderPlus database (for ARGs)
  - heavy metal resistance database BacMet (for MRGs)
- Parallelized execution on HPC clusters via job arrays.

### 2. Assembly- and MAG-Based Analysis (`Snakefile`)
- **Snakemake** workflow integrating several modules:
  - `prodigal`: ORF prediction on assemblies and MAGs.  
  - `DIAMOND for MRG`: alignment of predicted proteins against MRG databases BacMet.  
  - `mobileOG`: detection of mobile genetic elements.  
  - `hmmsearch`: screening against AMRFinderPlus HMM profiles (cut-GA).  
- Outputs summary tables for gene, contig, and sample-level annotations.

### 3. Integrative Analyses (`Python` and `R`)
- Python scripts for post-processing and merging HMM-based detections, creating `sample × HMM_accession` matrices, and integrating metadata.  
- R scripts and notebooks for:
  - mobileOG hits (by bitscore, identity, and coverage).  
  - MGE hits per contig.  
  - AMR and MGE annotations.  
  - Kaiju to Phyloseq - Merging with taxonomic assignments.  
