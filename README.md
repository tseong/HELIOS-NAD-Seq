# HELIOS-NAD-Seq

# HELIOS NAD-Seq Bioinformatics Pipeline

This repository contains the bioinformatics workflow for **HELIOS NAD-Seq**, a sequencing protocol to detect and quantify NAD-capped RNAs at the 5′ end of transcripts.  
The pipeline processes raw FASTQ files through demultiplexing, adapter trimming, alignment, and NAD-cap–specific analysis.  

---

## Features
- Demultiplexing of barcoded reads with **Cutadapt**
- Adapter and quality trimming (custom script, trimmomatic)
- Alignment to reference genome with **bowtie2**
- Filtering and counting of NAD-capped vs control libraries
- Differential analysis of NAD-capping enrichment
- Time-course and condition-specific analysis normalization and visualization (e.g., growth curve experiments)

---

The environment file (`envs/helios.yml`) will install:
  - python=3.10
  - cutadapt=4.4
  - bowtie2=2.4.5
  - trimmomatic=0.39
  - subread
  - samtools
  - pysam
  - regex
  - matplotlib
  - pandas
  - numpy
  - scipy                  # for z-score & other stats
  - scikit-learn           # KMeans, z-scoring utilities
  - tslearn                # dynamic time warping clustering
  - gffutils
  - biopython
  - pip
  - seqtk
  - weblogo
  - seaborn
  - matplotlib-venn
  - pip:
      - pydeseq2==0.4.10
      - venn


---


