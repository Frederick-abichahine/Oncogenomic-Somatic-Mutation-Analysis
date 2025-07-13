# Detection and Interpretation of Clinically Relevant Somatic Genomic Aberrations in an Oncologic Patient
### _Computational Human Genomics Report_

---
## Frederick , FAC, Abi Chahine
> M.Sc. in Quantitative & Computational Biology, University of Trento, Italy
## Hala, HA, Alshaar
> M.Sc. in Quantitative & Computational Biology, University of Trento, Italy 
## Anastasia, AB, Bertova
> M.Sc. in Quantitative & Computational Biology, University of Trento, Italy

---  

## Goal  

To detect and interpret the potential clinical relevance of somatic genomic aberrations harboured in
the genome of an oncologic patient.

## Introduction  

Genomic instability is a hallmark of cancer, arising from both inherited and somatic mutations that drive uncontrolled cell growth and malignant transformation [1]. These alterations range from single-nucleotide variants (SNVs) to large chromosomal rearrangements and identifying them is crucial for understanding tumor biology and designing targeted therapies [2]. The advent of high-throughput sequencing technologies has revolutionized cancer genome analysis by enabling researchers to distinguish between inherited (germline) and tumor-specific (somatic) mutations through comparison of matched c and tumor DNA samples [3]. This study employs an integrated bioinformatics pipeline to analyze paired control and tumor samples from a single patient, aiming to identify somatic mutations and structural genomic changes that may underlie tumor progression and reveal potential therapeutic targets.

## Materials and Methods  

### Data Preprocessing  

Raw sequencing data underwent a multi-step preprocessing pipeline to optimize alignment quality and variant detection. BAM files were first sorted and indexed using samtools (v1.21) for efficient genomic access [4]. Local realignment around indels was then performed with GATK (v3.8) RealignerTargetCreator and IndelRealigner, restricted to captured regions via the -L flag, to correct misalignments and reduce false positives [5]. Base quality scores, often affected by sequencing bias, were recalibrated using GATK BaseRecalibrator and PrintReads, with known HapMap (v3.3) variants used to model and adjust error profiles [6]. Recalibration was evaluated with AnalyzeCovariates to ensure score accuracy. PCR duplicates were removed using Picard (v2.27.4) MarkDuplicates with REMOVE_DUPLICATES=true, followed by reindexing with samtools [7]. Final alignment metrics were confirmed using samtools flagstat, ensuring high-quality, deduplicated BAM files for downstream analyses.

