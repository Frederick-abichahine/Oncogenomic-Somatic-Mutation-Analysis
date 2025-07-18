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

### Somatic Copy Number Calling  

Somatic copy number alterations were identified using VarScan2 (v2.3.9) copynumber and copyCaller modules on tumor and matched control samples processed with samtools mpileup with a minimum mapping quality of 1, using the human_g1k_v37 reference genome [8]. Segmentation analysis was performed in R using the DNAcopy package, where adjusted log2 ratios were smoothed and segmented by Circular Binary Segmentation (CBS) with a minimum segment width of 2 probes and an undo standard deviation threshold of 3 [9]. Final segmented copy number profiles were exported for downstream analyses.  

### Variant Calling and Annotation  

Germline variants were called from aligned sequencing data using bcftools (v.1.21) mpileup with the human_g1k_v37 reference genome, followed by bcftools call in consensus mode, and additionally with GATK UnifiedGenotyper [10]. Resulting VCF files were filtered using vcftools (v0.1.17) to retain variants with minimum quality score ≥20, mean read depth between 5 and 200, and excluding indels. Variants were annotated using snpEff (v4.3) to predict functional impacts, and snpSift Annotate to integrate data from HapMap and ClinVar Pathogenic databases, prioritizing known pathogenic variants [11]. High-impact variants with read depth >20 and existing IDs were filtered using snpSift filter for downstream analysis. Annotated variants of interest were further inspected using Integrative Genomics Viewer (IGV) (v2.19.4) to visually validate their presence and allelic composition in both tumor and control samples [12].  

### Somatic Variant Calling and Annotation  

Somatic mutations were identified by comparing tumor and matched control BAM files using samtools mpileup and VarScan2 in somatic mode, generating SNP and indel calls in VCF format. Somatic VCF files were annotated with snpSift Annotate using HapMap. Variants were filtered to separate known SNPs (with rsID) and novel mutations, and specific position-based filters were applied where relevant for analysis. Additionally, mutational signature analysis was performed using the COSMIC Signature Profiler web application by uploading somatic SNVs in VCF format and selecting the targeted sequencing mode [13].  

### Tumor Purity and Ploidy Estimation  

Heterozygous SNPs were identified in control samples using bcftools, and allele-specific read counts were obtained with GATK ASEReadCounter (min depth 20, mapping and base quality ≥20). CLONETv2 R package was used to compute beta tables, ploidy, admixture, and allele-specific copy number profiles from segmented SCNA data and allelic fractions [14]. Tumor purity was estimated with TPES library in R by integrating somatic SNV read counts (min coverage 10, min alternate reads 10, max allele frequency 0.6) with CLONET-derived ploidy, comparing observed allele fractions of heterozygous somatic SNVs to the expected 0.47 [15].  

## Results and Discussion  

### Data Preprocessing  

After completing the full preprocessing pipeline, both samples showed a substantial reduction in read count. The control sample decreased from ~19.7 million reads to 12.9 million after realignment, and to 11.1 million post-deduplication. The tumor sample followed a similar pattern, dropping from 15 million to 9.6 million, then to 8.4 million. The largest losses occurred during realignment, where ~6.8 million (control) and ~5.4 million (tumor) reads were discarded, as expected, as realignment around indels corrects misalignments and exposes low-quality reads that are removed to reduce false positives in downstream analyses. Despite the reduction, mapping rates improved slightly (99.75% to 99.81% for control; 99.96% to 99.98% for tumor), and proper pairing remained high (>99.6%), indicating improved alignment quality. The final BAM files also contained very few singleton reads (0.09% control; 0.01% tumor) or interchromosomal mate pairs (0.03% control; 0.025% tumor), further confirming the integrity of the processed data for downstream analysis.  

### Somatic Copy Number Calling  

Analysis of somatic copy number alterations using CBS revealed a series of genomic segments with altered copy number states. The segmentation plot (Figure 1) displays the log₂ ratios of tumor versus control read depths across genomic regions, with values centered around zero indicating diploid regions. Multiple contiguous regions exhibited reduced log₂ ratios, consistent with somatic copy number losses. In contrast, very few and small segments showed mild elevations, suggesting low-level gains that were not substantial enough to qualify as high-confidence amplifications. Red horizontal lines denote CBS-detected segments, representing the average log₂ ratio across each region. Color alternation between black and green reflects transitions between chromosomal intervals, which, based on the input data, correspond to chromosomes 15 through 18. Collectively, the segmentation results point to localized genomic instability, primarily in the form of deletions that may be implicated in tumor development.  

<img alt="figure1" src="figures/report_figure_1.png">

### Variant Calling and Annotation  

Variant calling in tumor and control samples revealed substantial differences, with 14,121 and 15,992 discordant sites between the two calling methods, reflecting tool-specific sensitivity and filtering criteria. After applying stringent filters based on read depth, quality, and indel exclusion, 14,123 high-confidence variants remained in the tumor sample. Functional annotation identified 39 high-impact mutations in cancer-relevant genes, including ANXA7 (a tumor suppressor in prostate and brain cancers) [16], PTGES3 (an Hsp90 co-chaperone supporting oncogenic protein stability) [17], and CLTC (involved in chromosomal translocations in lymphomas) [18], suggesting potential roles in tumorigenesis. Lastly, pathway enrichment analysis revealed several significantly affected pathways, including the Aryl Hydrocarbon Receptor Pathway (WP2873), linked to immune regulation and tumor progression [19], and the Omega-6 Fatty Acids in Senescence Pathway (WP5424), which connects lipid signaling to inflammation in cancer (Supp. Table 1) [20].  

### Visualization of BRCA1 Variant Using IGV  

To identify clinically relevant variants, the annotated dataset was filtered for pathogenic entries in ClinVar, revealing a single pathogenic SNV in BRCA1 on chromosome 17q21.31. This C>A substitution at chr17:41,246,494 results in a nonsense mutation (p.Glu352Ter; p.E352*), introducing a premature stop codon that truncates the protein and disrupts its role in homologous recombination-mediated DNA repair, potentially driving tumorigenesis.  

Visualization in IGV (Figure 2) shows 63% reference (C) and 37% alternate (A) allele support in the control sample, indicating heterozygosity. In the tumor, this shifts to 80% A and 20% C, suggesting loss of heterozygosity and enrichment of the pathogenic allele. This allele imbalance supports the somatic relevance of the variant. Known as rs80357472, the mutation is linked to hereditary breast and ovarian cancer (HBOC) and is considered clinically actionable, particularly in the context of PARP inhibitor therapies [21].  

<img alt="figure2" src="figures/report_figure_2.png">

### Somatic Variant Calling and Annotation  

Paired tumor–control analysis identified 220 high-confidence somatic variants from over 14 million genomic positions, using stringent thresholds for coverage and allele frequency. Additionally, 12,511 germline variants and 2,882 loss of heterozygosity (LOH) events were detected. Missense mutations on chromosome 15 were particularly notable (Table 1). Mutational signature analysis showed high concordance between observed and reconstructed spectra (cosine similarity 0.966, correlation 0.941), based on 14,251 single base substitutions. Three COSMIC SBS signatures were identified: SBS5 (67.5%), a common clock-like, age-related process [22]; SBS1 (16.2%), also age-related via 5-methylcytosine deamination; and SBS54 (16.4%), a rare, poorly characterized signature suggesting a novel mutational process (Supp. Figure 3). Overall, the mutational landscape is largely shaped by aging (SBS1, SBS5), with a secondary contribution from SBS54. Tumor mutational burden was 5.1 mutations per megabase using WGS-mode analysis.  

<img alt="table1" src="tables/report_table_1.png">

### Tumor Purity and Ploidy Analysis  

Tumor purity and ploidy were estimated using SCNA-based CLONETv2 and SNV-based TPES approaches. CLONETv2 integrated log₂ ratio–based segmentation with allelic imbalance from heterozygous SNPs, estimating a ploidy of 2.27 and admixture of 0.27, corresponding to ~73% tumor purity. Although ploidy exceeded diploid, most SCNA segments clustered around (1,0) and (1,1) in beta-log₂R space (Figure 3), consistent with hemizygous deletions and diploid regions. A few segments near (2,0) and (3,0) indicated copy-neutral and gain-associated LOH, contributing to the elevated ploidy. No segments were detected at (2,1), (2,2), (3,1), (4,0), likely due to insufficient heterozygous SNPs for beta estimation. While ploidy estimation incorporates all segments via log2 ratios, admixture estimation is limited to segments with sufficient heterozygous SNPs for beta calculation [14]. This asymmetry may bias estimates if gain regions are underrepresented in the beta table, possibly explaining ploidy >2.0 and underestimation of purity.  

TPES analysis under conservative thresholds (minAltReads=10, minCov=10) identified 4 somatic SNVs with a dominant allele frequency peak at 0.334, corresponding to ~71% purity (0.334/0.47), closely aligning with CLONETv2 (Figure 4). A secondary AF peak at 0.244 was supported by a single SNV. Relaxing filters (minAltReads=8, minCov=8) revealed a high-AF SNV at 0.419 (Supp. Figure 4), suggesting an upper purity bound of ~89%, though based on a single variant and potentially affected by noise.  

Together, CLONETv2 and TPES provided concordant purity estimates (~71-73%) under conservative settings, with TPES hinting at a higher upper limit. The elevated ploidy and loss-dominant profile may reflect gain regions excluded from beta analysis. These results highlight the importance of SNP density and segment representation in purity estimation and support the complementary use of SCNA- and SNV-based methods for tumor characterization [23].  

<img alt="figure3" src="figures/report_figure_3.png">  

<img alt="figure4" src="figures/report_figure_4.png">

## Conclusion  

This study applied an integrated genomic analysis pipeline to tumor and matched control samples, enabling the identification of somatic mutations, structural alterations, and clinically actionable variants. Key findings included a pathogenic BRCA1 nonsense mutation and multiple high-impact alterations in cancer-associated genes, supported by functional annotation and visual validation. Mutational signature analysis further contextualized the tumor’s evolutionary history. While some variability in purity estimation was observed, the combined analytical approaches provided a robust framework for characterizing tumor-specific genomic changes and informing potential therapeutic avenues.  

## References  

[1] Hanahan D, Weinberg RA. Hallmarks of cancer: the next generation. Cell. 2011;144(5):646–674. doi:10.1016/j.cell.2011.02.013

[2] Stratton MR, Campbell PJ, Futreal PA. The cancer genome. Nature. 2009;458(7239):719–724. doi:10.1038/nature07943

[3] Mardis ER. Genome sequencing and cancer. Curr Opin Genet Dev. 2012;22(3): 245–250. doi:10.1016/j.gde.2012.03.005

[4] Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. Bioinformatics. 2009;25(16):2078–2079. doi:10.1093/bioinformatics/btp352

[5] McKenna A, Hanna M, Banks E, et al. The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. Genome Res. 2010;20(9):1297–1303. doi:10.1101/gr.107524.110

[6] The International HapMap Consortium. A haplotype map of the human genome. Nature. 2005;437(7063):1299–1320. doi:10.1038/nature04226

[7] Picard Toolkit. Broad Institute. “Picard Toolkit.” 2019. Available at: https://broadinstitute.github.io/picard/ (accessed July 2025)

[8] Koboldt DC, Zhang Q, Larson DE, Shen D, McLellan MD, Lin L, et al. VarScan 2: somatic mutation and copy number alteration discovery in cancer by exome sequencing. Genome Res. 2012;22(3):568–576. https://doi.org/10.1101/gr.129684.111

[9] Seshan VE, Olshen A. DNAcopy: DNA copy number data analysis. R package version 1.72.0. Available at: https://bioconductor.org/packages/release/bioc/html/DNAcopy.html (accessed July 2025)

[10] Danecek P, Bonfield JK, Liddle J, Marshall J, Ohan V, Pollard MO, et al. Twelve years of SAMtools and BCFtools. Gigascience. 2021;10(2):giab008. https://doi.org/10.1093/gigascience/giab008

[11] Cingolani P, Platts A, Wang LL, Coon M, Nguyen T, Wang L, et al. A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of Drosophila melanogaster strain w1118; iso-2; iso-3. Fly (Austin). 2012;6(2):80–92. https://doi.org/10.4161/fly.19695

[12] Robinson JT, Thorvaldsdóttir H, Winckler W, Guttman M, Lander ES, Getz G, Mesirov JP. Integrative genomics viewer. Nat Biotechnol. 2011;29(1):24–26. https://doi.org/10.1038/nbt.1754

[13] COSMIC Mutational Signature Profiler. Wellcome Sanger Institute. Available at: https://cancer.sanger.ac.uk/signatures/assignment/app/ (accessed July 2025)

[14] Prandi D, Baca SC, Romanel A, Barbieri CE, Mosquera JM, Fontugne J, et al. Unraveling the clonal hierarchy of somatic genomic aberrations. Genome Biol. 2014;15(8):439. doi:10.1186/s13059-014-0439-6

[15] Locallo A, Prandi D, Fedrizzi T, Demichelis F. TPES: tumor purity estimation from SNVs. Bioinformatics. 2019;35(21):4433–4435. doi:10.1093/bioinformatics/btz406

[16] Srivastava M, Bubendorf L, Srikantan V, Fossom L, Nolan L, Glasman M, et al. ANX7, a candidate tumor suppressor gene for prostate cancer. Proc Natl Acad Sci USA. 2001;98(8):4575–4580. doi:10.1073/pnas.071055798

[17] Wang H, Sun P, Yao R, Zhang W, Zhou X, Yao J, He K. Comprehensive pan-cancer analysis of PTGES3 and its prognostic role in hepatocellular carcinoma. Front Oncol. 2023;13:1158490. doi:10.3389/fonc.2023.1158490

[18] Saito M, Novak U, Piovan E, Basso K, Sumazin P, Schneider C, et al. BCL6 suppression of BCL2 via Miz1 and its disruption in diffuse large B cell lymphoma. Proc Natl Acad Sci USA. 2009;106(27):11294–11299. https://doi.org/10.1073/pnas.0903854106

[19] Gutiérrez-Vázquez C, Quintana FJ. Regulation of the immune response by the aryl hydrocarbon receptor. Immunity. 2018;48(1):19–33. https://doi.org/10.1016/j.immuni.2017.12.012

[20] Innes JK, Calder PC. Omega-6 fatty acids and inflammation. Prostaglandins Leukot Essent Fatty Acids. 2018;132:41–48. doi: 10.1016/j.plefa.2018.03.004

[21] Ghadirian P, Robidoux A, Nassif E, Martin G, Potvin C, Patocskai E, et al. Screening for BRCA1 and BRCA2 mutations among French-Canadian breast cancer cases attending an outpatient clinic in Montreal. Clin Genet. 2014;85(1):31–35. doi:10.1111/cge.12174

[22] Alexandrov LB, Nik-Zainal S, Wedge DC, Aparicio SAJR, Behjati S, Biankin AV, et al. Signatures of mutational processes in human cancer. Nature. 2013;500(7463): 415–421. doi:10.1038/nature12477

[23] Ciani Y, Fedrizzi T, Prandi D, Lorenzin F, Locallo A, Gasperini P, et al. Allele-specific genomic data elucidate the role of somatic gain and copy-number neutral loss of heterozygosity in cancer. Cell Syst. 2022;13(2):183–193.e7. doi:10.1016/j.cels.2021.10.001