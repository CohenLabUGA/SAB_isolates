# Iron physiology and metabolism of model phytoplankton taxa in the South Atlantic Bight

## Project Background

This project paired physiology and transcriptomics to study the impact of iron availablility on phytoplankton isolates from the South Atlantic Bight (SAB). The SAB is a wide continental shelf ecosystem in the south east US that hosts diverse, productive phytoplankton communities. Iron availability decreases with distance from the shore, and physiochemical properties seperates the SAB into distinct shelf zones. The outer shelf is flanked by the Gulf Stream which regularly introduces nutrient rich North Atlantic Deep Water to the photic zone with latitudinal fluctuations in its position. These intrusions stimulate phytoplankton blooms that can grow into iron limitation because the retention time of SAB shelf water is ~9-11 weeks. 

A <i>C. closterium</i> and coccolithophores (<i>G. oceanica</i> and <i>G. huxleyi</i>) were isolated from the inner and outer shelf zone, respectively. Isolates were grown in steady state high and low iron treatments, with an addition 24 hr iron amendment treatment with half of the low iron treatment volume. The cultures were harvested at early exponention growth, after which physiological parameters were measured and RNA was extracted. 

This repository contains data and code used in the analysis of this project. 

## rawData 

All of the raw data collected for physiology is contained in **rawData**. 

`physio_exp.csv`:   
This file contains most of the physiological data. It contains the growth rate, Chl a concentration, cell size, cell count, pH of media before and after the experiment, as well as the maximum potential quantum yeild of PSII (Fv/Fm), the reoxidation time of the first quinone acceptor ($\tau Q_a$ ), and the functional absorption cross-section of PSII ($\sigma_{PSII}$) which were extracted from the FIRe output files. The full output from the FIRe is the `fire_exp.csv` file; these are results from dark-adapted samples. Results from the Actinic Light Source (ALS) run on the fire are in the `als_exp.csv` file; this contains the non-photochemical quenching (NPQ) data used in Figure 2F. physiology in the **figures** folder. 

#### rawData/histData

Maintainence cultures were kept in small volumes during the year prior to the experiment. Growth rate and FIRe measurements were also taken during this time and are found in the subfolder **histData**. Growth rates are in `histGrowth.csv`, FIRe output from dark-adapted samples such as Fv/Fm are in `histFire.csv`, and FIRe output using the ALS, such as NPQ, are in `histALS.csv`. These data were used in **Fig2_physiology** to make the dashed lines in the bar graphs.

>The sequences used in the transcriptomic analysis can be found at **link to ncbi**, and the coassemblies can be found at **link to zenodo**. 

## output

Summarized physiological data and descriptive data from sequencing and assembly used to make the tables in the **figures** folder. 
>Used in `makeTables.R` script under folder **r_scripts**. 

## r_scripts

 `physio_fig.R` :
 This file will clean up data, summarize, and run statistical tests on each physiology parameter, eventually creating the physiology figure.
 - Input: `physio_exp.csv` and `als_exp.csv`
 - Output: Physiology figure, raw csv files located in output folder

`sourmash.R`:
This script will create **Fig.4_sourmash** showing the Jaccard similarity coefficients between all samples in a clustered heat map. The correlations were calculated using the `sourmash.sh` bash script located in the **bash_scripts folder**.

`functions.R`: 
Functions repeatidly used in the 'physio_exp.R' script; these include loading/installing packages, statistical tests, or mutating dataframes. 

## bash_scripts

All bash scripts used to assemble a hybrid De Novo assembly for each organism, quantify translated genes, and functionally annotate the genes. 

#### Transcriptome assembly:

1. Initial quality assement with FastQC use `fastqc.sh`
2. Trim adapter sequences with Trimmomatic using `trimmomatic.sh`
3. Remove ribosomeal RNA (rRNA) with RiboDetector using `ribodetector.sh`
    - Although polyadenalated ends were selected for during sequencing to reduce rRNA, this was performed as an additional quality control step.
4. Quality check with FastQC and MultiQC using `fastqc.sh` and `multiqc.sh`
5. Reduce sequence redundancy with CD-HIT using `cd-hit.sh`
6. Check for potential contamination with sourmash using `sig.sh` then `gather.sh`
7. Assembly a hybrid De Novo transcriptome for each organism with rnaSPAdes usig `assemble.sh`
8. Assess transcriptome quality with rnaQUAST using `quast.sh`
9. Check transcriptome completeness with BUSCO using `busco.sh`

#### Quantification and functional annotation

1. Translate reads from nucleotide to open reading frames (ORFs) using TransDecoder with `transdecoder.sh`
2. Quantify gene expression with Salmon using `salmon.sh`
3. Functional annotation with eggNOG-mapper using `eggnog.sh`

#### Additional analysis

sourmash was used to calculate the Jaccard similarity between all samples and between transcriptomes, later used to create the figures Fig.3_sourmash and SFig.2_sourmashTranscriptome respectively. This process uses `sig.sh` and `compare.sh`. 

Lineage differences between diatoms and coccolithophores and among species of each taxa were found using `cox1.sh`. This script will output the COX1 gene for each isolate, which can later be used to construct the maximum likelihood phylogenetic tree in Fig.3_cox1Tree. 

`ncbi_blast.sh` was used to find the iron starvation induced proteins (ISIP) in each transcriptome. Kegg orthology does not include ISIPs although reference sequences are located on the NCBI website. These sequences were later added to the functional annotation output for differential expression analysis. 

## expression_analysis

Jupyter notebooks for the differential expression analysis, enrichment analysis, and making of heatmaps.

#### `Kegg_annotations.ipynb`
This notebook cleans the output of eggNOG-mapper and pulls out just the Kegg Ko annotations from each isolate. To make the counts of ORFs to a function/gene, this notebook creates a table for each isolate mapping the ORF ID to Kegg Ko's. Some ORFs were assigned multiple Ko's, so the final output creates a table with three columns, ORF ID, Kegg Ko, and Ko_iteration; where the ORFs are repeated for each unique assigned Ko, and Ko_iteration iterates from 1 to n unique Ko's assigned. These ORF-to-Ko tables are then used to create a list of all unique Kos assinged across all organism, including the name and symbol for that Ko. Additionally, this notebook uses the Kegg API to create tables for pathways of interest, where the ko, name, and symbol for each pathway is listed in the table. 
- Input:
  - eggNOG-mapper annotation files
- Output:
  - ORF-to-Ko table, with ORFs repeated for each unique Ko.
  - ko_def.csv, a table of all unique Kos assigned across organisms with name and symbol for each Ko
  - all_paths.csv a table of the Kegg pathways Photosynthesis, Nitrogen metabolism, and Carbon fixation in photosynthetic organisms, including the ko, name, and symbol for each pathway.

#### `DE_run.ipynb`
Differential expression anlysis with DESeq2. This notebook uses ORF counts from Salmon to compare expression of genes among iron treatments for each isolate. The analysis is first done on all ORFs quantified, then using Kegg annotations from eggNOG-mapper. The counts are split between the number of Kegg Ko's an ORF was assigned and summed to each Ko to compare expression at the Kegg Ko level.
- Input:
  - Salmon quant.sf file
  - Kegg notebook output, kegg_def.csv
- Output:
  - DESeq2 results comparing high and low iron, and iron amendment and low iron, with adj-pvalues, standard error, and log fold change of expression for each gene across the comparison.
  - Log2 fold shrinkage of DESeq2 results; these files are named lfc.4.Avl.csv where 4 is the organism and Avl is the comparison.
  - Variance transformed stabilized counts (VST counts); these files are named vsd.4.csv where 4 is the organism.
 
#### `DE_figs.ipynb`
Plotting the differential expression anlysis results. Using the output from `Kegg_annotations.ipynb` and `DE_run.ipynb`, this notebook creates the heatmaps **Fig.5_heatmap**, **SFig.3_orfHeatma**, and **SFig.4_carbonHeatmap** located in the **figures** folder.

### `clusterProfiler.ipynb`
Running a gene set enrichment analysis using the DESeq2 results. This notebook finds the pathways with are statistically over represented among significantly expressed genes and creates the enrichmnet dotplot **Fig6_gsea**



# Acknowledgements

Bioinformatic processing was performed with resources from the UGA Georgia Advanced Computing Resource Center (GACRC). We acknowledge funding support from UGA Skidaway Institute of Oceanography and National Science Foundation awards NSF-OPP #2240780 and NSF-OPP #2149071 to NRC during the time of writing.
