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











The sequences used in the transcriptomic analysis can be found at **link to ncbi**, and the coassemblies can be found at **link to zenodo**. 

