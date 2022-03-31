# CGI\_nature\_micro\_2022
Code and commands used in the Shuqi &amp; Poulton et al. Nature Micro 2022 paper


## Requirements


### Python code:

 - MAGeCK:  (Version 0.5.7): https://sourceforge.net/p/mageck/wiki/install/
 - pandas (Version 1.4.1)
 - ast (Version 3.10.4)
 - snippy (Version 3.2-dev): https://github.com/tseemann/snippy

### R code:

 - R (Version 4.1.0)
 - readr (Version 1.4.0)
 - ggplot2 (Version 3.3.3)
 - reshape2 (Version 1.4.4)
 - gplots (Version 3.1.1)
 - dendextend (Version 1.14.0)
 - RColorBrewer (Version 1.1-2)

<br>
<br>

## Directory Structure:


### ./ 

Contains the scripts used to analyze counts and plot results.


### Results

Contains results and temporary files


### data


 - Final\_PreDepletion\_Grouped\_Split.csv: Table summarizing MAGeCK analysis results
 - list\_of\_H37Rv\_orfs\_in\_library.txt: List of genes in M. Tuberculosis strain H37Rv
 - Weiz\_wLFC.csv: Cleaned copy of Supplementary Data from Weizhen et al 2017.
 - **counts:**
 - - Raw counts obtained from the sequencing pipeline described in Bosch &amp; DeJesus et al. 2021.
 - **mageck:**
 - - Contains MAGeCK result files

  <br>
  <br>

## MAGECK Analysis

All MAGECK analysis were done using this command with the following flags:

    mageck test -k <path to file with combined control and experimental counts> -c 0,1,2 -t 3,4,5 -n <output file name> --gene-lfc-method alphamedian --norm-method control --control-sgrna <file with negative control IDs>


Gene-level result outputs are provided in the repository, as well as files with read-counts for each of the conditions.

<br>
<br>

## Snippy Analysis

Raw reads were downloaded from NCBI and run through snippy using the following command:

    snippy --ref <NC_018143.2.gb> --prefix <sample ID> --outdir <path to output directory> --pe1 <path to read 1 FASTQ> --pe2 <path to read 2 FASTQ>

  
  <br>
  <br>

## Heatmap

Read through the MAGeCK result files to generate a table (./Results/Heatmap/htmp\_data\_gt1Tmt\_D5.csv) of log-fold change values as described in the paper.

    python Data_Prep_Heatmap_Script.py

Take the log-fold change data table and create the Heatmap (./Results/Heatmap/htmp\_Ward\_D\_gt1Tmt\_D5.png) from Supplementary material.
   
    Rscript Heatmap_Script_D5.R

  
  <br>
  <br>

## TnSeq Comparison

Read cleaned supplementary data from Weizhen et al 2017 (./data/Weiz\_wLFC.csv) and the table summarizing MAGeCK results (./data/Final\_PreDepletion\_Grouped\_Split.csv) to generate Supplementary Data 2 (./Results/Supplemental\_Data\_2/Supplemental\_Data\_2.xlsx).

    python Gen_Weizhen_Overlap_Supplemental_Fig2.py

  
  
