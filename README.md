# CGI\_nature\_micro\_2022
Code used in the Shuqi &amp; Poulton et al. Nature Micro 2022 paper



## Requirements


### Python code:

 - MAGeCK:  (Version 0.5.7) https://sourceforge.net/p/mageck/wiki/install/
 - pandas (Version 1.4.1)
 - ast (Version 3.10.4)


### R code:

 - R (Version 4.1.0)
 - readr (Version 1.4.0)
 - ggplot2 (Version 3.3.3)
 - reshape2 (Version 1.4.4)
 - gplots (Version 3.1.1)
 - dendextend (Version 1.14.0)
 - RColorBrewer (Version 1.1-2)



## Directory Structure:


### ./ 

Contains the scripts used to analyze counts and plot results.


### Results

Contains results and temporary files


### data


 - mageck:   Contains MAGeCK result files
 - counts:   Raw counts obtained from the sequencing pipeline described in Bosch &amp; DeJesus et al. 2021.
 - Final_PreDepletion_Grouped_Split.csv: Table summarizing MAGeCK analysis results
 - list_of_H37Rv_orfs_in_library.txt: List of genes in M. Tuberculosis strain H37Rv
 - Weiz_wLFC.csv: Cleaned copy of Supplementary Data from Weizhen et al 2017.


## MAGECK Analysis

All MAGECK analysis were done using this command with the following flags:

    mageck test -k <path to file with combined control and experimental counts> -c 0,1,2 -t 3,4,5 -n <output file name> --gene-lfc-method alphamedian --norm-method control --control-sgrna <file with negative control IDs>


Gene-level result outputs are provided in the repository, as well as files with read-counts for each of the conditions.



## Heatmap

Read through the MAGeCK result files to generate a table (~/Results/Heatmap/htmp_data_gt1Tmt_D5.csv) of log-fold change values as described in the paper.

    ~/python Data_Prep_Heatmap_Script.py

Take the log-fold change data table and create the Heatmap (~/Results/Heatmap/htmp_Ward_D_gt1Tmt_D5.png) from Supplementary material.
   
    ~/Rscript Heatmap_Script_D5.R



## TnSeq Comparison

Read cleaned supplementary data from Weizhen et al 2017 (~/data/Weiz_wLFC.csv) and the table summarizing MAGeCK results (~/data/Final_PreDepletion_Grouped_Split.csv) to generate Supplementary Data 2 (~/Results/Supplemental_Data_2/Supplemental_Data_2.xlsx).

    ~/python Gen_Weizhen_Overlap_Supplemental_Fig2.py
