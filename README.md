# CGI_nature_micro_2022
Code used in the Shuqi &amp; Poulton et al. Nature Micro 2022 paper



## Requirements


### Python code:

 - MAGeCK:  https://sourceforge.net/p/mageck/wiki/install/
 - pandas
 - numpy
 - scipy
 - matplotlib
 - seaborn


### R code:

 - ggplot2



## Directory Structure:


### ./ 

Contains the scripts used to analyze counts and plot results.


### data


 - mageck:   Contains MAGeCK result files
 - counts:   Raw counts obtained from the sequencing pipeline described in Bosch &amp; DeJesus et al. 2021.


## Analysis

    mageck test -k temp_combined_mageck_331_LEVO_D10_0_125X_vs_323_DMSO_D10_0X_alphamedian_control_control_lod100.txt -c 0,1,2 -t 3,4,5 -n result_331_LEVO_D10_0_125X_vs_323_DMSO_D10_0X_alphamedian_control_control_lod100.mageck --gene-lfc-method alphamedian --norm-method control --control-sgrna negatives_mageck_331_LEVO_D10_0_125X_vs_323_DMSO_D10_0X_alphamedian_control_control_lod100.txt


## Heatmap



