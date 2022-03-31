# CGI\_nature\_micro\_2022
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


### Results

Contains results and temporary files


### data


 - mageck:   Contains MAGeCK result files
 - counts:   Raw counts obtained from the sequencing pipeline described in Bosch &amp; DeJesus et al. 2021.


## MAGECK Analysis

All MAGECK analysis were done using this command with the following flags:

    mageck test -k <path to file with combined control and experimental counts> -c 0,1,2 -t 3,4,5 -n <output file name> --gene-lfc-method alphamedian --norm-method control --control-sgrna <file with negative control IDs>


Gene-level result outputs are provided in the repository, as well as files with read-counts for each of the conditions.



## Heatmap



