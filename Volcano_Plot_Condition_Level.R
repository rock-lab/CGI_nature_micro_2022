#######
BiocManager::install("MAGeCKFlute")
########
library(MAGeCKFlute)
library(ggplot2)
library(stringr)
library(parallel)

# ggplot theme Michael created
theme_set(theme_bw()) 
theme_update(axis.text = element_text(size=15), title = element_text(size=18),
             axis.title.y  = element_text(size=25, margin = margin(t = 0, r = 10, b = 0, l = 0)),
             axis.title.x  = element_text(size=25, margin = margin(t = 10, r = 0, b = 0, l = 0)),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(), 
             panel.border = element_rect(size=2),
             legend.text = element_text(size=15))

# MAGECK result files
read_path = 'data/mageck/'
prefix = 'result_'
suffix = '_alphamedian_control_control_lod100.mageck.gene_summary.txt'

# functions # 
make_plot <- function(gene, rra, Xlimit, Xtick, Drug, Day, Conc){
  print(paste('Generating a plot for condition', Drug, Day, Conc, sep=' '))
  cols <- c("0" = "gray", "1" = "red", "2" = "blue")
  border_cols <- c("0" = "gray", "1" = "red", "2" = "blue")

  rra$Hit_Status <- factor(rra$Hit_Status)
  gene <- str_replace(gene, ':', '_')
  
  # past together inputs from dictionary to finalize this path
  output_path_prefix = 'Results/Volcano_Plots/'
  output_path_full = paste(output_path_prefix, paste(Drug, Day, Conc, sep='_'), sep='')
  print(output_path_full)
  
  # .png would cause problems in path so save it for png()
  png(paste(output_path_full, '.png',sep=''), width=9, height=10, units="in", res=300)
  
  fig <- ggplot(rra, aes(x = L2FC, y = LogFDR, color=Hit_Status, size = Hit_Status,
                        fill = Hit_Status)) +
    
    geom_point(alpha = 0.5, show.legend=FALSE, shape=21) + 
    
    scale_colour_manual(values = border_cols) + 
    scale_fill_manual(values = cols) + 
    
    scale_x_continuous(limits = c(-Xlimit, Xlimit), breaks=seq(-Xlimit, Xlimit, Xtick)) +
    
    labs(y="-Log10FDR", x = "L2FC") +
    geom_hline(yintercept=-log10(0.01), linetype= "dashed") +
    geom_vline(xintercept=1, linetype="dashed") +
    geom_vline(xintercept=-1, linetype="dashed") +
    scale_size_discrete(range = c(2, 4)) +
    ggtitle(paste(Drug, Day, Conc, sep='_'))
  
  print(fig)
  dev.off() 
}

gen_volcano_plots <- function(treatment, Xlimit, Xtick, Drug, Day, Conc){
  
  # read in treatment specific mageck file
  filename = paste(prefix, treatment, suffix, sep='')
  data = read.table(paste(read_path, filename, sep=''), header = TRUE)
  
  # run Mageck flute and assign hit + size status to genes
  rra = ReadRRA(data)
  rra$LogFDR = -log10(rra$FDR)
  rra$L2FC <- rra$Score
  e_vec = c()
  for(i in 1:nrow(rra)) {
    
    fdr_val <- rra$LogFDR[i]
    l2fc_val <- rra$L2FC[i]
    
    if (fdr_val >= -log10(0.01) & l2fc_val < -1){
      e_vec <- c(e_vec, 2)
    } else if (fdr_val >= -log10(0.01) & l2fc_val > 1){
      e_vec <- c(e_vec, 1)
    } else {
      e_vec <- c(e_vec, 0)
    }
  }
  rra$Hit_Status <- e_vec
  rra$Hit_Status <- factor(rra$Hit_Status)

  gene_list <- rra$id
  ################################### subset of genes
  test_list <- gene_list[1]
  ################################### can just use gene_list when running final
  
  mclapply(test_list, make_plot, mc.cores=3, rra=rra, Xlimit=Xlimit, Xtick=Xtick, Drug=Drug, Day=Day, Conc=Conc)
  }

### MAKE INPUT TABLE ###
input_tab_path = 'data/Volcano_Plot_Axes_Table.csv'
input_tab = read.table(input_tab_path, header = TRUE ,sep = ',')

### ready to go

for (row in 1:nrow(input_tab)){

  gen_volcano_plots(
    treatment = toString(input_tab$File[row]),
    Xlimit = as.double(input_tab$Xmax[row]),
    Xtick = as.double(input_tab$Xtick[row]),
    Drug = toString(input_tab$Drug[row]),
    Day = toString(input_tab$Day[row]),
    Conc = toString(input_tab$Conc[row]) )

}





