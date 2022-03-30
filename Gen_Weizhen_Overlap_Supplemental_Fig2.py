import pandas as pd
import crispr_tools
import ast

WData = pd.read_csv('data/Weiz_wLFC.csv')
WData['Gene_name'] = WData['Gene_name'].replace('rv','RVBD',regex=True)
WData['Full_ID'] = WData[['Gene', 'Gene_name']].agg(':'.join, axis=1)

# Rif Results
W_Rif_Tab = WData[(WData.Rif_q.astype(float) <= 0.01) & (abs(WData.Rif_l2fc) >= 1)]
W_Rif_Tab = W_Rif_Tab.drop(columns = [u'Unnamed: 0', u'Gene', u'Gene_name', u'Vanc_l2fc', u'Vanc_q',
       u'INH_l2fc', u'INH_q', u'Emb_l2fc', u'Emb_q'])
W_Rif_Hits = W_Rif_Tab['Full_ID'].values

# Vanc Results
W_Vanc_Tab = WData[(WData.Vanc_q.astype(float) <= 0.01) & (abs(WData.Vanc_l2fc) >= 1)]
W_Vanc_Tab = W_Vanc_Tab.drop(columns = [u'Unnamed: 0', u'Gene', u'Gene_name',
       u'Rif_l2fc', u'Rif_q', u'INH_l2fc', u'INH_q', u'Emb_l2fc', u'Emb_q'])
W_Vanc_Hits = W_Vanc_Tab['Full_ID'].values

# INH Results
W_INH_Tab = WData[(WData.INH_q.astype(float) <= 0.01) & (abs(WData.INH_l2fc) >= 1)]
W_INH_Tab = W_INH_Tab.drop(columns = [u'Unnamed: 0', u'Gene', u'Gene_name', u'Vanc_l2fc', u'Vanc_q',
       u'Rif_l2fc', u'Rif_q', u'Emb_l2fc', u'Emb_q'])
W_INH_Hits = W_INH_Tab['Full_ID'].values

# EMB Results
W_Emb_Tab = WData[(WData.Emb_q.astype(float) <= 0.01) & (abs(WData.Emb_l2fc) >= 1)]
W_Emb_Tab = W_Emb_Tab.drop(columns=[u'Unnamed: 0', u'Gene', u'Gene_name', u'Vanc_l2fc', u'Vanc_q',
       u'Rif_l2fc', u'Rif_q', u'INH_l2fc', u'INH_q'])
W_Emb_Hits = W_Emb_Tab['Full_ID'].values

######

def InBoth(first, second):
    overlap_genes = []
    for x in first:
        if x in second:
            overlap_genes.append(x)
    return overlap_genes

# import Day1+Day5 data
dat_path = 'data/'
ourData = pd.read_csv(dat_path+'Final_PreDepletion_Grouped_Split.csv', index_col=0)
ourSmall = ourData[ourData.Comparison.isin(['INH_5', 'Emb_1', 'Rif_25','Vanc_25'])]
ourSmall['Essential Depletion Day1+Day5 List'] = ourSmall['Essential Depletion Day1+Day5 List'].map(ast.literal_eval)
ourSmall['Essential Enrichment Day1+Day5 List'] = ourSmall['Essential Enrichment Day1+Day5 List'].map(ast.literal_eval)
ourSmall['NonEssential Depletion Day1+Day5 List'] = ourSmall['NonEssential Depletion Day1+Day5 List'].map(ast.literal_eval)
ourSmall['NonEssential Enrichment Day1+Day5 List'] = ourSmall['NonEssential Enrichment Day1+Day5 List'].map(ast.literal_eval)

# INH CRISPRi Hits
inh_dep_non = ourSmall[ourSmall.Comparison == 'INH_5']['NonEssential Depletion Day1+Day5 List'].values[0]
inh_enr_non = ourSmall[ourSmall.Comparison == 'INH_5']['NonEssential Enrichment Day1+Day5 List'].values[0]
inh_dep_ess = ourSmall[ourSmall.Comparison == 'INH_5']['Essential Depletion Day1+Day5 List'].values[0]
inh_enr_ess = ourSmall[ourSmall.Comparison == 'INH_5']['Essential Enrichment Day1+Day5 List'].values[0]
inh_list = inh_dep_non + inh_enr_non + inh_dep_ess + inh_enr_ess
inh_set = set(inh_list)

# Rif CRISPRi Hits
rif_dep_non = ourSmall[ourSmall.Comparison == 'Rif_25']['NonEssential Depletion Day1+Day5 List'].values[0]
rif_enr_non = ourSmall[ourSmall.Comparison == 'Rif_25']['NonEssential Enrichment Day1+Day5 List'].values[0]
rif_dep_ess = ourSmall[ourSmall.Comparison == 'Rif_25']['Essential Depletion Day1+Day5 List'].values[0]
rif_enr_ess = ourSmall[ourSmall.Comparison == 'Rif_25']['Essential Enrichment Day1+Day5 List'].values[0]
rif_list = rif_dep_non + rif_enr_non + rif_dep_ess + rif_enr_ess
rif_set = set(rif_list)

# Vanc CRISPRi Hits
vanc_dep_non = ourSmall[ourSmall.Comparison == 'Vanc_25']['NonEssential Depletion Day1+Day5 List'].values[0]
vanc_enr_non = ourSmall[ourSmall.Comparison == 'Vanc_25']['NonEssential Enrichment Day1+Day5 List'].values[0]
vanc_dep_ess = ourSmall[ourSmall.Comparison == 'Vanc_25']['Essential Depletion Day1+Day5 List'].values[0]
vanc_enr_ess = ourSmall[ourSmall.Comparison == 'Vanc_25']['Essential Enrichment Day1+Day5 List'].values[0]
vanc_list = vanc_dep_non + vanc_enr_non + vanc_dep_ess + vanc_enr_ess
vanc_set = set(vanc_list)

# Emb CRISPRi Hits
emb_dep_non = ourSmall[ourSmall.Comparison == 'Emb_1']['NonEssential Depletion Day1+Day5 List'].values[0]
emb_enr_non = ourSmall[ourSmall.Comparison == 'Emb_1']['NonEssential Enrichment Day1+Day5 List'].values[0]
emb_dep_ess = ourSmall[ourSmall.Comparison == 'Emb_1']['Essential Depletion Day1+Day5 List'].values[0]
emb_enr_ess = ourSmall[ourSmall.Comparison == 'Emb_1']['Essential Enrichment Day1+Day5 List'].values[0]
emb_list = emb_dep_non + emb_enr_non + emb_dep_ess + emb_enr_ess
emb_set = set(emb_list)

###### OVERLAPS ########
# Rif
overlap_list_rif = InBoth(rif_set, W_Rif_Hits)
# INH
overlap_list_inh = InBoth(inh_set, W_INH_Hits)
# Vanc
overlap_list_vanc = InBoth(vanc_set, W_Vanc_Hits)
# Emb
overlap_list_emb = InBoth(emb_set, W_Emb_Hits)

print('Rif Weiz: ', len(W_Rif_Hits))
print('INH Weiz: ', len(W_INH_Hits))
print('Vanc Weiz: ', len(W_Vanc_Hits))
print('Emb Weiz: ', len(W_Emb_Hits))

print('Rif overlap: ', len(overlap_list_rif))
print('INH overlap: ', len(overlap_list_inh))
print('Vanc overlap: ', len(overlap_list_vanc))
print('Emb overlap: ', len(overlap_list_emb))

# adding CRISPRi_Hit column to Weizhen columns (p, lfc, ID)
def gen_overlap_column(table, olist):
    crisprHit_list = []
    for i in range(len(table)):
        if table.iloc[i][-1] in olist:
            crisprHit_list.append('True')
        else:
            crisprHit_list.append('False')
    table['CRISPRi_Hit'] = crisprHit_list
    return table








Rif_Tab = gen_overlap_column(W_Rif_Tab, overlap_list_rif)
Vanc_Tab = gen_overlap_column(W_Vanc_Tab, overlap_list_vanc)
INH_Tab = gen_overlap_column(W_INH_Tab, overlap_list_inh)
Emb_Tab = gen_overlap_column(W_Emb_Tab, overlap_list_emb)

output_path = 'Results/Supplement_Data_2/'
# Rif_Tab.to_csv(output_path+'Rif.csv')
# Vanc_Tab.to_csv(output_path+'Vanc.csv')
# INH_Tab.to_csv(output_path+'INH.csv')
# Emb_Tab.to_csv(output_path+'Emb.csv')

#emb = pd.read_csv('Supplement/Emb.csv')
emb = Emb_Tab
emb_gene_list = emb['Full_ID'].values
emb_file = 'result_304_EMB_D1_1X_vs_295_RLC0012-1_day-DMSO_alphamedian_control_control_lod100.mageck.gene_summary.txt'
emb_file5 = 'result_315_EMB_D5_1X_vs_308_DMSO_D5_0X_alphamedian_control_control_lod100.mageck.gene_summary.txt'

#inh = pd.read_csv('Supplement/INH.csv')
inh = INH_Tab
inh_gene_list = inh['Full_ID'].values
inh_file = 'result_299_RLC0012-1_day-INH0.5X_vs_295_RLC0012-1_day-DMSO_exp_1_2_alphamedian_control_control_lod100.mageck.gene_summary.txt'
inh_file5 = 'result_483_Pool1_INH_D5_0_5X_vs_308_DMSO_D5_0X_alphamedian_control_control_lod100.mageck.gene_summary.txt'

#rif = pd.read_csv('Supplement/Rif.csv')
rif = Rif_Tab
rif_gene_list = rif['Full_ID'].values
rif_file = 'result_298_RLC0012-1_day-Rif0.25_vs_295_RLC0012-1_day-DMSO_alphamedian_control_control_lod100.mageck.gene_summary.txt'
rif_file5 = 'result_311_RIF_D5_0_25X_vs_308_DMSO_D5_0X_alphamedian_control_control_lod100.mageck.gene_summary.txt'

#vanc = pd.read_csv('Supplement/Vanc.csv')
vanc = Vanc_Tab
vanc_gene_list = vanc['Full_ID'].values
vanc_file = 'result_1970_VAN_0_25X_D1_vs_1962_DMSO_D1_alphamedian_control_control_lod100.mageck.gene_summary.txt'
vanc_file5 = 'result_1975_VAN_0_25X_D5_vs_1972_DMSO_D5_alphamedian_control_control_lod100.mageck.gene_summary.txt'

# take in list of genes (weizhen hits) and pull lfc + fdr from mageck file
def gen_Ci_table(gene_list, treatment_file, prefix):
    AllDrug_Table = pd.DataFrame(data=gene_list, columns = ['Gene_ID'])

    m_path = 'data/mageck/'+treatment_file
    M_df = pd.read_csv(m_path, sep='\t')

    l2fc_list = []
    fdr_list = []

    for gene in gene_list:
        curTab = M_df[M_df.id == gene]
        #print(gene)
        try:
            if curTab['neg|fdr'].values[0] > curTab['pos|fdr'].values[0]:
                l2fc_list.append(curTab['pos|lfc'].values[0])
                fdr_list.append(curTab['pos|fdr'].values[0])
            else:
                l2fc_list.append(curTab['neg|lfc'].values[0])
                fdr_list.append(curTab['neg|fdr'].values[0])
        except IndexError:
            l2fc_list.append('NA')
            fdr_list.append('NA')

    AllDrug_Table[prefix+'CRISPRi_L2FC'] = l2fc_list
    AllDrug_Table[prefix+'CRISPRi_fdr'] = fdr_list

    return AllDrug_Table


## written to output intermediate tables (to check) and then load back in below ##
emb_Ci = gen_Ci_table(emb_gene_list, emb_file, prefix='D1_')
#emb_Ci.to_csv('Supplement/Emb_Ci.csv')
emb_Ci_5 = gen_Ci_table(emb_gene_list, emb_file5, prefix='D5_')
#emb_Ci_5.to_csv('Supplement/Emb_Ci_5.csv')

inh_Ci = gen_Ci_table(inh_gene_list, inh_file, prefix = 'D1_')
#inh_Ci.to_csv('Supplement/INH_Ci.csv')
inh_Ci_5 = gen_Ci_table(inh_gene_list, inh_file5, prefix = 'D5_')
#inh_Ci_5.to_csv('Supplement/INH_Ci_5.csv')

rif_Ci = gen_Ci_table(rif_gene_list, rif_file, prefix = 'D1_')
#rif_Ci.to_csv('Supplement/Rif_Ci.csv')
rif_Ci_5 = gen_Ci_table(rif_gene_list, rif_file5, prefix = 'D5_')
#rif_Ci_5.to_csv('Supplement/Rif_Ci_5.csv')

vanc_Ci = gen_Ci_table(vanc_gene_list, vanc_file, prefix = 'D1_')
vanc_Ci.to_csv('Supplement/Vanc_Ci.csv')
vanc_Ci_5 = gen_Ci_table(vanc_gene_list, vanc_file5, prefix = 'D5_')
vanc_Ci_5.to_csv('Supplement/Vanc_Ci_5.csv')



output_path = 'Results/Supplemental_Data_2/'
#### Have all 3 tables - Weizhen+Hit_Status, Mag_1 + Mag_5

###### Emb
#emb_weiz = pd.read_csv(output_path + 'Emb.csv', index_col=0)
emb_weiz = Emb_Tab
emb_weiz = emb_weiz.rename(columns={'Full_ID':'Gene_ID'})
#emb_mag_1 = pd.read_csv(output_path + 'Emb_Ci.csv', index_col=0)
emb_mag_1 = emb_Ci
#emb_mag_5 = pd.read_csv(output_path + 'Emb_Ci_5.csv', index_col=0)
emb_mag_5 = emb_Ci_5
emb_1 = pd.merge(emb_weiz, emb_mag_1, on='Gene_ID')
emb_15 = pd.merge(emb_1, emb_mag_5, on='Gene_ID')
emb_final = emb_15[['Gene_ID', 'Emb_l2fc', 'Emb_q', 'CRISPRi_Hit', 'D1_CRISPRi_L2FC', 'D1_CRISPRi_fdr', 'D5_CRISPRi_L2FC', 'D5_CRISPRi_fdr']]
emb_final.columns = ['ORF_ID', 'TnSeq_l2fc', 'TnSeq_p_adj', 'CRISPRi_Hit', 'CRISPRi_D1_l2fc', 'CRISPRi_D1_p_adj', 'CRISPRi_D5_l2fc', 'CRISPRi_D5_p_adj']
emb_final = emb_final.set_index('ORF_ID')

####### inh
#inh_weiz = pd.read_csv(output_path + 'INH.csv', index_col=0)
inh_weiz = Inh_Tab
inh_weiz = inh_weiz.rename(columns={'Full_ID':'Gene_ID'})
#inh_mag_1 = pd.read_csv(output_path + 'INH_Ci.csv', index_col=0)
inh_mag_1 = inh_Ci
#inh_mag_5 = pd.read_csv(output_path + 'INH_Ci_5.csv', index_col=0)
inh_mag_5 = inh_Ci_5
inh_1 = pd.merge(inh_weiz, inh_mag_1, on='Gene_ID')
inh_15 = pd.merge(inh_1, inh_mag_5, on='Gene_ID')
inh_final = inh_15[['Gene_ID', 'INH_l2fc', 'INH_q', 'CRISPRi_Hit', 'D1_CRISPRi_L2FC', 'D1_CRISPRi_fdr', 'D5_CRISPRi_L2FC', 'D5_CRISPRi_fdr']]
inh_final.columns = ['ORF_ID', 'TnSeq_l2fc', 'TnSeq_p_adj', 'CRISPRi_Hit', 'CRISPRi_D1_l2fc', 'CRISPRi_D1_p_adj', 'CRISPRi_D5_l2fc', 'CRISPRi_D5_p_adj']
inh_final = inh_final.set_index('ORF_ID')

####### rif
#rif_weiz = pd.read_csv(output_path + 'Rif.csv', index_col=0)
rif_weiz = Rif_Tab
rif_weiz = rif_weiz.rename(columns={'Full_ID':'Gene_ID'})
#rif_mag_1 = pd.read_csv(output_path + 'Rif_Ci.csv', index_col=0)
rif_mag_1 = rif_Ci
#rif_mag_5 = pd.read_csv(output_path + 'Rif_Ci_5.csv', index_col=0)
rif_mag_5 = rif_Ci_5
rif_1 = pd.merge(rif_weiz, rif_mag_1, on='Gene_ID')
rif_15 = pd.merge(rif_1, rif_mag_5, on='Gene_ID')
rif_final = rif_15[['Gene_ID', 'Rif_l2fc', 'Rif_q', 'CRISPRi_Hit', 'D1_CRISPRi_L2FC', 'D1_CRISPRi_fdr', 'D5_CRISPRi_L2FC', 'D5_CRISPRi_fdr']]
rif_final.columns = ['ORF_ID', 'TnSeq_l2fc', 'TnSeq_p_adj', 'CRISPRi_Hit', 'CRISPRi_D1_l2fc', 'CRISPRi_D1_p_adj', 'CRISPRi_D5_l2fc', 'CRISPRi_D5_p_adj']
rif_final = rif_final.set_index('ORF_ID')

####### vanc
#vanc_weiz = pd.read_csv(output_path + 'Vanc.csv', index_col=0)
vanc_weiz = Vanc_Tab
vanc_weiz = vanc_weiz.rename(columns={'Full_ID':'Gene_ID'})
#vanc_mag_1 = pd.read_csv(output_path + 'Vanc_Ci.csv', index_col=0)
vanc_mag_1 = vanc_Ci
#vanc_mag_5 = pd.read_csv(output_path + 'Vanc_Ci_5.csv', index_col=0)
vanc_mag_5 = vanc_Ci_5
vanc_1 = pd.merge(vanc_weiz, vanc_mag_1, on='Gene_ID')
vanc_15 = pd.merge(vanc_1, vanc_mag_5, on='Gene_ID')
vanc_final = vanc_15[['Gene_ID', 'Vanc_l2fc', 'Vanc_q', 'CRISPRi_Hit', 'D1_CRISPRi_L2FC', 'D1_CRISPRi_fdr', 'D5_CRISPRi_L2FC', 'D5_CRISPRi_fdr']]
vanc_final.columns = ['ORF_ID', 'TnSeq_l2fc', 'TnSeq_p_adj', 'CRISPRi_Hit', 'CRISPRi_D1_l2fc', 'CRISPRi_D1_p_adj', 'CRISPRi_D5_l2fc', 'CRISPRi_D5_p_adj']
vanc_final = vanc_final.set_index('ORF_ID')

writer = pd.ExcelWriter(output_path+'Supplemental_Data_2.xlsx')
emb_final.to_excel(writer, sheet_name='EMB')
inh_final.to_excel(writer, sheet_name='INH')
rif_final.to_excel(writer, sheet_name='RIF')
vanc_final.to_excel(writer, sheet_name='VAN')
writer.save()
writer.close()

# add legend tab
# only left would be switch to arial and column width **** AND THE TRUE/FALSE to Yes/No




