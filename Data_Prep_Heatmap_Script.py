import pandas as pd

import numpy as np

def gen_lfc_table(comparison_list):
    sample_Res_file = pd.read_csv('data/list_of_H37Rv_orfs_in_library.txt', names=['Orf'], header=None)
    all_gene_list = sample_Res_file.Orf.values

    AllDrug_Table = pd.DataFrame(data=all_gene_list, columns = ['id'])

    for comp in comparison_list:
        print(comp)

        m_path = 'data/mageck/result_'+comp+'_alphamedian_control_control_lod100.mageck.gene_summary.txt'

        M_df = pd.read_csv(m_path, sep='\t')

        curTable = pd.DataFrame(index = M_df['id'])
        curCol = []
        for i in range(len(M_df)):
            if M_df['neg|fdr'][i] > M_df['pos|fdr'][i]:
                curCol.append(M_df['pos|lfc'][i])
            else:
                curCol.append(M_df['neg|lfc'][i])
        comp_name = comp.split('_vs')[0]
        curTable[comp_name] = curCol

        AllDrug_Table = AllDrug_Table.merge(curTable, on='id', how='left')

    return AllDrug_Table


path = 'Results/Heatmap/'

############# SETTING WHAT TREATMENTS TO LOOK AT ####################
Comparison_List =['1794_pool_BDQ_0625_5day_vs_308_DMSO_D5_0X',
'1981_CLR_4X_D5_vs_1972_DMSO_D5',
'662_Pool4_0_125X_INH_Day5_vs_308_DMSO_D5_0X',
'310_RIF_D5_0_125X_vs_308_DMSO_D5_0X',
'319_BDQ_D5_0_125X_vs_308_DMSO_D5_0X',
'1979_CLR_1X_D5_vs_1972_DMSO_D5',
'315_EMB_D5_1X_vs_308_DMSO_D5_0X',
'1980_CLR_2X_D5_vs_1972_DMSO_D5',
'311_RIF_D5_0_25X_vs_308_DMSO_D5_0X',
'663_Pool5_0_25X_INH_Day5_vs_308_DMSO_D5_0X',
'313_EMB_D5_0_25X_vs_308_DMSO_D5_0X',
'316_LEVO_D5_0_125X_vs_308_DMSO_D5_0X',
'428_Pool3_Strept_D5_1X_vs_363_Pool6_0X_DMSO_Day5',
'368_Pool8_0_25X_Linez_Day5_vs_363_Pool6_0X_DMSO_Day5',
'483_Pool1_INH_D5_0_5X_vs_308_DMSO_D5_0X',
'366_Pool7_0_0625X_Linez_Day5_vs_363_Pool6_0X_DMSO_Day5',
'321_BDQ_D5_0_25X_vs_308_DMSO_D5_0X',
'318_LEVO_D5_0_5X_vs_308_DMSO_D5_0X',
'1975_VAN_0_25X_D5_vs_1972_DMSO_D5',
'1974_VAN_0_125X_D5_vs_1972_DMSO_D5',
'426_Pool2_Strept_D5_0_25X_vs_363_Pool6_0X_DMSO_Day5',
'1973_VAN_0_0625X_D5_vs_1972_DMSO_D5',
'314_EMB_D5_0_5X_vs_308_DMSO_D5_0X',
'317_LEVO_D5_0_25X_vs_308_DMSO_D5_0X',
'427_Pool2_Strept_D5_0_5X_vs_363_Pool6_0X_DMSO_Day5',
'309_RIF_D5_0_0625X_vs_308_DMSO_D5_0X',
'367_Pool8_0_125X_Linez_Day5_vs_363_Pool6_0X_DMSO_Day5']

############# MAKE APPROPRIATELY LABELLED LFC TABLE ####################

lfc_table = gen_lfc_table(Comparison_List)

def gen_fdr_table(comparison_list):

    sample_Res_file = pd.read_csv('data/list_of_H37Rv_orfs_in_library.txt', names=['Orf'], header=None)
    all_gene_list = sample_Res_file.Orf.values

    AllDrug_Table = pd.DataFrame(data=all_gene_list, columns = ['id'])

    for comp in comparison_list:
        print(comp)

        m_path = 'data/mageck/result_'+comp+'_alphamedian_control_control_lod100.mageck.gene_summary.txt'

        M_df = pd.read_csv(m_path, sep='\t')

        curTable = pd.DataFrame(index = M_df['id'])
        curCol = []
        for i in range(len(M_df)):
            if M_df['neg|fdr'][i] > M_df['pos|fdr'][i]:
                curCol.append(M_df['pos|fdr'][i])
            else:
                curCol.append(M_df['neg|fdr'][i])

        comp_name = comp.split('_vs')[0]
        curTable[comp_name] = curCol

        AllDrug_Table = AllDrug_Table.merge(curTable, on='id', how='left')

    return AllDrug_Table

############# MAKING APPROPRIATELY LABELLED FDR TABLE ####################
fdr_table = gen_fdr_table(Comparison_List)

fdr_table = fdr_table.fillna(0)

# making array of fdr values then boolean for which ones are hits
fdr_array = np.array(fdr_table.drop(columns='id').values)
fdr_hits_array = fdr_array < 0.01

############# NOW LFC ####################

lfc_table = lfc_table.fillna(0)
lfc_array = np.array(lfc_table.drop(columns='id').values)

# two separate boolean arrays, one for depletion and one for enrichment
lfc_depl_hits_array = lfc_array < -1
lfc_enrich_hits_array = lfc_array > 1

# now multiply fdr and respective lfc arrays to get boolean hit array
Depletion_hits_array = fdr_hits_array * lfc_depl_hits_array
Enrichment_hits_array = fdr_hits_array * lfc_enrich_hits_array

##### Find out how many treatments a gene is a hit under (in either enr or depletion) #####
dep_freq_vector = Depletion_hits_array.sum(axis=1)
bool_dep_freq_vector_gt1 = dep_freq_vector > 1

# do same for enrichment
enr_freq_vector = Enrichment_hits_array.sum(axis=1)
bool_enr_freq_vector_gt1 = enr_freq_vector > 1

# boolean vector to only include genes that are hits in 2 or more tmts
either_freq_vector = dep_freq_vector + enr_freq_vector
bool_either_freq_vector_gt1 = either_freq_vector > 1

# Need to give heatmap lfc table of only the 'or' genes, also replace Nan with zero
lfc_table = lfc_table.set_index('id')
or_table = lfc_table[bool_either_freq_vector_gt1]


############# REMAKE FDR TABLE ####################

fdr_table = gen_fdr_table(Comparison_List)
fdr_table = fdr_table.set_index('id')

fdr_hits_table = fdr_table[bool_either_freq_vector_gt1]
fdr_hits_bool = fdr_hits_table > 0.01 # this boolean will help assign zero or not depending on if gene is hit in that condition

or_table[fdr_hits_bool] = 0
or_table = or_table.fillna(0)

############# SAVE OUTPUT USED FOR HEATMAP #############

or_table.to_csv('Results/Heatmap/htmp_data_gt1Tmt_D5.csv')



