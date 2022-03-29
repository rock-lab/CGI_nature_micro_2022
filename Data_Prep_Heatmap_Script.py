import pandas as pd
#import crispr_tools
import numpy as np

def gen_lfc_table(comparison_list):
    sample_Res_file = pd.read_csv('result_296_RLC0012-1_day-Rif0.0625_vs_295_RLC0012-1_day-DMSO_lod20.resampling copy.txt', sep='\t')
    all_gene_list = sample_Res_file.gene.values
    AllDrug_Table = pd.DataFrame(data=all_gene_list, columns = ['id'])

    for comp in comparison_list:
        print(comp)
        m_path = '/Users/zacharyazadian/RockLab Dropbox/Projects/chemical_genomics/data/MAGECK/result_'+comp+'_alphamedian_control_control_lod100.mageck.gene_summary.txt'
        #M = crispr_tools.data.results.open_file(m_path)
        #M_df = M.data
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

path = '/Users/zacharyazadian/RockLab Dropbox/Zachary Azadian/ChemGen/New_Data_Results/AAA_Final_Data/Heatmap/'

############# SETTING WHAT TREATMENTS TO LOOK AT ####################
# Comparison_DF = pd.read_csv(path+'D5_List.csv')
# Comparison_List = Comparison_DF['0'].values
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
path = '/Users/zacharyazadian/RockLab Dropbox/Zachary Azadian/ChemGen/New_Data_Results/AAA_Final_Data/Heatmap/'
lfc_table = gen_lfc_table(Comparison_List)
#direction_lfc_table.to_csv(path+'LFC_table_D5.csv')

def gen_fdr_table(comparison_list):

    sample_Res_file = pd.read_csv('result_296_RLC0012-1_day-Rif0.0625_vs_295_RLC0012-1_day-DMSO_lod20.resampling copy.txt', sep='\t')
    all_gene_list = sample_Res_file.gene.values
    AllDrug_Table = pd.DataFrame(data=all_gene_list, columns = ['id'])

    for comp in comparison_list:
        print(comp)
        m_path = '/Users/zacharyazadian/RockLab Dropbox/Projects/chemical_genomics/data/MAGECK/result_'+comp+'_alphamedian_control_control_lod100.mageck.gene_summary.txt'
        #M = crispr_tools.data.results.open_file(m_path)
        #M_df = M.data
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
#fdr_table.to_csv(path+'FDR_table_D5.csv')

############# LOAD IN CORRECT TABLES ####################
#fdr_table = pd.read_csv('FDR_table_D5.csv', index_col=0)
fdr_table = fdr_table.fillna(0)

# making array of fdr values then boolean for which ones are hits
fdr_array = np.array(fdr_table.drop(columns='id').values)
fdr_hits_array = fdr_array < 0.01

############# NOW LFC ####################
#lfc_table = pd.read_csv('LFC_table_D5.csv', index_col=0)
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

either_freq_vector = dep_freq_vector + enr_freq_vector
bool_either_freq_vector_gt1 = either_freq_vector > 1

# print('This many genes were a hit (enr or depl) in more than one treatment',
#       np.sum(bool_either_freq_vector_gt1))
#
# print('How many deplete in more than one tmt:',np.sum(bool_dep_freq_vector_gt1))
# print('How many enrich in more than one tmt:',np.sum(bool_enr_freq_vector_gt1))
# print('There could be a discrepancy for cases where a gene depletes in one tmt and enriches in one tmt. That gene would'
#       'only begin to make the cut when enr and dep vectors are added together')

# Need to give heatmap lfc table of only the 'or' genes, also replace Nan with zero
lfc_table = lfc_table.set_index('id')
or_table = lfc_table[bool_either_freq_vector_gt1]

############# LOAD IN CORRECT FDR TABLE ####################
#fdr_table = pd.read_csv('FDR_table_D5.csv', index_col=0)
fdr_table = gen_fdr_table(Comparison_List)
fdr_table = fdr_table.set_index('id')

fdr_hits_table = fdr_table[bool_either_freq_vector_gt1]
fdr_hits_bool = fdr_hits_table > 0.01 # this boolean will help assign zero or not depending on if gene is hit in that condition

or_table[fdr_hits_bool] = 0
or_table = or_table.fillna(0)

############# SAVE OUTPUT USED FOR HEATMAP #############
#or_table.to_csv('htmp_data_gt1Tmt_D5.csv')


