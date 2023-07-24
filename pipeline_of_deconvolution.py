
from Bio import SeqIO
import pandas as pd
records = []
for seqrecord in SeqIO.parse("all_sequences.fastq", "fastq"):
    #print(seqrecord.seq)
    record = []
    record.append(seqrecord.id)
    record.append(len(seqrecord.seq))
    records.append(record)
all_seq_data = pd.DataFrame(records)
all_seq_data.columns = ['seqID', 'read_length']

index_mismatch0=pd.read_csv("index_result.txt",delimiter="\t") # for index barcode without mismatch
reverse_dataset=pd.read_csv("reverse_result.txt",delimiter="\t") # for barcode reverse 16s, 18s, its 
barcode_dataset=pd.read_csv("forward_result.txt",delimiter="\t") # for barcode forward 16s, 18s, its

index_forward0 = pd.merge(index_mismatch0, barcode_dataset, on ="seqID")

tx_data = index_forward0
primer_filter_1 = tx_data.drop(tx_data[(tx_data['start_x'] >=80) & (tx_data['start_x'] - 1 != tx_data['end_y'])].index)

primer_filter_2 = primer_filter_1.drop(primer_filter_1[(primer_filter_1['start_x'] <80) & (primer_filter_1[ 'start_y'] !=primer_filter_1['end_x']+1)].index)
  
w_primer2 = tx_data.loc[(tx_data['start_x'] >=80) & (tx_data['start_x'] - 2 == tx_data['end_y'])]
w_primer = tx_data.loc[(tx_data['start_x'] <80) & (tx_data['start_y'] == tx_data['end_x']+2)]
liste_filt = [primer_filter_2, w_primer, w_primer2]
index0forward_data = pd.concat(liste_filt)

index0_forward_reverse = pd.merge(index0forward_data, reverse_dataset, on ="seqID")

d_new = index0_forward_reverse.drop(index0_forward_reverse[(index0_forward_reverse ['patternName_y'] == '16S_341F')&(index0_forward_reverse ['patternName']=='18S_1626r')].index)

data_filter_2 = d_new.drop(d_new[(d_new['patternName_y'] == '18S_391f')&(d_new['patternName']=='16S_1391R')].index) 
#all_data_filter_2
data_filter_3 = data_filter_2.drop(data_filter_2[(data_filter_2['patternName_y'] == '16S_341F')&(data_filter_2['patternName']=='ITS_u4')].index) 
#all_data_filter_3
data_filter_4 = data_filter_3.drop(data_filter_3[( data_filter_3['patternName_y'] == '18S_391f')&( data_filter_3['patternName']=='ITS_u4')].index) 
#all_data_filter_4
data_filter_5 = data_filter_4.drop(data_filter_4[(data_filter_4['patternName_y'] == 'ITS_u1')&(data_filter_4['patternName']=='16S_1391R')].index) 
#all_data_filter_5 
data_filter_6 = data_filter_5.drop(data_filter_5[(data_filter_5['patternName_y'] == 'ITS_u1')&(data_filter_5['patternName']=='18S_1626r')].index)
#all_data_filter_6
#data_filter_7 = data_filter_6.drop(data_filter_6[(data_filter_6['patternName_y'] == 'ITS_u1')&(data_filter_6['patternName']=='ITS_u4')&(data_filter_6['read_length']>900)].index)
#all_data_filter_7                                                                                                               
#data_filter_8 = data_filter_7.drop(data_filter_7[(data_filter_7['patternName_y'] == '16S_341F')&(data_filter_7['patternName']=='16S_1391R')&(data_filter_7['read_length']<1000)].index)
#all_data_filter_8
#data_filter_index0 = data_filter_8.drop(data_filter_8[(data_filter_8['patternName_y'] == '18S_391f')&(data_filter_8['patternName']=='18S_1626r')&(data_filter_8['read_length']<1000)].index)  

for group, dataframe in data_filter_6.groupby('patternName_x'):
    
    # save the dataframe for each group to a csv
    dataframe.to_csv(f'{group}.csv', index=False)
