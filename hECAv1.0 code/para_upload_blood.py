import sys
import numpy as np
import pandas as pd
import tablestore
import time
import os
from queue import Queue
import multiprocessing
import warnings
warnings.filterwarnings('ignore')

sys.path.append(sys.argv[1])

# Please enter the information of the server here.
endpoint = "https://HCAd-Datasets3.cn-beijing.ots.aliyuncs.com"
access_key_id = sys.argv[3]
access_key_secret = sys.argv[4]
instance_name = "HCAd-Datasets3"
table_name = "HCA_d"

# Connect to the server
ots_client = tablestore.OTSClient(endpoint, access_key_id, access_key_secret, instance_name, logger_name = "table_store.log", retry_policy = tablestore.WriteRetryPolicy())
threads_num = int(sys.argv[5])
if threads_num < 0:
    threads_num = multiprocessing.cpu_count() - 1
threads_num = min(multiprocessing.cpu_count() - 1, threads_num)


def Cell2Row(sample):
    # build the attribute columns part
    attribute_columns = []
    for i in range(sample.shape[0]):
        if sample.index[i]!= "cid":
            if sample.index[i] == 'user_id':
                attribute_columns.append((sample.index[i],int(sample[i])))
            elif isinstance(sample[i],np.int64):
                attribute_columns.append((sample.index[i],float(sample[i])))
            else:
                attribute_columns.append((sample.index[i],sample[i]))
    
    # build the primary_key part 
    primary_key = [('cid',int(sample['cid']))]
    
    # the maxium number of attribute columns in each writing operation is 1024, so we split a row into blocks
    row_blocks = []
    for i in range(len(attribute_columns)//1024+1): # the maxium number of attribute columns in each writing operation is 1024
        if i==0:
            row_blocks.append(tablestore.Row(primary_key,attribute_columns[i*1024:min(i*1024+1024,len(attribute_columns))]))
        else:
            row_blocks.append(tablestore.Row(primary_key,{'PUT':attribute_columns[i*1024:min(i*1024+1024,len(attribute_columns))]}))
            
    return row_blocks


def insert_row(row):
    condition = tablestore.Condition(tablestore.RowExistenceExpectation.IGNORE)
    
    for i in range(len(row)):
        try :
            if (i==0):
                consumed, return_row = ots_client.put_row(table_name, row[i], condition)
            else:
                consumed, return_row = ots_client.update_row(table_name, row[i], condition)
        except tablestore.OTSClientError as e:
            print (e.get_error_message())
            return(-1)
        except tablestore.OTSServiceError as e:
            print (e.get_error_message())
            return(-1)
    return(0)


def para_upload(bag_arg):
    i = bag_arg
    row = Cell2Row(df.iloc[i,:])
    insert_stat = insert_row(row)
    return(insert_stat)

def check_df(genenum_chk=True):
    if genenum_chk and df.shape[1] != 43896:
        print("Gene number error.")
        return False
    
    if sum([x in df.columns for x in ['user_id', 'study_id', 'cell_id', 'organ', 'region', 'subregion','seq_tech', 'sample_status', 'donor_id', 'donor_gender', 'donor_age', 'original_name', 'cl_name', 'hcad_name', 'tissue_type', 'cell_type', 'marker_gene','cid']]) != 18:
        print("Metadata names error.")
        return False
    
    if set(np.unique(df['donor_gender'])).union({'Female', 'Male', 'NA'}) != {'Female', 'Male', 'NA'}:
        print("donor_gender error.")
        return False

    if set(np.unique(df['donor_gender'])).union({'Female', 'Male', 'NA'}) != {'Female', 'Male', 'NA'}:
        print("donor_gender error.")
        return False

    quota = df['user_id'][0]
    if min(df['cid'])<quota*1000000 or max(df['cid'])>quota*1000000 + 999999:
        print("cid error.")
        return False
    
    print("Dataset pass the check.")
    return True


# Load data
df_path = sys.argv[2]
df = pd.read_csv(df_path,index_col = 0)
for i in range(20):
    df.iloc[:,-i][df.iloc[:,-i]=="To be Change"]="NA"
if check_df()==False:
    sys.exit(0)

# Parallel Upload
print("Upload start.")
start_time = time.time()

nrow = df.shape[0]
indexes = range(nrow)
arg_bag = [(i) for i in indexes]
with multiprocessing.Pool(threads_num) as pool:
    b = pool.map(para_upload,arg_bag)
    pool.close()
    pool.join()
indexes = np.array(indexes)[np.array(b)==-1]

while sum(b)!= 0:
    arg_bag = [(start_row,i) for i in indexes]
    with multiprocessing.Pool(threads_num) as pool:
        b = pool.map(para_upload,arg_bag)
        pool.close()
        pool.join()
    indexes = np.array(indexes)[np.array(b)==-1]

print(nrow, "cells uploaded.",time.time()-start_time,"seconds used.")
print("next cid should be:",max(df['cid'])+1)
