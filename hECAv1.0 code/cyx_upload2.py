import sys
import pandas as pd
sys.path.append(sys.argv[1])
from DB_Uploader import DB_Uploader
import time


# load data
df_exp_path = sys.argv[2]
df_anno_path = sys.argv[3]
df_expression = pd.read_csv(df_exp_path,index_col = 0)
df_annotation = pd.read_csv(df_anno_path,index_col = 0)
df_annotation = df_annotation.fillna("NA")

# connected to the server
dbup = DB_Uploader.DB_Uploader()

#Please enter the information of the server here.
endpoint = "https://HCAd-Datasets.cn-beijing.ots.aliyuncs.com"
access_key_id = sys.argv[4]
access_key_secret = sys.argv[5]
instance_name = "HCAd-Datasets"

dbup.Setup_Client(endpoint, access_key_id, access_key_secret, instance_name, 'HCAd')

# Upload matrix
start_row = 0
try_times = 1
while(start_row>=0 and try_times<500)
    print("This is the %dth time, start at the %th row."%(try_times,start_row))
    start_row = -1* dbup.insert_matrix(df_expression,df_annotation, True, start_row)
    t = t+1
    time.sleep(5)
if t == 500 and start_row>=0:
    print("Reach the max retry times,failed at the %th row."%start_row)