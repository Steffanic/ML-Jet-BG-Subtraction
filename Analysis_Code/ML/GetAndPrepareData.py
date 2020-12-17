import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split



def get_data(filename, rows_=10000000000000, names_=['p_T', 'Eta', 'Phi', 'Area', 'Rho', 'p_T-corr', 'N-Trk', 'Angularity', 'Angularity-NW', 'Mean-p_T', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5','distmatch','XMatch', 'X_tru']):
       '''
       Takes a filename and the names of the columns in the CSV file. Returns a DataFrame.
       '''
       dat = pd.read_csv(filename, names=names_,nrows=rows_, low_memory=False)
       dat = dat[dat.p_T!='p_T']
       dat = dat.rename(columns={'Rho': 'Epsilon'})

       return dat

def add_feats(dat):
       dat_copy = dat.copy()
       dat_copy["pythia-mom"] = dat_copy['X_tru']*dat_copy['N-Trk']
       return dat_copy

def label_data(dat, label_metric = 'pythia-mom',label_bins=[-1,0.00000001 ,10,10000000], labels_=[1,2,3]):
       dat['Label'] = pd.cut(dat[label_metric], bins=label_bins, labels=labels_)
       
       dat['Label']=dat['Label'].astype(np.int64)
       
       return dat

def type_data(dat):
       dat_typed = dat.astype({'p_T': np.float64, 'Eta': np.float64, 'Phi': np.float64, 'Area': np.float64, 'Epsilon':np.float64, 'p_T-corr': np.float64, 'N-Trk': np.int64, 'Angularity': np.float64, 'Angularity-NW':np.float64, 'Mean-p_T': np.float64, 'p_T_1': np.float64, 'p_T_2': np.float64, 'p_T_3': np.float64, 'p_T_4': np.float64, 'p_T_5': np.float64, 'distmatch':np.float64, 'XMatch':np.float64, 'X_tru': np.float64})
       return dat_typed
       
def split_data(dat):
       train, test = train_test_split(dat, test_size=0.2, random_state=42)
       return train, test

def drop_feat(dat, feat=['p_T_3', 'p_T_4', 'p_T_5']):
       dat_drop = dat.drop(feat, 1)
       return dat_drop

def balance_classes(dat, label='Label'):
       dat_bal = dat.loc[dat['Label']==2].copy()
       dat_bal = dat_bal.append(dat.loc[dat['Label']==1].sample(n=len(dat.loc[dat['Label']==2]), replace=True, random_state=42))
       dat_bal = dat_bal.append(dat.loc[dat['Label']==0].sample(n=len(dat.loc[dat['Label']==2]), replace=True, random_state=42))
       return dat_bal
