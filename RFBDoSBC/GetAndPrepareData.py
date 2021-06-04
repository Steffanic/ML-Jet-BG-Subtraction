import  pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split

def DataPipeline(filepath, ptm, rows=10000000000, keep_pT=False):
       '''
       Get and Prepare the data from `filepath` 
       '''
       dat = get_data(filepath, rows) #Get data from `filepath`; also cleans repeated headers and renames Eps to Epsilon
       dat_typed = type_data(dat) #Types every column to either an integer or float
       dat_feat_added = add_feats(dat_typed) 
       dat_labeled = label_data(dat_feat_added, 'pythia-mom', [-1,0.000000001, ptm, 10000000], [1,2,3], ["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"])
       dat_drop = drop_feat(dat_labeled, ['Eta', 'Phi', 'p_T', 'Angularity-NW', 'N-Trk', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5', 'distmatch', 'XMatch', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd'] \
              if not keep_pT else ['Eta', 'Phi', 'Angularity-NW', 'N-Trk', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5', 'distmatch', 'XMatch', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd'])
       dat_bal = balance_classes(dat_drop)
       train, test = split_data(dat_bal)
       return train, test

def DataPipelineBatch(filepath, ptm, keep_pT=False, keep_all=False):
       '''
       Get and Prepare a batch of the data from `filepath` 
       '''
       dat_it = get_data_from_iterator(filepath) #Get data from `filepath`; also cleans repeated headers and renames Eps to Epsilon
       for data_batch in dat_it:
              dat = data_batch.drop(data_batch[data_batch.Eps==' Eps'].index)
              dat = dat.rename(columns={'Eps': 'Epsilon'})     
              dat_typed = type_data(dat) #Types every column to either an integer or float
              dat_feat_added = add_feats(dat_typed) 
              dat_labeled = label_data(dat_feat_added, 'pythia-mom', [-1,0.000000001, ptm, 10000000], [1,2,3], ["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"])
              drop_list = ['Eta', 'Phi', 'p_T', 'Angularity-NW', 'N-Trk', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5', 'distmatch', 'XMatch', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd'] \
                     if not keep_pT else ['Eta', 'Phi', 'Angularity-NW', 'N-Trk', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5', 'distmatch', 'XMatch', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd'] \
                            if not keep_all else []
              dat_drop = drop_feat(dat_labeled, drop_list)
              dat_bal = balance_classes(dat_drop)
              train, test = split_data(dat_bal)
              yield train, test

def get_data(filename, rows_=10000000000, names_=['p_T', 'Eta', 'Phi', 'Area', 'Eps', 'p_T-corr', 'N-Trk', 'Angularity', 'Angularity-NW', 'Mean-p_T', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5','distmatch','XMatch', 'X_tru', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd']):
       '''
       Takes a filename and the names of the columns in the CSV file. Returns a DataFrame.
       '''
       dat = pd.read_csv(filename, names=names_,nrows=rows_, low_memory=False)
       dat = dat.drop(dat[dat.Eps==' Eps'].index)
       dat = dat.rename(columns={'Eps': 'Epsilon'})
       return dat

def get_data_from_iterator(filename, chunksize_=80000, names_=['p_T', 'Eta', 'Phi', 'Area', 'Eps', 'p_T-corr', 'N-Trk', 'Angularity', 'Angularity-NW', 'Mean-p_T', 'p_T_1', 'p_T_2', 'p_T_3', 'p_T_4', 'p_T_5','distmatch','XMatch', 'X_tru', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd']):
       '''
       Takes a filename and the names of the columns in the CSV file. Returns a DataFrame.
       '''
       dat = pd.read_csv(filename, names=names_,chunksize=chunksize_, low_memory=False)
       return dat

def add_feats(dat):
       dat_copy = dat.copy()
       dat_copy["pythia-mom"] = dat_copy['X_tru']*dat_copy['N-Trk']
       return dat_copy

def label_data(dat, label_metric = 'pythia-mom',label_bins=[-1,0.00000001 ,10,10000000], labels_=[1,2,3], label_names_=["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"]):
       dat['Label'] = pd.cut(dat[label_metric], bins=label_bins, labels=labels_)
       
       dat['Label']=dat['Label'].astype(np.int64)

       dat['Label_Name'] = [label_names_[i] for i in dat['Label'].apply(lambda x:x-1)]
       
       return dat


def type_data(dat):
       dat_typed = dat.astype({'p_T': np.float32, 'Eta': np.float32, 'Phi': np.float32, 'Area': np.float32, 'Epsilon':np.float32, 'p_T-corr': np.float32, 'N-Trk': np.uint32, 'Angularity': np.float32, 'Angularity-NW':np.float32, 'Mean-p_T': np.float32, 'p_T_1': np.float32, 'p_T_2': np.float32, 'p_T_3': np.float32, 'p_T_4': np.float32, 'p_T_5': np.float32, 'distmatch':np.float32, 'XMatch':np.float32, 'X_tru': np.float32, 'Y_quark': np.float32, 'Y_gluon': np.float32, 'Y_beam':np.float32, 'Y_bkgd':np.float32})
       return dat_typed
       
def split_data(dat):
       train, test = train_test_split(dat, test_size=0.2, random_state=42)
       return train, test

def drop_feat(dat, feat=['p_T_3', 'p_T_4', 'p_T_5', 'Y_quark', 'Y_gluon', 'Y_beam', 'Y_bkgd']):
       dat_drop = dat.drop(feat, 1)
       return dat_drop

def balance_classes(dat, label='Label'):
       dat_bal = dat.loc[dat['Label']==2].copy()
       dat_bal = dat_bal.append(dat.loc[dat['Label']==1].sample(n=len(dat.loc[dat['Label']==2]), replace=True, random_state=42))
       dat_bal = dat_bal.append(dat.loc[dat['Label']==3].sample(n=len(dat.loc[dat['Label']==2]), replace=True, random_state=42))
       return dat_bal

def split_feat_label(data, drop_labels=['p_T', 'X_tru', 'pythia-mom', 'Label', 'Label_Name'], label_='Label'):
    X = data.drop(drop_labels, 1)
    Y=data[label_]
    return X, Y