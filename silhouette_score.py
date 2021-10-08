import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.decomposition import PCA
from RFBDoSBC.GetAndPrepareData import *
from RFBDoSBC.modelPreparation import *
from RFBDoSBC.modelEvaluation import *
from RFBDoSBC.plotData import *
from RFBDoSBC.utility import *
from sklearn.preprocessing import MinMaxScaler

matplotlib.rcParams.update({"font.size":16})

rad_mom_gen = res_pT_iterator()

while True:

    rad, ptm = next(rad_mom_gen)
    msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
    dataGen = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm, True, False)

    batch_num = 0

    fig, ax = plt.subplots(1, 4, figsize=(15,5))

    while batch_num<1:
        try:
            print(f"Doing batch_num {batch_num}")
            
            train, test = next(dataGen)
            X, Y = split_feat_label(train)

            X = X.loc[Y!=2]
            Y = Y.loc[Y!=2]

            print(X.columns)
            #scaler = MinMaxScaler() 
            Xpt = X.pop('p_T')
            Xptcorr = X.pop('p_T-corr')
            Xpythiapt = X.pop('pythia-mom')
            #Xnorm = pd.DataFrame(scaler.fit_transform(X))


            #ss_vals = silhouette_samples(Xnorm, Y, n_jobs=-1)

            fake_inds = true_indices_for_class(1, Y)
            real_inds = true_indices_for_class(3, Y)

            all_inds = np.concatenate([fake_inds, real_inds])

            #real_ss_vals = ss_vals[real_inds]
            #fake_ss_vals = ss_vals[fake_inds]


            ax[0].hexbin(Xptcorr.iloc[all_inds], X["Area"].iloc[all_inds], gridsize=20, bins='log')
            
            ax[1].hexbin(Xptcorr.iloc[all_inds], X["Angularity"].iloc[all_inds], gridsize=20, bins='log')

            ax[2].hexbin(Xptcorr.iloc[all_inds], X["Mean-p_T"].iloc[all_inds], gridsize=20, bins='log')

            ax[3].hexbin(Xptcorr.iloc[all_inds], X["p_T_1"].iloc[all_inds], gridsize=20, bins='log')

            batch_num+=1
        except StopIteration:
            print(f"No more data in batch num: {batch_num}, trimming forest.")
            break

    ax[0].set_title("Area vs. Jet p_T")
    ax[0].set_xlabel("Corrected Jet p_T")
    ax[0].set_ylabel("Area")
    ax[1].set_title("Angularity vs. Jet p_T")
    ax[1].set_xlabel("Corrected Jet p_T")
    ax[1].set_ylabel("Angularity")
    ax[2].set_title("Mean p_T vs. Jet p_T")
    ax[2].set_xlabel("Corrected Jet p_T")
    ax[2].set_ylabel("Mean p_T")
    ax[3].set_title("Leading Hadron p_T vs. Jet p_T")
    ax[3].set_xlabel("Corrected Jet p_T")
    ax[3].set_ylabel("Leading Hadron p_T")

    
    plt.suptitle(f"R = {rad}, p_T hard min = {ptm}")
    #plt.show()
    plt.tight_layout()
    plt.savefig(f"distributions/inclusive/R{rad}_pt{ptm}.png")

