import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
import sys

from RFBDoSBC.GetAndPrepareData import *
from RFBDoSBC.modelPreparation import *
from RFBDoSBC.modelEvaluation import *
from RFBDoSBC.plotData import *
from RFBDoSBC.utility import *

DO_FIT = False
DO_PLOTTING = False
LOW_PT = False

# python3 -m train_and_evaluate.py [ --fit | --plot | --lowpt ] 

def parse_command_line_args():
    if '--fit' in sys.argv:
        msg("Fitting is on.")
        DO_FIT = True
    if '--plot' in sys.argv:
        msg("Plotting is on.")
        if DO_FIT:
            msg("Plotting and fitting are on. Watch your RAM go 🚀 🌕")
        DO_PLOTTING = True
    if '--lowpt' in sys.argv:
        msg("Low-p_T sub-dataset activated.")
        LOW_PT = True

if __name__=='__main__':

    parse_command_line_args()

    # First, grab the iterator for our (resolution parameter, p_T hard min) pairs
    #
    res_mom_it = res_pT_iterator()

    for res, ptm in res_mom_it:

        msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(res, ptm))

        # I have saved a set of parameters that produces decent quality predictions. 
        # This is worth revisiting
        #
        best_params = loadBestParameters(res, ptm)
        if LOW_PT:
            best_params_lowpt = loadBestParameters(res, ptm, low_pt=True)


        # Initialize 10 RF models. The data is split into N<10 batches 
        # and each forest is trained on a different batch. The untrained models are removed
        #
        rfModel = [makeRandomForest(best_params) for _ in range(15)]
        if LOW_PT:
            rfModel_lowpT = [makeRandomForest(best_params_lowpt) for _ in range(15)]

        
        current_batch = 0
        feature_importances = {}

        

        # Returns an iterator over batches of pre-processed data.
        # If we want a low_pT dataset, we need to keep the pT column to split by.
        #
        data_batch_pT_it = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{res}-pThardmin-{ptm}.0.csv"
                                                , ptm
                                                , keep_pT=LOW_PT
                                                , keep_all=False)


        # We want to plot all the variables, so we don't drop them in DataPipelineBatch
        #
        if DO_PLOTTING:
            data_batch_all_it = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{res}-pThardmin-{ptm}.0.csv"
                                                    , ptm
                                                    , keep_pT=LOW_PT
                                                    , keep_all=True)
            
            fig_all, ax_all = plt.subplots(8+15, 1, figsize=(15,40))
            fig_lowpT, ax_lowpT = plt.subplots(8+15, 1, figsize=(15,40))

        while current_batch<15:
            try:
                print(f"Processing current_batch #{current_batch}")

                i=0

                if DO_FIT:
                    train, test = next(data_batch_pT_it)
                    X, Y = split_feat_label(train)
                    rfModel[current_batch] = doRandomForestFit(X, Y, rfModel[current_batch], current_batch)

                    if LOW_PT:
                        train_lowpT = train.loc[train['p_T']<100]
                        X_lowpT, Y_lowpT = split_feat_label(train_lowpT)
                        rfModel_lowpT[current_batch] = doRandomForestFit(X_lowpT, Y_lowpT, rfModel_lowpT[current_batch], current_batch)
                
                if DO_PLOTTING:
                    train_all_features, test_all_features = next(data_batch_all_it)
                    X_all, Y_all = split_feat_label(train)

                    if LOW_PT:
                        train_all_features_lowpT = train_all_features.loc[train_all_features['p_T']<100]
                        X_all_lowpT, Y_all_lowpT = split_feat_label(train_lowpT)

                if DO_PLOTTING:
                    print("Plotting...")
                    for feat in train_all_features.columns:

                        if feat in ['p_T', 'Label_Name']:
                            # Don't plot p_T vs p_T or vs Label_Name
                            continue

                        sns.scatterplot(data=train_all_features, x='p_T',y=feat, hue='Label_Name',hue_order=["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"], alpha=0.005, ax=ax_all[i], s=1)
                        ax_all[i].legend([],[])

                        if LOW_PT:
                            sns.scatterplot(data=train_all_features_lowpT , x='p_T',y=feat, hue='Label_Name',hue_order=["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"], alpha=0.005, ax=ax_lowpT[i], s=1)
                            ax_lowpT[i].legend([],[])  
                    
                    i+=1

                current_batch+=1
            except StopIteration:
                print(f"No more data in batch num: {current_batch}, trimming forest.")
                rfModel = rfModel[:current_batch]
                if LOW_PT:
                    rfModel_lowpT = rfModel_lowpT[:current_batch]
                break

        if DO_PLOTTING:
            print("Saving...")

            fig_all.legend(["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"])
            fig_all.savefig(f'all_jet_features_R{res}_pt{ptm}.png', format='png', dpi=600)
           
            if LOW_PT:
                fig_lowpT.legend(["Fake", "< Hard Scattering p_T", ">= Hard Scattering p_T"])
                fig_lowpT.savefig(f"lowp_T_jet_features_R{res}_pt{ptm}.png", format='png', dpi=600)
    if DO_FIT:
        # Now get the full dataset and train the oracle on the mode of the predictions from our random forests
        #
        train, test = DataPipeline(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{res}-pThardmin-{ptm}.0.csv", ptm, 200000000, True)
        X, Y = split_feat_label(train)

        if LOW_PT:
            train_lowpT = train.loc[train['p_T']<100]
            X_lowpT, Y_lowpT = split_feat_label(train_lowpT)
