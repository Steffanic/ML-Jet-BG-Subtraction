# ~~~~~~~~~~~ Patrick Steffanic ~~~~~~~~~~~~~~~~~~
# =================IMPORTS and FRONT MATTER=======
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib
import gc
import os
import sys
from threading import Thread
from multiprocessing import Pool
from matplotlib import cm, colors
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import silhouette_score
from GetAndPrepareData import *
from modelPreparation import *
from modelEvaluation import *
from plotData import *
from utility import msg, logs

import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.family':'DejaVu Sans'})
matplotlib.use('Agg')

#====================================================

doGridSearch=False

if __name__=='__main__':

    # *************************Setting Grid Search****************************
    GridArg=0
    if(len(sys.argv)>=2):
        GridArg = int(sys.argv[1])
        msg("Doing Grid Search" if GridArg else "Not Doing Grid Search")
    global doGridSearch=True if GridArg else False

    # *************************Setting batch***********************************
    doBatch = 0
    if(sys.argv[2]=='-b'){
        doBatch=1
        msg("Processing data in batches!")
    }

    # *************************Setting number of Rows**************************
    if(len(sys.argv)>=4):
        rows = int(sys.argv[3])
        msg("Only reading %d rows from file."%rows)
    else:
        rows=1000000000000000000

    log("Number of rows read: %d"%rows)


    # Select the jet radii and pthardmin that correspond to the data file you wish to look at
    R = [0.2, 0.3, 0.4, 0.5, 0.6]
    pThardmin = [10,20,30,40]
    Rpt = [(0.2, 10), (0.2, 20), (0.2, 30), (0.2, 40), (0.5, 10), (0.5, 20), (0.5, 30), (0.5, 40), (0.3, 10), (0.3, 20), (0.3, 30), (0.3, 40), (0.4, 10), (0.4, 20), (0.4, 30), (0.4, 40), (0.6, 10), (0.6, 20), (0.6, 30), (0.6, 40)]
    
    feat_imp = {}
    train = None
    test = None
    batch_num=0

    #*************************Main Loop through Jet Radii and pT hardmin*************************

    for rad in R:
        for ptm in pThardmin:

            # Clean up from the previous data file
            if(train is not None):
                del train, test
                msg("Deleting train, test from the previous run.")
            gc.collect()

            msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
            log("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))

            #Load data files into train, test
            train, test = DataPipeline("Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-%1.1f-pThardmin-%d.0.csv"%(rad, ptm), ptm, rows)


            # In order to use the sklearn random forest we need a feature vector X
            # and a label vector Y
            X, Y = split_feat_label(train)
            Xtest, Ytest = split_feat_label(test)


            doDataExploration(X, Y, train, rad, ptm)
            
            best_params = doGridSearchOrLoadBestParams(X, Y, rad, ptm)

            rfModel = doRandomForestFit(X, Y, best_params)

            oracle = doOracleFit(X, rfModel, rad, ptm)
            
            doModelEvaluation(X, Y, Xtest, Ytest, rfModel, rad, ptm)


            log(" ")
            log(" ")
            log(" ")
            log(" ")
        with open("Objects/R=%1.1f/feat_imp.pickle"%rad, 'wb') as fil:
            pickle.dump(feat_imp, fil)

        plt.close('all')



def doDataExploration(X, Y, train, rad, ptm, batch_num):
    #===============Logging Results============
    log(f"Batch Number: {batch_num}")
    log("Number of Features: %d"%len(X.columns))
    log("Features: "+str(X.columns))
    log("Number of Fake Jets: %d"%len(train.loc[train['Label']==1]))
    log("Number of Squishy Jets: %d"%len(train.loc[train['Label']==2]))
    log("Number of Real Jets: %d"%len(train.loc[train['Label']==3]))
    log("Silhouette Scores:")
    log(f"Overall Score: {silhouette_score(X, Y)}")
    log(f"Per Column Scores: ")
    for col in X.columns:
        log(f"Column: {col}  Score: {silhouette_score(np.array(X[col]).reshape((-1,1)), Y)}")
    #==========================================

    #===============INFO============
    msg("Number of Features: %d"%len(X.columns))
    msg("Number of Fake Jets: %d"%len(train.loc[train['Label']==1]))
    msg("Number of Squishy Jets: %d"%len(train.loc[train['Label']==2]))
    msg("Number of Real Jets: %d"%len(train.loc[train['Label']==3]))
    msg(train.describe())
    #===============================

    msg("Saving out plots to ./Plots")
    plot_everything(train, rad, ptm)

def doGridSearchOrLoadBestParams(X, Y, rad, ptm):
    if global doGridSearch:
        msg("Doing GridSearch")
        best_params = GridSearchHandler(X, Y, rad, ptm)
    else:
        msg("Loading pre-built parameters")
        best_params = loadBestParameters(rad, ptm)
    return best_params

def doRandomForestFit(X, Y, best_params):
    msg("Initializing model and fitting to data.")
    rfModel = RandomForestClassifier(**best_params)
    rfModel.fit(X,Y)
    return rfModel

def doOracleFit(X, rfModel, rad, ptm):
    msg("Initializing Oracle and fitting to rfModel's predictions")
    oracle = DecisionTreeClassifier(max_depth=3)
    rfModel_guess = pd.Series(rfModel.predict(X))
    oracle.fit(X, rfModel_guess)
    display_single_tree(oracle, X, rfModel_guess, "Oracle_%1.1f_%d"%(rad, ptm))
    return oracle

def doModelEvaluation(X, Y, Xtest, Ytest, rfModel, rad, ptm):
    msg("Computing feature importances and other model statistics.")
    importances, indices, std, quant_75 = compute_model_statistics(rfModel)

    log("Feature importances: "+str(list(zip(X.columns, importances))))

    msg("Feature rankings:")
    print_feature_importances(X, importances, indices)

    msg("Plotting feature importances.")
    plot_feature_importances(X, importances, indices, std, rad, ptm)

    msg("Plot distribution of importances.")
    plot_feature_importance_distributions(rfModel, X, rad, ptm)

    feat_imp["pthard=%d"%ptm] = (X.columns, importances, std, quant_75)
    
    #display_single_tree(rfModel, X, Y, rad, ptm)
    msg("Computing performance metrics")
    compute_performance_metrics(rfModel, X, Y, Xtest, Ytest, rad, ptm)