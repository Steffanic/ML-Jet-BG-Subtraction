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
from pandas_profiling import ProfileReport
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.inspection import permutation_importance
from GetAndPrepareData import *
from modelPreparation import *
from modelEvaluation import *
from plotData import *
from utility import msg

import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.family':'DejaVu Sans'})

if __name__=='__main__':

    #*************************Setting Grid Search****************************
    GridArg=0
    if(len(sys.argv)>=2):
        GridArg = int(sys.argv[1])
        msg("Doing Grid Search" if GridArg else "Not Doing Grid Search")
    doGridSearch=True if GridArg else False

    #*************************Setting number of Rows**************************
    if(len(sys.argv)>=3):
        rows = int(sys.argv[2])
        msg("Only reading %d rows from file."%rows)
    else:
        rows=1000000000000000000

    #Multithreading for the Gridsearch?
    #Thread(traget=GridSearchHandler, args=(X, Y, rad, ptm))


    # Select the jet radii and pthardmin that correspond to the data file you wish to look at
    R = [0.2, 0.3, 0.4, 0.5, 0.6]
    pThardmin = [10,20,30,40]
    Rpt = [(0.2, 10), (0.2, 20), (0.2, 30), (0.2, 40), (0.5, 10), (0.5, 20), (0.5, 30), (0.5, 40), (0.3, 10), (0.3, 20), (0.3, 30), (0.3, 40), (0.4, 10), (0.4, 20), (0.4, 30), (0.4, 40), (0.6, 10), (0.6, 20), (0.6, 30), (0.6, 40)]

    '''
    with Pool(8) as p:
        p.starmap(GridSearchHandler, [(X, Y, r, pt) for r, pt in Rpt])
    '''

    feat_imp = {}
    train = None
    test = None

    #*************************Main Loop through Jet Radii and pT hardmin*************************

    for rad in R:
        for ptm in pThardmin:

            # Clean up from the previous data file
            if(train is not None):
                del train, test
                msg("Deleting train, test from the previous run.")
            gc.collect()

            msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))

            #Load data files into train, test
            train, test = DataPipeline("Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-%1.1f-pThardmin-%d.0.csv"%(rad, ptm), ptm, rows)

            #===============INFO============
            msg("Number of Features: %d"%len(train.columns))
            msg("Number of Fake Jets: %d"%len(train.loc[train['Label']==1]))
            msg("Number of Squishy Jets: %d"%len(train.loc[train['Label']==2]))
            msg("Number of Real Jets: %d"%len(train.loc[train['Label']==3]))
            msg(train.describe())
            #===============================

            # In order to use the sklearn random forest we need a feature vector X
            # and a label vector Y
            X, Y = split_feat_label(train)
            Xtest, Ytest = split_feat_label(test)

            msg("Saving out plots to ./Plots")
            plot_everything(train, rad, ptm)

            if doGridSearch:
                msg("Doing GridSearch")
                best_params = GridSearchHandler(X, Y, rad, ptm)
            else:
                msg("Loading pre-built parameters")
                best_params = loadBestParameters(rad, ptm)


            msg("Initializing model and fitting to data.")
            clf = RandomForestClassifier(**best_params)
            clf.fit(X,Y)
            
            msg("Computing feature importances and other model statistics.")
            importances, indices, std, quant_75 = compute_model_statistics(clf)

            msg("Feature rankings:")
            print_feature_importances(X, importances, indices)

            msg("Plotting feature importances.")
            plot_feature_importances(X, importances, indices, std, rad, ptm)

            msg("Plot distribution of importances.")
            plot_feature_importance_distributions(clf, X, rad, ptm)

            feat_imp["pthard=%d"%ptm] = (X.columns, importances, std, quant_75)

            '''
            i=0
            for col in X.columns:
                msg("Feature Importances:")
                msg(col, clf.estimators_[3].feature_importances_[i])
                i+=1
            '''
            
            #display_single_tree(clf, X, Y, rad, ptm)
            msg("Computing performance metrics")
            compute_performance_metrics(clf, X, Y, Xtest, Ytest, rad, ptm)

        with open("Objects/R=%1.1f/feat_imp.pickle"%rad, 'wb') as fil:
            pickle.dump(feat_imp, fil)

        plt.close('all')