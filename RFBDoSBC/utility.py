import pandas as pd

import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib
import gc
import os
import sys

from matplotlib import cm, colors
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import silhouette_score
from scipy.stats import mode
from .GetAndPrepareData import *
from .modelPreparation import *
from .modelEvaluation import *
from .plotData import *

def msg(_msg):
    print("JEB: %s"%_msg)

def log(_msg, filename="info.txt"):
    with open(filename, "a") as f:
        f.write(_msg+'\n')

def doDataExploration(X, Y, train, rad, ptm, batch_num, do_silhouette_score=False):
    #===============Logging Results============
    log(f"Batch Number: {batch_num}")
    log("Number of Features: %d"%len(X.columns))
    log("Features: "+str(X.columns))
    log("Number of Fake Jets: %d"%len(train.loc[train['Label']==1]))
    log("Number of Squishy Jets: %d"%len(train.loc[train['Label']==2]))
    log("Number of Real Jets: %d"%len(train.loc[train['Label']==3]))
    if do_silhouette_score:
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
    plot_everything(train, rad, ptm, batch_num)
    #Ready for batch

def doGridSearchOrLoadBestParams(X, Y, rad, ptm):
    #This was written to find the best parameters for the randomforest. In batch mode it will find the "best parameters" for each batch.
    if doGridSearch:
        msg("Doing GridSearch")
        best_params = GridSearchHandler(X, Y, rad, ptm)
    else:
        msg("Loading pre-built parameters")
        best_params = loadBestParameters(rad, ptm)
    return best_params
    #Ready for batch

def makeRandomForest(best_params):
    rfModel = RandomForestClassifier(**best_params)
    return rfModel

def doRandomForestFit(X, Y, rfModel, batch_num):
    msg(f"Fitting to data for batch {batch_num}.")
    rfModel.fit(X,Y)
    return rfModel

def doOracleFit(X, rfModel, rad, ptm, low_pT=False):
    msg("Initializing Oracle and fitting to rfModel's predictions")
    oracle = DecisionTreeClassifier(max_depth=3)
    predictions = np.array(list(map(lambda mod: mod.predict(X), rfModel)))
    print(predictions.shape, predictions)
    rfModel_guess = pd.Series(mode(predictions).mode[0])
    print(rfModel_guess.shape, rfModel_guess)
    oracle.fit(X, rfModel_guess)
    display_single_tree(oracle, X, rfModel_guess, "Oracle%s_%1.1f_%d"%( "" if not low_pT else "Lowpt",rad, ptm))
    return oracle

def doModelEvaluation(X, Y, Xtest, Ytest, rfModel, rad, ptm, feat_imp):
    msg("Computing feature importances and other model statistics.")
    importances, indices, std, quant_75 = compute_importance_statistics(rfModel)

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
