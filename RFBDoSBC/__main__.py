

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
from utility import *

import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.family':'DejaVu Sans'})
matplotlib.use('Agg')

#====================================================

doGridSearch=False
feat_imp={}

if __name__=='__main__':

    # *************************Setting Grid Search****************************
    GridArg=0
    if(len(sys.argv)>=2):
        GridArg = int(sys.argv[1])
        msg("Doing Grid Search" if GridArg else "Not Doing Grid Search")
    doGridSearch=True if GridArg else False

    # *************************Setting batch***********************************
    doBatch = 0
    if(sys.argv[2]=='-b'):
        doBatch=1
        msg("Processing data in batches!")
    

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
    
    train = None
    test = None
    batch_num=0

    #*************************Main Loop through Jet Radii and pT hardmin*************************

    for rad in R:
        for ptm in pThardmin:
            makeNewModel=True
            # Clean up from the previous data file

            msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
            log("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))

            best_params = loadBestParameters(rad, ptm)
            rfModel = makeRandomForest(best_params)

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Batch loop~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            while True: # This while loop emulates a do...while ;;;; Matches if(not doBatch)
                #Load data files into train, test
                if(doBatch):
                    batch_num+=1
                    dataGen = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm)
                    train, test = next(dataGen)
                else:
                    train, test = DataPipeline("Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-%1.1f-pThardmin-%d.0.csv"%(rad, ptm), ptm, rows)
                    best_params = doGridSearchOrLoadBestParams(X, Y, rad, ptm)
                # In order to use the sklearn random forest we need a feature vector X
                # and a label vector Y
                X, Y = split_feat_label(train)
                Xtest, Ytest = split_feat_label(test)

                doDataExploration(X, Y, train, rad, ptm, batch_num)
                               
                rfModel = doRandomForestFit(X, Y, rfModel, batch_num)

                oracle = doOracleFit(X, rfModel, rad, ptm)
                
                doModelEvaluation(X, Y, Xtest, Ytest, rfModel, rad, ptm)


                log(" ")
                log(" ")
                log(" ")
                log(" ")
                if(not doBatch):
                    break

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Batch loop~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~            

        with open("Objects/R=%1.1f/feat_imp.pickle"%rad, 'wb') as fil:
            pickle.dump(feat_imp, fil)

        plt.close('all')



