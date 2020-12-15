import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib
import gc
import os
import sys
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

import warnings
warnings.filterwarnings("ignore")


matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.family':'DejaVu Sans'})

def msg(_msg):
    print("JEB: %s"%_msg)

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

    

    # Select the jet radii and pthardmin that correspond to the data file you wish to look at
    R = [0.2, 0.3, 0.4, 0.5, 0.6]
    pThardmin = [10,20,30,40]

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
            msg("Number of Fake Jets: %d"%len(train.loc[train['Label']==0]))
            msg("Number of Squishy Jets: %d"%len(train.loc[train['Label']==1]))
            msg("Number of Real Jets: %d"%len(train.loc[train['Label']==2]))
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

            plt.figure()
            plt.title("Feature importances")
            plt.bar(range(X.shape[1]), importances[indices],
                    color="r", yerr=std[indices], align="center")
            plt.xticks(range(X.shape[1]), X.columns[indices])
            plt.xlim([-1, X.shape[1]])
            plt.show()

            #Plot distribution of importances
            plt.figure(figsize=(10,10))
            feat_imp_dist = [tree.feature_importances_ for tree in clf.estimators_]
            
            phi = np.array(list(zip(*feat_imp_dist)))[1]*3.1415

            rgb_cycle = np.vstack((            # Three sinusoids
                .5*(1.+np.cos(phi          )), # scaled to [0,1]
                .5*(1.+np.cos(phi+2*np.pi/3)), # 120° phase shifted.
                .5*(1.+np.cos(phi-2*np.pi/3)))).T 
            i=0
            for feats in zip(*feat_imp_dist):
                plt.scatter(np.ones(len(feats))*i+(np.random.randn(len(feats))*0.05), feats, c=rgb_cycle, s=50, alpha=0.2)
                i+=1
            plt.xticks(np.linspace(0, len(X.columns)-1, len(X.columns)), X.columns)
            plt.title("Distribution of Feature Importances in Random Forest", fontsize=16)
            plt.xlabel("Features", fontsize=14)
            plt.ylabel("Feature Importance", fontsize=14)
            plt.xlim(-0.5, len(X.columns)+0.5)
            plt.ylim(-0.01, 1.01)
            plt.savefig("Plots/R=%1.1f/FeatureImportancesDistribution_pthard=%d.png"%(rad, ptm))
            #plt.show()

            feat_imp["pthard=%d"%ptm] = (X.columns, importances, std, quant_75)

            '''
            i=0
            for col in X.columns:
                msg("Feature Importances:")
                msg(col, clf.estimators_[3].feature_importances_[i])
                i+=1
            '''
            
            #msg(gcv.best_estimator_.get_params())
            # Extract single tree
            estimator = clf.estimators_[3]
            #estimator.fit(X, Y)

            from sklearn.tree import export_graphviz
            # Export as dot file
            export_graphviz(estimator, out_file='tree.dot', 
                            feature_names = X.columns,
                            class_names = list(map(lambda x:str(x), train['Label'].unique())),
                            rounded = True, proportion = False, 
                            precision = 2, filled = True)

            # Convert to png using system command (requires Graphviz)
            from subprocess import call
            call(['dot', '-Tpng', 'tree.dot', '-o', 'Plots/R=%1.1f/tree_pTmin%d.png'%(rad, ptm), '-Gdpi=600'])


            #Plot performance metrics for each radius and pthard

            perf_metrics = {}

            y_pred_train = np.array(clf.predict(X))
            msg(y_pred_train)
            msg(len(np.array(Y==2).nonzero()[0]))
            real_jet_ind_train = np.array(Y==2).nonzero()[0]
            fake_jet_ind_train = np.array(Y==0).nonzero()[0]
            squish_jet_ind_train = np.asarray(Y==1).nonzero()[0]
            pred_for_real_jets_train = y_pred_train[real_jet_ind_train]
            pred_for_fake_jets_train = y_pred_train[fake_jet_ind_train]
            pred_for_squish_jets_train = y_pred_train[squish_jet_ind_train] 

            msg(len(real_jet_ind_train),  np.asarray(pred_for_real_jets_train==2).nonzero()[0])

            real_jet_rate_train = float(len(np.asarray(pred_for_real_jets_train==2).nonzero()[0]))/float(len(real_jet_ind_train))
            fake_jet_rate_train = float(len(np.asarray(pred_for_fake_jets_train==0).nonzero()[0]))/float(len(fake_jet_ind_train))
            squish_jet_rate_train = float(len(np.asarray(pred_for_squish_jets_train==1).nonzero()[0]))/float(len(squish_jet_ind_train))

            fall_out_train = float(len(np.asarray(pred_for_fake_jets_train==2).nonzero()[0]))/float(len(fake_jet_ind_train))
            fall_in_train = float(len(np.asarray(pred_for_real_jets_train==0).nonzero()[0]))/float(len(real_jet_ind_train))

            perf_metrics['real_jet_rate_train'] = real_jet_rate_train
            perf_metrics['fake_jet_rate_train'] = fake_jet_rate_train
            perf_metrics['squish_jet_rate_train'] = squish_jet_rate_train
            perf_metrics['fall_out_train'] = fall_out_train
            perf_metrics['fall_in_train'] = fall_in_train

            msg("Real Jet Rate train: ", real_jet_rate_train, "\nFake Jet Rate train: ", fake_jet_rate_train, "\nSquish Jet Rate train: ", squish_jet_rate_train, "\nFake predicted real train: ", fall_out_train, "\nReal predicted fake train: ", fall_in_train)

            X_test, Y_test = split_feat_label(test)

            #Real Jet rate = real jets pred/real jets
            y_pred_test = np.array(clf.predict(X_test))
            msg(y_pred_test)
            real_jet_ind_test = np.asarray(Y_test==2).nonzero()[0]
            fake_jet_ind_test = np.asarray(Y_test==0).nonzero()[0]
            squish_jet_ind_test = np.asarray(Y_test==1).nonzero()[0]
            pred_for_real_jets_test = y_pred_test[real_jet_ind_test]
            pred_for_fake_jets_test = y_pred_test[fake_jet_ind_test]
            pred_for_squish_jets_test = y_pred_test[squish_jet_ind_test] 

            msg(len(real_jet_ind_test), len(np.asarray(pred_for_real_jets_test==2).nonzero()[0]))

            real_jet_rate_test = float(len(np.asarray(pred_for_real_jets_test==2).nonzero()[0]))/float(len(real_jet_ind_test))
            fake_jet_rate_test = float(len(np.asarray(pred_for_fake_jets_test==0).nonzero()[0]))/float(len(fake_jet_ind_test))
            squish_jet_rate_test = float(len(np.asarray(pred_for_squish_jets_test==1).nonzero()[0]))/float(len(squish_jet_ind_test))

            fall_out_test = float(len(np.asarray(pred_for_fake_jets_test==2).nonzero()[0]))/float(len(fake_jet_ind_test))
            fall_in_test = float(len(np.asarray(pred_for_real_jets_test==0).nonzero()[0]))/float(len(real_jet_ind_test))


            perf_metrics['real_jet_rate_test'] = real_jet_rate_test
            perf_metrics['fake_jet_rate_test'] = fake_jet_rate_test
            perf_metrics['squish_jet_rate_test'] = squish_jet_rate_test
            perf_metrics['fall_out_test'] = fall_out_test
            perf_metrics['fall_in_test'] = fall_in_test


            msg("Real Jet Rate test: ", real_jet_rate_test, "\nFake Jet Rate test: ", fake_jet_rate_test, "\nSquish Jet Rate test: ", squish_jet_rate_test, "\nFake predicted real test: ", fall_out_test, "\nReal predicted fake test: ", fall_in_test)


            if(not os.path.isdir("Objects/R=%1.1f"%rad)):
                os.mkdir("Objects/R=%1.1f"%rad)
            with open("Objects/R=%1.1f/perf_metric_pthard=%d"%(rad, ptm), 'wb') as perffile:
                pickle.dump(perf_metrics, perffile)

            #plt.show()
        
        with open("Objects/R=%1.1f/feat_imp.pickle"%rad, 'wb') as fil:
            pickle.dump(feat_imp, fil)

        plt.close('all')