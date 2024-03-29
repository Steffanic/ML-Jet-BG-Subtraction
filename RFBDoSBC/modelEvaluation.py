# -*- coding: utf-8 -*-
import numpy as np
from sklearn import tree
import matplotlib.pyplot as plt

import os
import pickle

def log(_msg, filename="info.txt"):
    with open(filename, "a") as f:
        f.write(_msg+'\n')

def compute_importance_statistics(clf):
    importances = clf.feature_importances_
    #permutation_importance(clf, Xtest, Ytest, n_repeats=10, random_state=42)
    std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)
    quant_75 = np.quantile([tree.feature_importances_ for tree in clf.estimators_],0.75, axis=0)
    indices = np.argsort(importances)[::-1]
    return importances, indices, std, quant_75

def print_feature_importances(X, importances, indices):
    for f in range(X.shape[1]):
        print("%d. %s (%f)" % (f + 1, X.columns[indices[f]], importances[indices[f]]))

def plot_feature_importances(X, importances, indices, std, rad, ptm, low_pt=False, model_no=0):
    plt.figure()
    plt.title("Feature importances")
    plt.bar(range(X.shape[1]), importances[indices],
            color="r", yerr=std[indices], align="center")
    plt.xticks(range(X.shape[1]), X.columns[indices])
    plt.xlim([-1, X.shape[1]])
    plt.savefig("Plots/R=%1.1f/FeatureImportances_modelno%d_pthard=%d_%s.png"%(rad, model_no, ptm, "" if not low_pt else "low_pt"))

def plot_feature_importance_distributions(clf, X, rad, ptm, low_pt=False, model_no=0):
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
    plt.savefig("Plots/R=%1.1f/FeatureImportancesDistribution_modelno%d_pthard=%d_%s.png"%(rad, model_no, ptm, "" if not low_pt else "low_pt"))


def display_single_tree(_tree, X, Y, name):
    #msg(gcv.best_estimator_.get_params())
    # Extract single tree
    estimator = _tree
    #estimator.fit(X, Y)
    classes = np.array(["", "Fake", "Squish", "Real"])
    print(np.sort(Y.unique()))

    from sklearn.tree import export_graphviz
    # Export as dot file
    export_graphviz(estimator, out_file='Oracles/'+name+'.dot', 
                    feature_names = X.columns,
                    class_names=classes[np.sort(Y.unique())],
                    rounded = True, proportion = False, 
                    precision = 2, filled = True)

    # Convert to png using system command (requires Graphviz)
    from subprocess import call
    call(['dot', '-Tpng', 'Oracles/'+name+'.dot', '-o', 'Oracles/'+name+'.png', '-Gdpi=600'])


def compute_performance_metrics(clf, X, Y, X_test, Y_test, rad, ptm):
    perf_metrics = {}

    # Get the predictions on the training set
    y_pred_train = np.array(clf.predict(X))

    predictions_train = []

    for lbl in Y.unique():
        # This returns the predicted labels for the pure populations of Y==lbl
        predictions_train.append(predictions_for_class(lbl, y_pred_train, Y))
    print(f"{predictions_train=}")

    population_rates_train = []

    calc_rate = lambda cls_num, df, indices: float(len(np.asarray(df==cls_num).nonzero()[0]))/float(len(indices))

    for lbl, pop in zip(Y.unique(),predictions_train):
        population_rates_train.append(calc_rate(lbl, pop, true_indices_for_class(lbl, Y))) # Sensitivity and Specifity in a two-class system

    print(f"{population_rates_train=}")

    perf_metrics['population_rates_train'] = population_rates_train


    y_pred_test = np.array(clf.predict(X_test))

    predictions_test = []

    for lbl in Y_test.unique():
        # This returns the predicted labels for the pure populations of Y==lbl
        predictions_test.append(predictions_for_class(lbl, y_pred_test, Y_test))
    print(f"{predictions_test=}")

    population_rates_test = []

    calc_rate = lambda cls_num, df, indices: float(len(np.asarray(df==cls_num).nonzero()[0]))/float(len(indices))

    for lbl, pop in zip(Y_test.unique(), predictions_test):
        population_rates_test.append(calc_rate(lbl, pop, true_indices_for_class(lbl, Y_test))) # Sensitivity and Specifity in a two-class system

    print(f"{population_rates_test=}")

    perf_metrics['population_rates_test'] = population_rates_test


    '''
    # Partition the data by classes
    print(f"Number of Real Jets train: {len(np.array(Y==3).nonzero()[0])}")
    print(f"Number of Fake Jets train: {len(np.array(Y==1).nonzero()[0])}")
    real_jet_ind_train = np.array(Y==3).nonzero()[0]
    fake_jet_ind_train = np.array(Y==1).nonzero()[0]
    squish_jet_ind_train = np.asarray(Y==2).nonzero()[0]
    pred_for_real_jets_train = y_pred_train[real_jet_ind_train]
    pred_for_fake_jets_train = y_pred_train[fake_jet_ind_train]
    pred_for_squish_jets_train = y_pred_train[squish_jet_ind_train] 
    
    #print(str(len(real_jet_ind_train))+' '+str(np.asarray(pred_for_real_jets_train==3).nonzero()[0]))

    real_jet_rate_train = float(len(np.asarray(pred_for_real_jets_train==3).nonzero()[0]))/float(len(real_jet_ind_train))
    fake_jet_rate_train = float(len(np.asarray(pred_for_fake_jets_train==1).nonzero()[0]))/float(len(fake_jet_ind_train))
    squish_jet_rate_train = float(len(np.asarray(pred_for_squish_jets_train==2).nonzero()[0]))/float(len(squish_jet_ind_train))
    
    fall_out_train = float(len(np.asarray(pred_for_fake_jets_train==3).nonzero()[0]))/float(len(fake_jet_ind_train))
    fall_in_train = float(len(np.asarray(pred_for_real_jets_train==1).nonzero()[0]))/float(len(real_jet_ind_train))
    '''

    '''
    perf_metrics['real_jet_rate_train'] = real_jet_rate_train
    perf_metrics['fake_jet_rate_train'] = fake_jet_rate_train
    perf_metrics['squish_jet_rate_train'] = squish_jet_rate_train
    perf_metrics['fall_out_train'] = fall_out_train
    perf_metrics['fall_in_train'] = fall_in_train

    print("Real Jet Rate train: "+str(real_jet_rate_train)+"\nFake Jet Rate train: "+str(fake_jet_rate_train)+"\nSquish Jet Rate train: "+str(squish_jet_rate_train)+"\nFake predicted real train: "+str(fall_out_train)+"\nReal predicted fake train: "+str(fall_in_train))
    log("Real Jet Rate train: "+str(real_jet_rate_train)+"\nFake Jet Rate train: "+str(fake_jet_rate_train)+"\nSquish Jet Rate train: "+str(squish_jet_rate_train)+"\nFake predicted real train: "+str(fall_out_train)+"\nReal predicted fake train: "+str(fall_in_train))
    '''


    #Real Jet rate = real jets pred/real jets

    '''
    print(f"Number of Real Jets test: {len(np.array(Y_test==3).nonzero()[0])}")
    print(f"Number of Fake Jets test: {len(np.array(Y_test==1).nonzero()[0])}")
    real_jet_ind_test = np.asarray(Y_test==3).nonzero()[0]
    fake_jet_ind_test = np.asarray(Y_test==1).nonzero()[0]
    squish_jet_ind_test = np.asarray(Y_test==2).nonzero()[0]
    pred_for_real_jets_test = y_pred_test[real_jet_ind_test]
    pred_for_fake_jets_test = y_pred_test[fake_jet_ind_test]
    pred_for_squish_jets_test = y_pred_test[squish_jet_ind_test] 

    #print(str(len(real_jet_ind_test))+' '+str(len(np.asarray(pred_for_real_jets_test==3).nonzero()[0])))

    real_jet_rate_test = float(len(np.asarray(pred_for_real_jets_test==3).nonzero()[0]))/float(len(real_jet_ind_test))
    fake_jet_rate_test = float(len(np.asarray(pred_for_fake_jets_test==1).nonzero()[0]))/float(len(fake_jet_ind_test))
    squish_jet_rate_test = float(len(np.asarray(pred_for_squish_jets_test==2).nonzero()[0]))/float(len(squish_jet_ind_test))

    fall_out_test = float(len(np.asarray(pred_for_fake_jets_test==3).nonzero()[0]))/float(len(fake_jet_ind_test))
    fall_in_test = float(len(np.asarray(pred_for_real_jets_test==1).nonzero()[0]))/float(len(real_jet_ind_test))


    perf_metrics['real_jet_rate_test'] = real_jet_rate_test
    perf_metrics['fake_jet_rate_test'] = fake_jet_rate_test
    perf_metrics['squish_jet_rate_test'] = squish_jet_rate_test
    perf_metrics['fall_out_test'] = fall_out_test
    perf_metrics['fall_in_test'] = fall_in_test


    print("Real Jet Rate test: "+str(real_jet_rate_test)+"\nFake Jet Rate test: "+str(fake_jet_rate_test)+"\nSquish Jet Rate test: "+str(squish_jet_rate_test)+"\nFake predicted real test: "+str(fall_out_test)+"\nReal predicted fake test: "+str(fall_in_test))
    log("Real Jet Rate test: "+str(real_jet_rate_test)+"\nFake Jet Rate test: "+str(fake_jet_rate_test)+"\nSquish Jet Rate test: "+str(squish_jet_rate_test)+"\nFake predicted real test: "+str(fall_out_test)+"\nReal predicted fake test: "+str(fall_in_test))
    '''

    if(not os.path.isdir("Objects/R=%1.1f"%rad)):
        os.mkdir("Objects/R=%1.1f"%rad)
    with open("Objects/R=%1.1f/perf_metric_pthard=%d"%(rad, ptm), 'wb') as perffile:
        pickle.dump(perf_metrics, perffile)

    return perf_metrics

def predictions_for_class(class_num, df, label):
    '''This extracts the predicted values for class class_num. 
    
    Arguments:
    class_num -- integer label for the class you want.
    df -- dataframe containing predictions for all classes.
    label -- Series which only contains class labels

    Returns:
    class_exclusive_df -- dataframe containing only the rows of df which are in class, class_num
    '''

    class_indices = true_indices_for_class(class_num, label)
    class_exclusive_df = df[class_indices]
    return class_exclusive_df

def true_indices_for_class(class_num, label):
    '''This extracts the indices for class class_num in label. It can be used as a mask.
    
    Arguments:
    class_num -- integer label for the class you want.
    label -- Series which only contains class labels

    Returns:
    class_indices -- dataframe containing only the indices which are in class, class_num
    '''
    return np.array(label==class_num).nonzero()[0]