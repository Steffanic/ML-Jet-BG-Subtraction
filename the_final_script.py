import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
import pprint
import pickle
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import silhouette_score, confusion_matrix, accuracy_score
from sklearn.decomposition import PCA
from RFBDoSBC.GetAndPrepareData import *
from RFBDoSBC.modelPreparation import *
from RFBDoSBC.modelEvaluation import *
from RFBDoSBC.plotData import *
from RFBDoSBC.utility import *


rad_mom_gen = res_pT_iterator()
text_results = {}


while True:

    rad, ptm = next(rad_mom_gen)
    msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))

    best_params = loadBestParameters(rad, ptm)

    rfModel = [makeRandomForest(best_params) for _ in range(10)]

    batch_num = 0
    feat_imp = {}

    dataGen = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm, False, False)

    text_results[f'R{rad}pt{ptm}'] = {}

    while batch_num<12:
        try:
            print(f"Doing batch_num {batch_num}")
            
            train, test = next(dataGen)
            
            X, Y = split_feat_label(train)
     
            # Test PCA

            pca = PCA()

            pca.fit(X)

            text_results[f'R{rad}pt{ptm}']['pca'] = [(pca.explained_variance_ratio_[i], dict(zip(X.columns,pca.components_[i]))) for i in range(len(X.columns))]
            
            rfModel[batch_num] = doRandomForestFit(X, Y, rfModel[batch_num], batch_num)

            batch_num+=1

        except StopIteration:
            print(f"No more data in batch num: {batch_num}, trimming forest.")
            rfModel = rfModel[:batch_num]
            break

    
    train, test = DataPipeline(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm, 200000000, False)

    X, Y = split_feat_label(train)

    oracle = doOracleFit(X, rfModel, rad, ptm)


    Xtest, Ytest = split_feat_label(test)
    oracle_test_predictions = oracle.predict(Xtest)
    oracle_acc = accuracy_score(Ytest, oracle_test_predictions)
    print(f"Oracle accuracy: {oracle_acc}")
    text_results[f'R{rad}pt{ptm}']['oracle_acc'] = oracle_acc

    ave_rf_acc=0

    for i,rf in enumerate(rfModel):
        rf_predictions = rf.predict(Xtest)
        rf_acc = accuracy_score(Ytest, rf_predictions)
        ave_rf_acc+=rf_acc
        print(f"Random Forest accuracy, model {i}: {rf_acc}")
    ave_rf_acc/=(i+1)

    text_results[f'R{rad}pt{ptm}']['ave_rf_acc'] = ave_rf_acc

    feat_imps = []
    model_no=0

    ave_perf_metrics = {}

    for model in rfModel:
        feat_imp = {}
        perf_metrics = doModelEvaluation(X, Y, Xtest, Ytest, model, rad, ptm, feat_imp,False, model_no)
        
        if model_no==0:
            ave_perf_metrics=perf_metrics
        else:
            for key in ave_perf_metrics.keys():
                ave_perf_metrics[key] = [old+new for old, new in zip(ave_perf_metrics[key], perf_metrics[key])]
        feat_imps.append(feat_imp[f'pthard={ptm}'][1])
        model_no+=1
    ave_perf_metrics = {key:list(map(lambda x:x/(model_no+1), value)) for (key, value) in zip(ave_perf_metrics.keys(),ave_perf_metrics.values())}
    print(f"Average Performance Metrics:{ave_perf_metrics}")
    text_results[f"R{rad}pt{ptm}"]['ave_perf_metrics'] = ave_perf_metrics
    text_results[f"R{rad}pt{ptm}"]['ave_feat_imps'] = np.mean(feat_imps, axis=0)

    names = X.columns

    top_feat = names[oracle.tree_.feature[0]]
    top_thresh = oracle.tree_.threshold[0]
    print(f"{top_feat} <= {top_thresh}")
    key_names = train['Label_Name'].unique()
    jet_vals = dict(zip(key_names, oracle.tree_.value[0][0]))
    lt_jet_vals = dict(zip(key_names, oracle.tree_.value[oracle.tree_.children_left[0]][0]))
    rt_jet_vals = dict(zip(key_names, oracle.tree_.value[oracle.tree_.children_right[0]][0]))



    print(f"Vals: {jet_vals}, lt_vals: {lt_jet_vals}, rt_vals: {rt_jet_vals}")
    print()

    rejection_rates = [100*lt_jet_vals[lbl]/jet_vals[lbl] for lbl in key_names]
    retention_rates = [100*rt_jet_vals[lbl]/jet_vals[lbl] for lbl in key_names]

    print(f"{rejection_rates=}, {retention_rates=}")

    rejection_string = ", ".join([f"{rt:2.1f}% of our {name} jets" for rt, name in zip(rejection_rates, key_names)])
    retention_string = ", ".join([f"{rt:2.1f}% of our {name} jets" for rt, name in zip(retention_rates, key_names)])

    lt_string = f"{top_feat} <= {top_thresh:2.3f} keeps {rejection_string}"

    rt_string = f"{top_feat} > {top_thresh:2.3f} keeps {retention_string}"


    text_results[f"R{rad}pt{ptm}"]["top_cut"] = (top_feat, top_thresh)
    text_results[f"R{rad}pt{ptm}"]["top_cut_strings"] = (lt_string, rt_string)

    pp = pprint.PrettyPrinter(depth=6)
    pp.pprint(text_results)
    if rad==0.6 and ptm==40:
        print("We made it! Let's save the data.")
        pickle.dump(text_results, open("NumericalData.pickle", "wb"))