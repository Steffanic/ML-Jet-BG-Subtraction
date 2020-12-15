import numpy as np
from sklearn import tree


def compute_model_statistics(clf):
    importances = clf.feature_importances_
    #permutation_importance(clf, Xtest, Ytest, n_repeats=10, random_state=42)
    std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)
    quant_75 = np.quantile([tree.feature_importances_ for tree in clf.estimators_],0.75, axis=0)
    indices = np.argsort(importances)[::-1]
    return importances, indices, std, quant_75

def print_feature_importances(X, importances, indices):
    for f in range(X.shape[1]):
        print("%d. %s (%f)" % (f + 1, X.columns[indices[f]], importances[indices[f]]))

