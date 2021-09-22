from sys import float_repr_style
from sklearn.ensemble import  RandomForestClassifier
from sklearn.model_selection import GridSearchCV
import numpy as np
import pickle

def loadBestParameters(rad, ptm, low_pt=False, print_params=False):
    with open("Objects/best_params_%1.1f_%d%s.pickle"%(rad, ptm, "" if not low_pt else "_lowpt"), 'rb') as f:
        best_params = pickle.load(f)
        best_params['max_features']='auto'
        if low_pt:
            best_params['class_weight'] = {1:1, 2:1, 3:4}
        else:
            best_params.pop('class_weight',None)
        if print_params:
            print(best_params)
    return best_params

def GridSearchHandler(X, Y, rad, ptm):
    best_params, cv_results = DoGridSearch(X, Y)
    #print("Best Params:\n"+str(best_params)+ "\nCV Results:\n"+ str(cv_results))

    best_params["oob_score"]=True

    best_params['random_state']=42
    best_params['n_jobs']=-1
    with open("Objects/best_params_%1.1f_%d.pickle"%(rad, ptm), 'wb') as fil:
        pickle.dump(best_params, fil)
    return best_params

def DoGridSearch(X, Y):
    import itertools
    clf = RandomForestClassifier(oob_score=True, random_state=42, n_jobs=-1)
    best_params = {}
    cv_results = {}

    n_est_s = (100,200) # Number of Estimators in the forest

    crit_s = ('gini', 'entropy') # Criterion for evaluating purity

    max_depth = (1, 2,3) # Number of decision layers available to the decision trees

    min_samples_split = (2, 4, 20, 40, 100) # Minimum number of samples required to split an internal node

    min_samples_leaf = (100, 1000, 2000) # Minimum number of samples allowed left after a split

    min_weight_fraction_leaf = (0.1, 0.2, 0.3, 0.45) # ??

    max_features = (1,2,3,4, 5)

    max_leaf_nodes = (10, 100, 1000)

    min_impurity_decrease = (0.0, 0.3, 0.8)

    class_weight = ({1:1, 2:1, 3:10}, {1:1,2:1,3:100}, {1:1, 2:1, 3:1000})

    max_samples = (0.1, 0.2, 0.4, 0.8, 0.9)

    names=np.array(['n_estimators', 'criterion', 'max_depth', 'min_samples_leaf', 'min_samples_split', 'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 'min_impurity_decrease', 'class_weight', 'max_samples'])

    all_combinations = []
    all_combinations.append(n_est_s)
    all_combinations.append(crit_s)
    all_combinations.append(max_depth)
    all_combinations.append(min_samples_leaf)
    all_combinations.append(min_samples_split)
    all_combinations.append(min_weight_fraction_leaf)
    all_combinations.append(max_features)
    all_combinations.append(max_leaf_nodes)
    all_combinations.append(min_impurity_decrease)
    all_combinations.append(class_weight)
    all_combinations.append(max_samples)

    tot = sum(map(len, all_combinations))
    run = 0
    for name, feat in zip(names,all_combinations):
        print("%d%% complete"%((run/tot)*100))
        run = run + len(feat)
        params = {name:feat}
        gcv = GridSearchCV(clf, params, scoring='accuracy')
        gcv.fit(X,Y)
        best_params[name]=gcv.best_params_[name]
        cv_results[name]=gcv.cv_results_
    
    return best_params, cv_results