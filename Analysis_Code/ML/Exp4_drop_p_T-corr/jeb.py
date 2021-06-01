import  pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import pickle
import matplotlib
import gc
from pandas_profiling import ProfileReport
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from GetAndPrepareData import get_data, add_feats, label_data, type_data, split_data, drop_feat, balance_classes

matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'font.family':'Liberation Serif'})

def DataPipeline(filepath, ptm):
    dat = get_data(filepath)#, rows_=100000)#../Generator Output/LOWSTATS/Rparam-=-0.2-pThardmin-=10.0/merged-ML-output-LOWSTATS-Rparam-0.2-pThardmin-10.0_MINI_00.csv")
    dat_typed = type_data(dat)
    dat_feat_added = add_feats(dat_typed)
    dat_labeled = label_data(dat_feat_added, 'pythia-mom', [-1,0.000000001, ptm, 10000000], [0,1,2])
    dat_drop = drop_feat(dat_labeled, ['Eta', 'Phi', 'p_T', 'Angularity-NW', 'N-Trk', 'p_T-corr', 'p_T_3', 'p_T_4', 'p_T_5', 'distmatch', 'XMatch'])
    dat_bal = balance_classes(dat_drop)
    train, test = split_data(dat_bal)
    return train, test

def split_feat_label(data, drop_labels=['X_tru', 'pythia-mom', 'Label'], label_='Label'):
    X = data.drop(drop_labels, 1)
    Y=data[label_]
    return X, Y

def plot_population(feat, rad, ptm):
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    fig=plt.figure()
    minx = train[feat].min()
    maxx = train[feat].max()
    train.hist(feat,bins=25, density=True,ax=plt.gca(), label="Total",histtype=u'step', linewidth=2, color='black', linestyle='dashed')
    train.loc[train['Label']==0].hist(feat,bins=25, density=True,ax=plt.gca(), label="Fake",histtype=u'step', linewidth=2, color='blue')
    train.loc[train['Label']==1].hist(feat,bins=25, density=True,ax=plt.gca(), label="0-10GeV Pythia",histtype=u'step', linewidth=2, color='orange')
    train.loc[train['Label']==2].hist(feat,bins=25, density=True,ax=plt.gca(), label=">10GeV Pythia",histtype=u'step', linewidth=2, color='red')
    plt.xlim(minx, maxx)
    plt.xlabel(feat, fontsize=14)
    plt.ylabel("Jet Yield", fontsize=14)
    plt.title(feat+", R=%1.1f, p_T hardmin=%d"%(rad, ptm), fontsize=14)
    plt.tight_layout()
    plt.legend(loc='best', fontsize=14)
    plt.grid(0)
    #plt.semilogy()
    with open("Plots/Feature_plots/R=%1.1f/%s_pTmin%d.pickle"%(rad, feat, ptm), 'wb') as fil:
        pickle.dump(fig, fil)
    plt.close()


def plot_population_log(feat):
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    minx = train[feat].min()
    maxx = train[feat].max()
    train.hist(feat, density=True,ax=plt.gca(), label="Total",histtype=u'step', linewidth=2, color='black')
    train.loc[train['Label']==0].hist(feat,bins=30, density=True,ax=plt.gca(), label="Fake",histtype=u'step', linewidth=2, color='blue')
    train.loc[train['Label']==1].hist(feat,bins=30, density=True,ax=plt.gca(), label="0-10GeV Pythia",histtype=u'step', linewidth=2, color='orange')
    train.loc[train['Label']==2].hist(feat,bins=30, density=True,ax=plt.gca(), label=">10GeV Pythia",histtype=u'step', linewidth=2, color='red')
    plt.xlim(minx, maxx)
    plt.legend(loc='best')
    plt.title(feat, fontsize=14, fontname="Liberation Serif")
    plt.xticks(fontsize=14, fontname="Liberation Serif")
    plt.grid(0)
    plt.semilogy()
    plt.close()


def factor_scatter_matrix(df, factor, rad, ptm,palette=None):
    '''Create a scatter matrix of the variables in df, with differently colored
    points depending on the value of df[factor].
    inputs:
        df: pandas.DataFrame containing the columns to be plotted, as well 
            as factor.
        factor: string or pandas.Series. The column indicating which group 
            each row belongs to.
        palette: A list of hex codes, at least as long as the number of groups.
            If omitted, a predefined palette will be used, but it only includes
            9 groups.
    '''
    fig=plt.figure(figsize=(15,15))
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    import matplotlib.colors
    from pandas.plotting import scatter_matrix
    from scipy.stats import gaussian_kde

    if isinstance(factor, str):
        factor_name = factor #save off the name
        factor = df[factor] #extract column
        df = df.drop(factor_name,axis=1) # remove from df, so it 
        # doesn't get a row and col in the plot.

    classes = list(set(factor))

    if palette is None:
        palette = ['#e41a1c', '#377eb8', '#4eae4b', 
                '#994fa1', '#ff8101', '#fdfc33', 
                '#a8572c', '#f482be', '#999999']

    color_map = dict(zip(classes,palette))

    if len(classes) > len(palette):
        raise ValueError('''Too many groups for the number of colors provided.
We only have {} colors in the palette, but you have {}
groups.'''.format(len(palette), len(classes)))
    colors = factor.apply(lambda group: color_map[group])
    axarr = scatter_matrix(df, ax = plt.gca(),grid=True,marker='o',s=4,c=colors,diagonal='hist', hist_kwds={'bins':30, 'histtype':u'step'})
    plt.suptitle("2-D Feature Scatter Matrix", fontsize=24)
    for ax in axarr.flatten():
        ax.set_xlabel(ax.get_xlabel(), fontsize=14, fontname='Liberation Serif')
        ax.set_ylabel(ax.get_ylabel(), fontsize=14, fontname='Liberation Serif')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(True, which='both')
        


    '''
    for rc in range(len(df.columns)):
        for group in classes:
            y = df[factor == group].iloc[rc].values
            gkde = gaussian_kde(y)
            ind = np.linspace(y.min(), y.max(), 1000)
            axarr[rc][rc].plot(ind, gkde.evaluate(ind),c=color_map[group], label=group)
    '''
    plt.text(350, 200,'Fake - Blue\nSquitch - Orange\nReal - Red', fontsize=16)
    
    #plt.show()
    plt.tight_layout()
    plt.savefig("Plots/R=%1.1f/scatter_mat_pTmin%d.png"%(rad, ptm))
    with open("Plots/Feature_plots/R=%1.1f/scatter_mat_pTmin%d.pickle"%(rad, ptm), 'wb') as fil:
        pickle.dump(axarr, fil)
    plt.close()
    return axarr, color_map

from matplotlib import cm, colors
def plot_corr_mat(train, rad, ptm):
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    fig,ax = plt.subplots(1,figsize=(15,15))
    ax.matshow(train.corr())
    cbar = fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=-1, vmax=1)), ax=ax)
    cols = train.columns
    cbar.ax.tick_params(labelsize=16)
    plt.xticks(np.linspace(0, len(cols)-2, len(cols)-1), cols.drop('Label'), rotation=30, fontsize=16)
    plt.yticks(np.linspace(0, len(cols)-2, len(cols)-1), cols.drop('Label'), rotation=30, fontsize=16)
    plt.ylim(-0.5, len(cols)-2+0.5)
    plt.xlim(-0.5, len(cols)-2+0.5)
    plt.title("Feature Correlations", pad=50, fontsize=20)
    #plt.tight_layout()
    
    with open("Plots/Feature_plots/R=%1.1f/corr_mat_pTmin%d.pickle"%(rad, ptm), 'wb') as fil:
        pickle.dump(fig, fil)
    plt.close(fig)

def DoGridSearch(X, Y):
    import itertools
    clf = RandomForestClassifier(oob_score=True, random_state=42, 
    n_jobs=-1)
    best_params = {}
    cv_results = {}

    n_est_s = (100,200) # Number of Estimators in the forest

    crit_s = ('gini', 'entropy') # Criterion for evaluating purity

    max_depth = (1, 2,3) # Number of decision layers available to the decision trees

    min_samples_split = (2, 4, 20, 40, 100) # Minimum number of samples required to split an internal node

    min_samples_leaf = (100, 1000, 2000) # Minimum number of samples allowed left after a split

    min_weight_fraction_leaf = (0.1, 0.2, 0.3, 0.45) # ??

    max_features = (1,2,3,4, 5, 6)

    max_leaf_nodes = (10, 100, 1000)

    min_impurity_decrease = (0.0, 0.3, 0.8)

    class_weight = ({0:1, 1:1, 2:10}, {0:1,1:1,2:100}, {0:1, 1:1, 2:1000})

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

    print(len(all_combinations))
    for name, feat in zip(names,all_combinations):
        params = {name:feat}
        print(params)
        gcv = GridSearchCV(clf, params, scoring='accuracy')
        gcv.fit(X,Y)
        best_params[name]=gcv.best_params_[name]
        cv_results[name]=gcv.cv_results_
    
    return best_params, cv_results




if __name__=='__main__':
    doGridSearch=False

    R = [0.2, 0.3, 0.4, 0.5, 0.6]
    pThardmin = [10,20,30,40]

    feat_imp = {}

    train = None
    test = None
    for rad in R:
        for ptm in pThardmin:
            if(train is not None):
                del train, test
                print("Deleting train, test from the previous run.")
            gc.collect()
            print("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
            train, test = DataPipeline("../Generator Output/merged-ML-output-LOWSTATS-Rparam-%1.1f-pThardmin-%d.0.csv"%(rad, ptm), ptm)


            print("Number of Features:",len(train.columns))
            print("Number of Fake Jets:", len(train.loc[train['Label']==0]))
            print("Number of Squishy Jets:", len(train.loc[train['Label']==1]))
            print("Number of Real Jets:", len(train.loc[train['Label']==2]))
            print(train.describe())

            X, Y = split_feat_label(train)

            #Plotting

            print("Plotting feature distributions.")
            plt.style.use("seaborn-dark")
            for cols in X.columns:
                plot_population(cols, rad, ptm)
                #plot_population_log(cols)

            print("Plotting scatter matrix.")
            axarr, color_map = factor_scatter_matrix(train,'Label', rad, ptm, ['blue', 'orange', 'red'])
            
            plt.style.use("classic")
            print("Plotting correlation matrix.")
            plot_corr_mat(train, rad, ptm)

            if doGridSearch:
                best_params, cv_results = DoGridSearch(X, Y)
                print("Best Params:\n",best_params, "\nCV Results:\n", cv_results)

                best_params["oob_score"]=True

                best_params['random_state']=42
                best_params['n_jobs']=-1
                with open("Objects/best_params.pickle", 'wb') as fil:
                    pickle.dump(best_params, fil)
                doGridSearch=False
            else:
                with open("Objects/best_params.pickle", 'rb') as f:
                    best_params = pickle.load(f)


            print("Initializing model and fitting to data.")
            clf = RandomForestClassifier(**best_params)
            clf.fit(X,Y)


            
            importances = clf.feature_importances_
            feat_names = X.columns
            std = np.std([tree.feature_importances_ for tree in clf.estimators_],axis=0)
            quant_75 = np.quantile([tree.feature_importances_ for tree in clf.estimators_],0.75, axis=0)
            indices = np.argsort(importances)[::-1]

            # Print the feature ranking
            print("Feature ranking:")

            for f in range(X.shape[1]):
                print("%d. %s (%f)" % (f + 1, X.columns[indices[f]], importances[indices[f]]))

            plt.figure()
            plt.title("Feature importances")
            plt.bar(range(X.shape[1]), importances[indices],
                    color="r", yerr=std[indices], align="center")
            plt.xticks(range(X.shape[1]), X.columns[indices])
            plt.xlim([-1, X.shape[1]])
            #plt.show()

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
                print("Feature Importances:")
                print(col, clf.estimators_[3].feature_importances_[i])
                i+=1
            '''
            
            #print(gcv.best_estimator_.get_params())
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
            print(y_pred_train)
            print(len(np.array(Y==2).nonzero()[0]))
            real_jet_ind_train = np.array(Y==2).nonzero()[0]
            fake_jet_ind_train = np.array(Y==0).nonzero()[0]
            squish_jet_ind_train = np.asarray(Y==1).nonzero()[0]
            pred_for_real_jets_train = y_pred_train[real_jet_ind_train]
            pred_for_fake_jets_train = y_pred_train[fake_jet_ind_train]
            pred_for_squish_jets_train = y_pred_train[squish_jet_ind_train] 

            print(len(real_jet_ind_train),  np.asarray(pred_for_real_jets_train==2).nonzero()[0])

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

            print("Real Jet Rate train: ", real_jet_rate_train, "\nFake Jet Rate train: ", fake_jet_rate_train, "\nSquish Jet Rate train: ", squish_jet_rate_train, "\nFake predicted real train: ", fall_out_train, "\nReal predicted fake train: ", fall_in_train)

            X_test, Y_test = split_feat_label(test)

            #Real Jet rate = real jets pred/real jets
            y_pred_test = np.array(clf.predict(X_test))
            print(y_pred_test)
            real_jet_ind_test = np.asarray(Y_test==2).nonzero()[0]
            fake_jet_ind_test = np.asarray(Y_test==0).nonzero()[0]
            squish_jet_ind_test = np.asarray(Y_test==1).nonzero()[0]
            pred_for_real_jets_test = y_pred_test[real_jet_ind_test]
            pred_for_fake_jets_test = y_pred_test[fake_jet_ind_test]
            pred_for_squish_jets_test = y_pred_test[squish_jet_ind_test] 

            print(len(real_jet_ind_test), len(np.asarray(pred_for_real_jets_test==2).nonzero()[0]))

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


            print("Real Jet Rate test: ", real_jet_rate_test, "\nFake Jet Rate test: ", fake_jet_rate_test, "\nSquish Jet Rate test: ", squish_jet_rate_test, "\nFake predicted real test: ", fall_out_test, "\nReal predicted fake test: ", fall_in_test)


            with open("Objects/R=%1.1f/perf_metric_pthard=%d"%(rad, ptm), 'wb') as perffile:
                pickle.dump(perf_metrics, perffile)

            #plt.show()
        
        with open("Objects/R=%1.1f/feat_imp.pickle"%rad, 'wb') as fil:
            pickle.dump(feat_imp, fil)

        plt.close('all')