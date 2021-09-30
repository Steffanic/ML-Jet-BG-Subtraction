import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import time
from sklearn.tree import DecisionTreeClassifier
from sklearn import tree
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.ensemble import  RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.metrics import silhouette_score, silhouette_samples
from sklearn.decomposition import PCA
from RFBDoSBC.GetAndPrepareData import *
from RFBDoSBC.modelPreparation import *
from RFBDoSBC.modelEvaluation import *
from RFBDoSBC.plotData import *
from RFBDoSBC.utility import *
from sklearn.preprocessing import MinMaxScaler
from matplotlib.ticker import ScalarFormatter, LogFormatter, LogLocator


matplotlib.rcParams.update({"font.size":12, "text.usetex": True})


feats = ['Area', 'p_T', 'p_T-corr', 'Mean-p_T', 'p_T_1', 'Angularity', 'Var-p_T']
feat_labels = ['Area', '$p_T$', '$p_T^{corr}$', '$\left<p_T\\right>$', '$p_T^1$', 'Angularity', '$\left<(\Delta p_T)^2\\right>$']
feat_units = ['', '[GeV]', '[GeV]', '[GeV]', '[GeV]', '', '[GeV]']

master_xticks = {'Area': [ 0.1,  0.45, 0.7]
                ,'p_T': [  50, 125, 200.]
                ,'p_T-corr': [-20., 35,  90]
                ,'Mean-p_T':  [0.5,1., 1.5 ]
                ,'p_T_1': [3,  10., 40]
                ,'Angularity': [ 0.15, 0.25, 0.35]
                ,'Var-p_T': [1, 5, 20]}
master_yticks = {'Area': [10e-4, 10e-2, 10e0]
                ,'p_T': [10e-6, 10e-4, 10e-2]
                ,'p_T-corr': [10e-6, 10e-4, 10e-2]
                ,'Mean-p_T': [10e-5, 10e-3, 10e-1]
                ,'p_T_1': [10e-6, 10e-4, 10e-2]
                ,'Angularity': [10e-4, 10e-2, 10e0]
                ,'Var-p_T': [10e-4, 10e-2, 10e0]}

masked_plots = [(0,0), (1,0), (2, 0),(3,0),(0,1), (1,1), (2, 1), (3, 1), (0,2), (0,3), (0,4)]

real_markers = ["v", "o", "8", "s", "p", "P"]
fake_markers = ["v", "o", "8", "s", "p", "P"]

for feat, label, unit in zip(feats, feat_labels, feat_units):
    i=0
    rad_mom_gen = res_pT_iterator()
    fig, ax = plt.subplots(1,5, figsize=(8, 1.5))
    fig.subplots_adjust(hspace=0, wspace=0)
    print(f"Making plot for {feat}")

    ext_ylim = [0.01, 1]
    try:
        while True:

            rad, ptm = next(rad_mom_gen)
            if rad in [] or ptm in []:
                continue
            msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
            dataGen = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm, True, False)

            batch_num = 0

            while batch_num<1:
                try:
                    print(f"Doing batch_num {batch_num}")
                    
                    train, test = next(dataGen)
                    print(f"{len(train)} samples")
                    X, Y = split_feat_label(train)

                    X = X.loc[Y!=2] # Remove the squishy jets.
                    Y = Y.loc[Y!=2]

                    print(X.columns)

                    #ss_vals = silhouette_samples(Xnorm, Y, n_jobs=-1)

                    fake_inds = true_indices_for_class(1, Y)
                    real_inds = true_indices_for_class(3, Y)

                    all_inds = np.concatenate([fake_inds, real_inds])

                    #real_ss_vals = ss_vals[real_inds]
                    #fake_ss_vals = ss_vals[fake_inds]
                    nr, binsr, patches=ax[i//6].hist(X[feat].iloc[real_inds], density=True, log=True, histtype="step", linewidth=0, ls='-')
                    nc, binsc, patches=ax[i//6].hist(X[feat].iloc[fake_inds], density=True, log=True, histtype="step", linewidth=0, ls='--')
                    ax[i//6].scatter(binsr[:-1], nr,  marker=real_markers[i%6], c='blue', linewidths=.5, s=10, alpha=1, label=f"Real {ptm} GeV/c")
                    ax[i//6].plot(binsr[:-1], nr, linewidth=.2, c='blue', alpha=1)
                    ax[i//6].scatter(binsc[:-1], nc, marker=fake_markers[i%6], c='none', edgecolors='red', linewidths=.5, s=10, alpha=1, label=f"Combinatorial {ptm} GeV/c")
                    ax[i//6].plot(binsc[:-1], nc, linewidth=.2, c='red', alpha=1)
                    #ax[i//2].fill_between(binsr[:-1], nr,color='blue', alpha=0.05)
                    #ax[i//2].fill_between(binsc[:-1], nc,color='red', alpha=0.05)
                    #print(ax[i//4, i%4].get_xticks(),ax[i//4, i%4].get_yticks())
                    ext_ylim[1] = ax[i//6].get_ylim()[1] if ax[i//6].get_ylim()[1]>ext_ylim[1] else ext_ylim[1]
                    ext_ylim[0] = ax[i//6].get_ylim()[0] if (ax[i//6].get_ylim()[0]<ext_ylim[0] and ax[i//6].get_ylim()[0]>0) else ext_ylim[0]
                    if feat=='p_T_1':
                        ax[i//6].set_xscale('log')
                        ax[i//6].set_yscale('log')
                        ax[i//6].set_xticks([t for t in master_xticks[feat] if t < ax[i//6].get_xlim()[1]])
                        ax[i//6].set_yticks(master_yticks[feat])
                        ax[i//6].xaxis.set_major_formatter(ScalarFormatter())
                        #ax[i//2].yaxis.set_major_locator(LogLocator())
                        #ax[i//2].yaxis.set_minor_locator(LogLocator())
                        #ax[i//2].yaxis.set_major_formatter(LogFormatter())
                        #ax[i//2].yaxis.set_minor_formatter(LogFormatter())


                    else:
                        ax[i//6].set_yscale('log')
                        ax[i//6].set_xticks([t for t in master_xticks[feat] if t < ax[i//6].get_xlim()[1]])
                        ax[i//6].set_yticks(master_yticks[feat])
                        #ax[i//2].yaxis.set_major_locator(LogLocator())
                        #ax[i//2].yaxis.set_minor_locator(LogLocator())
                        #ax[i//2].yaxis.set_major_formatter(LogFormatter())
                        #ax[i//2].yaxis.set_minor_formatter(LogFormatter())
                    ax[i//6].set_zorder(i//2)
                    #if i%2==0:
                    #    ax[i%2, 0].set_title(f"R={rad}")
                    '''
                    if i==0:
                        master_xticks = ax[i//4, i%4].get_xticks()
                        master_yticks = ax[i//4, i%4].get_yticks()
                    '''
                    
                    
                    
                    batch_num+=1
                except StopIteration:
                    print(f"No more data in batch num: {batch_num}, trimming forest.")
                    break
            i+=1

        
    except StopIteration:
        #ax[4,0]._shared_x_axes.remove(ax[4,0])
        #ax[4,0]._shared_y_axes.remove(ax[4,0])
        #ax[4,0].xaxis.set_tick_params(tick1On=True)
        #ax[4,0].yaxis.set_tick_params(tick1On=True)
        #[ax[i].set_xticklabels([]) for i in range(4) if i!=0]
        #ax[4,0].set_xticklabels(master_xticks)
        [ax[i].set_yticklabels([]) for i in range(5) if i!=0]
        [ax[i].set_ylim(ext_ylim) for i in range(5)]
        #ax[4,0].set_yticklabels(master_yticks)
        ax[0].set_ylabel("$\\frac{1}{N_{jets}} \\frac{dN_{jets}}{d%s}$"%(label if '$' not in label else "".join(label.split("$")).replace('$', '')))
        ax[0].set_xlabel(f"{label} {unit}", loc='center')
        ax[4].legend(bbox_to_anchor=(2.62, 0.7), loc='upper right', prop={'size':10})
        fig.subplots_adjust(top=0.995, bottom=0.34, right=0.75, left=0.125)
        #plt.tight_layout()
        #plt.show()
        plt.savefig(f"paper_figs/{feat}_dist.png", dpi=300)
