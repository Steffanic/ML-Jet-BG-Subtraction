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


matplotlib.rcParams.update({"font.size":14, "text.usetex": True})


feats = ['Area', 'p_T', 'p_T-corr', 'Mean-p_T', 'p_T_1', 'Angularity']
feat_labels = ['Area', '$p_T$', '$p_T$-corr', 'Mean-$p_T$', '$p_T^1$', 'Angularity']
feat_units = ['', '[GeV]', '[GeV]', '[GeV]', '[GeV]', '']

master_xticks = {'Area': [0.,  0.35, 0.65, 1.]
                ,'p_T': [  0., 50, 125, 200., 300.]
                ,'p_T-corr': [-40., 0., 40, 80, 120.]
                ,'Mean-p_T':  [0.3,1., 1.75, 2.5 ]
                ,'p_T_1': [0.1, 1,  10., 100]
                ,'Angularity': [0.,  0.15, 0.3, 0.45]}
master_yticks = {'Area': [10e-4, 10e-2, 10e0]
                ,'p_T': [10e-6, 10e-4, 10e-2]
                ,'p_T-corr': [10e-6, 10e-4, 10e-2]
                ,'Mean-p_T': [10e-5, 10e-3, 10e-1]
                ,'p_T_1': [10e-6, 10e-4, 10e-2]
                ,'Angularity': [10e-4, 10e-2, 10e0]}

masked_plots = [(0,0), (1,0), (2, 0),(3,0),(0,1), (1,1), (2, 1), (3, 1), (0,2), (0,3), (0,4)]

real_markers = ["v", "+"]
fake_markers = ["^", "x"]

for feat, label, unit in zip(feats, feat_labels, feat_units):
    i=0
    rad_mom_gen = res_pT_iterator()
    fig, ax = plt.subplots(2,4, figsize=(8, 1.5))
    fig.subplots_adjust(hspace=0, wspace=0)
    print(f"Making plot for {feat}")

    ext_ylim = [0.01, 1]
    try:
        while True:

            rad, ptm = next(rad_mom_gen)
            if rad in [0.6] or ptm in [10, 30]:
                continue
            msg("Analyzing jets with R=%1.1f and p_T hardmin=%d"%(rad, ptm))
            dataGen = DataPipelineBatch(f"Analysis_Code/Generator Output/merged-ML-output-LOWSTATS-Rparam-{rad}-pThardmin-{ptm}.0.csv", ptm, True, False)

            batch_num = 0

            while batch_num<1:
                try:
                    print(f"Doing batch_num {batch_num}")
                    
                    train, test = next(dataGen)
                    X, Y = split_feat_label(train)

                    X = X.loc[Y!=2] # Remove the squishy jets.
                    Y = Y.loc[Y!=2]

                    print(X.columns)
                    scaler = MinMaxScaler() 

                    Xnorm = pd.DataFrame(scaler.fit_transform(X))

                    ss_vals = silhouette_samples(Xnorm, Y, n_jobs=-1)

                    fake_inds = true_indices_for_class(1, Y)
                    real_inds = true_indices_for_class(3, Y)

                    all_inds = np.concatenate([fake_inds, real_inds])

                    real_ss_vals = ss_vals[real_inds]
                    fake_ss_vals = ss_vals[fake_inds]

                    x_max, x_min = X[feat].max(), X[feat].min()

                    n_bins = 20
                    bin_width = (x_max-x_min)/n_bins
                    bin_left = [x for x in np.linspace(x_min, x_max, n_bins)]
                    bin_right = [x+bin_width for x in np.linspace(x_min, x_max, n_bins)]
                    x_bin_contents_real = [[ss_val for f_val, ss_val in zip(X[feat].iloc[real_inds], real_ss_vals) if f_val<right_edge and f_val>left_edge] for left_edge, right_edge in zip(bin_left, bin_right)]
                    ax[i%2,i//2].hist2d(X[feat].iloc[real_inds], real_ss_vals,bins=(20, 20), cmap=plt.cm.Blues , alpha=0.75)
                    ax[i%2,i//2].hist2d(X[feat].iloc[fake_inds], fake_ss_vals,bins=(20, 20), cmap=plt.cm.Reds, alpha=0.75)
                    
                    #nr, binsr, patches=ax[i//2].hist(X[feat].iloc[real_inds], density=True, log=True, histtype="step", linewidth=0, ls='-')
                    #nc, binsc, patches=ax[i//2].hist(X[feat].iloc[fake_inds], density=True, log=True, histtype="step", linewidth=0, ls='--')
                    #ax[i//2].scatter(binsr[:-1], nr, marker=real_markers[i%2], c='blue', linewidths=.5, s=10, alpha=1, label=f"Real $p_T$ hard min={ptm}")
                    #ax[i//2].plot(binsr[:-1], nr, linewidth=.2, c='blue', alpha=1)
                    #ax[i//2].scatter(binsc[:-1], nc, marker=fake_markers[i%2], c='red', linewidths=.5, s=10, alpha=1, label=f"Combinatorial $p_T$ hard min={ptm}")
                    #ax[i//2].plot(binsc[:-1], nc, linewidth=.2, c='red', alpha=1)
                    #ax[i//2].fill_between(binsr[:-1], nr,color='blue', alpha=0.05)
                    #ax[i//2].fill_between(binsc[:-1], nc,color='red', alpha=0.05)
                    #print(ax[i//4, i%4].get_xticks(),ax[i//4, i%4].get_yticks())
                    #ext_ylim[1] = ax[i//2].get_ylim()[1] if ax[i//2].get_ylim()[1]>ext_ylim[1] else ext_ylim[1]
                    #ext_ylim[0] = ax[i//2].get_ylim()[0] if (ax[i//2].get_ylim()[0]<ext_ylim[0] and ax[i//2].get_ylim()[0]>0) else ext_ylim[0]
                    #if feat=='p_T_1':
                    #    ax[i//2].set_xscale('log')
                    #    ax[i//2].set_yscale('log')
                    #    ax[i//2].set_xticks(master_xticks[feat])
                    #    ax[i//2].set_yticks(master_yticks[feat])
                    #    ax[i//2].xaxis.set_major_formatter(ScalarFormatter())
                    #    #ax[i//2].yaxis.set_major_locator(LogLocator())
                    #    #ax[i//2].yaxis.set_minor_locator(LogLocator())
                    #    #ax[i//2].yaxis.set_major_formatter(LogFormatter())
                    #    #ax[i//2].yaxis.set_minor_formatter(LogFormatter())


                    #else:
                    #    ax[i//2].set_yscale('log')
                    #    ax[i//2].set_xticks(master_xticks[feat])
                    #    ax[i//2].set_yticks(master_yticks[feat])
                    #    #ax[i//2].yaxis.set_major_locator(LogLocator())
                    #    #ax[i//2].yaxis.set_minor_locator(LogLocator())
                    #    #ax[i//2].yaxis.set_major_formatter(LogFormatter())
                    #    #ax[i//2].yaxis.set_minor_formatter(LogFormatter())
                    #ax[i//2].set_zorder(i//2)
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
        #[ax[i].set_yticklabels([]) for i in range(4) if i!=0]
        #[ax[i].set_ylim(ext_ylim) for i in range(4)]
        #ax[4,0].set_yticklabels(master_yticks)
        #ax[0].set_ylabel("$\\frac{1}{N_{jets}} \\frac{dN_{jets}}{d%s}$"%(label if '$' not in label else "".join(label.split("$")).replace('$', '')))
        #ax[0].set_xlabel(f"{label} {unit}")
        #ax[3].legend(bbox_to_anchor=(2.62, 0.5), loc='upper right', prop={'size':8})
        fig.subplots_adjust(top=0.995, bottom=0.34, right=0.75, left=0.125)
        #plt.tight_layout()
        #plt.show()
        plt.savefig(f"paper_figs/{feat}_silhouettevalues.png", dpi=300)
