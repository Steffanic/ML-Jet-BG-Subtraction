import matplotlib.pyplot as plt
import os
import numpy as np
import pickle
from matplotlib import cm, colors
from pandas.plotting import scatter_matrix
from scipy.stats import gaussian_kde
from sklearn.metrics import silhouette_samples

def msg(_msg):
    print("JEB: %s"%_msg)

def plot_everything(train, rad, ptm, batch_num, do_silhouette_score=False):
    msg("Plotting feature distributions.")
    plt.style.use("seaborn-dark")
    plot_all_columns(train, rad, ptm, batch_num)

    msg("Plotting scatter matrix.")
    _, _ = factor_scatter_matrix(train,'Label', rad, ptm, batch_num, ['blue', 'orange', 'red'])
    
    plt.style.use("classic")
    msg("Plotting correlation matrix.")
    plot_corr_mat(train, rad, ptm, batch_num)

    from RFBDoSBC.GetAndPrepareData import split_feat_label

    if do_silhouette_score:
        plt.style.use("seaborn-dark")
        X, y = split_feat_label(train)
        msg("Plotting silhouette distributions.")
        plot_silhouette_score_distributions(X, y, rad, ptm, batch_num)

def plot_all_columns(train, rad, ptm, batch_num):
    for cols in train.columns:
        plot_population(train, cols, rad, ptm, batch_num)
        #plot_population_log(cols)

def plot_population(train, feat, rad, ptm, batch_num):
    with np.errstate(divide='ignore', invalid='ignore'):
        plt.rc('xtick',labelsize=14)
        plt.rc('ytick',labelsize=14)
        plt.rc('font', family = "Liberation Serif")
        fig=plt.figure()
        minx = train[feat].min()
        maxx = train[feat].max()
        train.hist(feat,bins=25, density=True,ax=plt.gca(), label="Total",histtype=u'step', linewidth=2, color='black', linestyle='dashed')
        train.loc[train['Label']==1].hist(feat,bins=25, density=True,ax=plt.gca(), label="Fake",histtype=u'step', linewidth=2, color='blue')
        train.loc[train['Label']==2].hist(feat,bins=25, density=True,ax=plt.gca(), label="0-10GeV Pythia",histtype=u'step', linewidth=2, color='orange')
        train.loc[train['Label']==3].hist(feat,bins=25, density=True,ax=plt.gca(), label=">10GeV Pythia",histtype=u'step', linewidth=2, color='red')
        plt.xlim(minx, maxx)
        plt.xlabel(feat, fontsize=14)
        plt.ylabel("Jet Yield", fontsize=14)
        plt.title(feat+", R=%1.1f, p_T hardmin=%d"%(rad, ptm), fontsize=14)
        plt.tight_layout()
        plt.legend(loc='best', fontsize=14)
        plt.grid(0)
        #plt.semilogy()
        if(not os.path.isdir("Plots")):
            os.mkdir("Plots")
        if(not os.path.isdir("Plots/Feature_Plots")):
            os.mkdir("Plots/Feature_Plots")
        if(not os.path.isdir("Plots/Feature_Plots/R=%1.1f"%rad)):
            os.mkdir("Plots/Feature_Plots/R=%1.1f"%rad)
        with open(f"Plots/Feature_Plots/R={rad}/{feat}_pTmin{ptm}_batch{batch_num}.pickle", 'wb') as fil:
            pickle.dump(fig, fil)
        plt.close()

def plot_population_log(train, feat, batch_num):
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    minx = train[feat].min()
    maxx = train[feat].max()
    train.hist(feat, density=True,ax=plt.gca(), label="Total",histtype=u'step', linewidth=2, color='black')
    train.loc[train['Label']==1].hist(feat,bins=30, density=True,ax=plt.gca(), label="Fake",histtype=u'step', linewidth=2, color='blue')
    train.loc[train['Label']==2].hist(feat,bins=30, density=True,ax=plt.gca(), label="0-10GeV Pythia",histtype=u'step', linewidth=2, color='orange')
    train.loc[train['Label']==3].hist(feat,bins=30, density=True,ax=plt.gca(), label=">10GeV Pythia",histtype=u'step', linewidth=2, color='red')
    plt.xlim(minx, maxx)
    plt.legend(loc='best')
    plt.title(feat, fontsize=14, fontname="Liberation Serif")
    plt.xticks(fontsize=14, fontname="Liberation Serif")
    plt.grid(0)
    plt.semilogy()
    plt.close()

def plot_corr_mat(train, rad, ptm, batch_num):
    plt.rc('xtick',labelsize=14)
    plt.rc('ytick',labelsize=14)
    plt.rc('font', family = "Liberation Serif")
    fig,ax = plt.subplots(1,figsize=(15,15))
    ax.matshow(train.corr())
    cmap = cm.ScalarMappable(norm=colors.Normalize(vmin=-1, vmax=1))
    cmap.set_array([])
    cbar = fig.colorbar(cmap, ax=ax)
    cols = train.columns
    cbar.ax.tick_params(labelsize=16)
    plt.xticks(np.linspace(0, len(cols)-2, len(cols)-1), cols.drop('Label'), rotation=30, fontsize=16)
    plt.yticks(np.linspace(0, len(cols)-2, len(cols)-1), cols.drop('Label'), rotation=30, fontsize=16)
    plt.ylim(-0.5, len(cols)-2+0.5)
    plt.xlim(-0.5, len(cols)-2+0.5)
    plt.title("Feature Correlations", pad=50, fontsize=20)
    #plt.tight_layout()
    
    with open(f"Plots/Feature_Plots/R={rad}/corr_mat_pTmin{ptm}_batch{batch_num}.pickle", 'wb') as fil:
        pickle.dump(fig, fil)
    plt.close(fig)

def factor_scatter_matrix(df, factor, rad, ptm,batch_num, palette=None):
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
    

    factor_name = factor #save off the name
    factor = df[factor] #extract column
    df = df.drop(factor_name,axis=1) # remove from df, so it 
    # doesn't get a row and col in the plot.

    classes = [1,2,3]

    palette = ['blue', 'orange', 'red']

    color_map = dict(zip(classes,palette))

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
    if(not os.path.isdir("Plots/R=%1.1f"%rad)):
        os.mkdir("Plots/R=%1.1f"%rad)
    if(not os.path.isdir("Plots/Feature_Plots")):
        os.mkdir("Plots/Feature_Plots")
    plt.savefig("Plots/R=%1.1f/scatter_mat_pTmin%d.png"%(rad, ptm))
    with open(f"Plots/Feature_Plots/R={rad}/scatter_mat_pTmin{ptm}_batch{batch_num}.pickle", 'wb') as fil:
        pickle.dump(axarr, fil)
    plt.close()
    return axarr, color_map



def plot_silhouette_score_distributions(X, y, rad, ptm, batch_num):
    for col in X.columns:
        fig=plt.figure(figsize=(15,15))
        plt.rc('xtick',labelsize=24)
        plt.rc('ytick',labelsize=24)
        plt.rc('font', family = "Liberation Serif")
        plt.hist(silhouette_samples(np.array(X[col]).reshape((-1,1)), y))
        plt.title(f"Silhouette Score Distribution for {col}", fontsize=24)
        plt.xlabel("Silhouette Score", fontsize=24)
        plt.ylabel("Frequency", fontsize=24)    
        plt.tight_layout()

        
        if(not os.path.isdir("Plots")):
            os.mkdir("Plots")
        if(not os.path.isdir("Plots/Feature_Plots")):
            os.mkdir("Plots/Feature_Plots")
        if(not os.path.isdir("Plots/Feature_Plots/R=%1.1f"%rad)):
            os.mkdir("Plots/Feature_Plots/R=%1.1f"%rad)
        with open(f"Plots/Feature_Plots/R={rad}/Silhouette_Score_{col}_pTmin{ptm}_batch{batch_num}.pickle", 'wb') as fil:
            pickle.dump(fig, fil)
        plt.close()