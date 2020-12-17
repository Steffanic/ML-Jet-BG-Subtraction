import matplotlib.pyplot as plt
import os
import numpy as np
import pickle
from matplotlib import cm, colors
from pandas.plotting import scatter_matrix
from scipy.stats import gaussian_kde

def msg(_msg):
    print("JEB: %s"%_msg)

def plot_everything(train, rad, ptm):
    msg("Plotting feature distributions.")
    plt.style.use("seaborn-dark")
    plot_all_columns(train, rad, ptm)

    msg("Plotting scatter matrix.")
    _, _ = factor_scatter_matrix(train,'Label', rad, ptm, ['blue', 'orange', 'red'])
    
    plt.style.use("classic")
    msg("Plotting correlation matrix.")
    plot_corr_mat(train, rad, ptm)

def plot_all_columns(train, rad, ptm):
    for cols in train.columns:
        plot_population(train, cols, rad, ptm)
        #plot_population_log(cols)

def plot_population(train, feat, rad, ptm):
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
        with open("Plots/Feature_Plots/R=%1.1f/%s_pTmin%d.pickle"%(rad, feat, ptm), 'wb') as fil:
            pickle.dump(fig, fil)
        plt.close()

def plot_population_log(train, feat):
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
    
    with open("Plots/Feature_Plots/R=%1.1f/corr_mat_pTmin%d.pickle"%(rad, ptm), 'wb') as fil:
        pickle.dump(fig, fil)
    plt.close(fig)

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
    if(not os.path.isdir("Plots/R=%1.1f"%rad)):
        os.mkdir("Plots/R=%1.1f"%rad)
    if(not os.path.isdir("Plots/Feature_Plots")):
        os.mkdir("Plots/Feature_Plots")
    plt.savefig("Plots/R=%1.1f/scatter_mat_pTmin%d.png"%(rad, ptm))
    with open("Plots/Feature_Plots/R=%1.1f/scatter_mat_pTmin%d.pickle"%(rad, ptm), 'wb') as fil:
        pickle.dump(axarr, fil)
    plt.close()
    return axarr, color_map

