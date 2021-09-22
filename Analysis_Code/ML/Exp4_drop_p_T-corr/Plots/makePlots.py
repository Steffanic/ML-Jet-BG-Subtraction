import pickle
import matplotlib.pyplot as plt
import numpy as np
import os



for plott in os.listdir("Plots/Feature_plots/R=0.2"):
    with open("Plots/Feature_plots/R=0.2/%s"%plott, 'rb') as f:
        print("Saving R=0.2 plots: %s"%plott)
        fig = pickle.load(f)
        plt.savefig("Plots/R=0.2/%s"%plott[:-7]+".png")
        plt.close()
        

for plott in os.listdir("Plots/Feature_plots/R=0.3"):
    with open("Plots/Feature_plots/R=0.3/%s"%plott, 'rb') as f:
        print("Saving R=0.3 plots: %s"%plott)
        fig = pickle.load(f)
        plt.savefig("Plots/R=0.3/%s"%plott[:-7]+".png")
        plt.close()


for plott in os.listdir("Plots/Feature_plots/R=0.4"):
    with open("Plots/Feature_plots/R=0.4/%s"%plott, 'rb') as f:
        print("Saving R=0.4 plots: %s"%plott)
        fig = pickle.load(f)
        plt.savefig("Plots/R=0.4/%s"%plott[:-7]+".png")
        plt.close()


for plott in os.listdir("Plots/Feature_plots/R=0.5"):
    with open("Plots/Feature_plots/R=0.5/%s"%plott, 'rb') as f:
        print("Saving R=0.5 plots: %s"%plott)
        fig = pickle.load(f)
        plt.savefig("Plots/R=0.5/%s"%plott[:-7]+".png")
        plt.close()


for plott in os.listdir("Plots/Feature_plots/R=0.6"):
    with open("Plots/Feature_plots/R=0.6/%s"%plott, 'rb') as f:
        print("Saving R=0.6 plots: %s"%plott)
        fig = pickle.load(f)
        plt.savefig("Plots/R=0.6/%s"%plott[:-7]+".png")
        plt.close()