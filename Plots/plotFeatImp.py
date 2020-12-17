import numpy as np 
import pickle
import matplotlib.pyplot as plt
import pandas
import os
import matplotlib

matplotlib.rcParams.update({'font.size': 16})
matplotlib.rcParams.update({'font.family':'DejaVu Sans'})

for fil in os.listdir("Objects/"):
    if fil[0]=='R':
        with open("Objects/%s/feat_imp.pickle"%fil, "rb") as f:
            fignum = int(fil[-1])-1
            print("Figure %d"%fignum)
            plt.figure(fignum, figsize=(10, 5))
            feat_imp = pickle.load(f)
            print(feat_imp)
            i=2
            plt.plot(np.linspace(-1, 20, 2), [0,0])
            #plt.figure(1)
            
            plt.errorbar(np.linspace(0.2*i, 10+len(feat_imp['pthard=10'][1])+0.2*i,len(feat_imp['pthard=10'][1])), feat_imp['pthard=10'][1], yerr=feat_imp['pthard=10'][2], label="pthardmin=10", linewidth=i*0.3)
            plt.xticks(np.linspace(0.5, 10+len(feat_imp['pthard=10'][1])+0.5,len(feat_imp['pthard=10'][1])), feat_imp['pthard=10'][0])
            plt.title(fil)
            plt.xlabel("Feature")
            plt.ylabel("Feature Importance")
            plt.tight_layout()
            plt.legend()
            i+=1
            #plt.figure(2)
            plt.errorbar(np.linspace(0.2*i, 10+len(feat_imp['pthard=20'][1])+0.2*i,len(feat_imp['pthard=20'][1])), feat_imp['pthard=20'][1], yerr=feat_imp['pthard=20'][2], label="pthardmin=20", linewidth=i*0.3)
            plt.xticks(np.linspace(0.5, 10+len(feat_imp['pthard=20'][1])+0.5,len(feat_imp['pthard=20'][1])), feat_imp['pthard=20'][0])
            plt.title(fil)
            plt.xlabel("Feature")
            plt.ylabel("Feature Importance")
            plt.tight_layout()
            plt.legend()
            i+=1

            #plt.figure(3)
            plt.errorbar(np.linspace(0.2*i, 10+len(feat_imp['pthard=30'][1])+0.2*i,len(feat_imp['pthard=30'][1])), feat_imp['pthard=30'][1], yerr=feat_imp['pthard=30'][2], label="pthardmin=30", linewidth=i*0.3)
            plt.xticks(np.linspace(0.5, 10+len(feat_imp['pthard=30'][1])+0.5,len(feat_imp['pthard=30'][1])), feat_imp['pthard=30'][0])
            plt.title(fil)
            plt.xlabel("Feature")
            plt.ylabel("Feature Importance")
            plt.tight_layout()
            plt.legend()
            i+=1

            #plt.figure(4)
            plt.errorbar(np.linspace(0.2*i, 10+len(feat_imp['pthard=40'][1])+0.2*i,len(feat_imp['pthard=40'][1])), feat_imp['pthard=40'][1], yerr=feat_imp['pthard=40'][2], label="pthardmin=40", linewidth=i*0.3)
            plt.xticks(np.linspace(0.5, 10+len(feat_imp['pthard=40'][1])+0.5,len(feat_imp['pthard=40'][1])), feat_imp['pthard=40'][0])
            plt.title(fil)
            plt.xlabel("Feature")
            plt.ylabel("Feature Importance")
            plt.tight_layout()
            plt.legend()

            i+=1
        
plt.figure(1)
plt.savefig("Plots/R=0.2/feat_imp.png")
plt.figure(2)
plt.savefig("Plots/R=0.3/feat_imp.png")
plt.figure(3)
plt.savefig("Plots/R=0.4/feat_imp.png")
plt.figure(4)
plt.savefig("Plots/R=0.5/feat_imp.png")
plt.figure(6)
plt.savefig("Plots/R=0.6/feat_imp.png")

