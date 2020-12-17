import numpy as np 
import pickle
import matplotlib.pyplot as plt
import matplotlib
import pandas
import os


rjr_train = []
fjr_train = []
sjr_train = []
rjr_test = []
fjr_test = []
sjr_test = []
fo_train = []
fi_train = []
fo_test = []
fi_test = []
matplotlib.rcParams.update({'font.size': 14})
matplotlib.rcParams.update({'font.family':'Liberation Serif'})
for fil in os.listdir("Objects/"):
    if fil[0]=='R':
        rjr_train = []
        fjr_train = []
        sjr_train = []
        rjr_test = []
        fjr_test = []
        sjr_test = []
        fo_train = []
        fi_train = []
        fo_test = []
        fi_test = []
        for ptm in [10,20,30,40]:
            with open("Objects/%s/perf_metric_pthard=%d"%(fil, ptm), "rb") as metfile:
                perf_metric = pickle.load(metfile)
                rjr_test.append(perf_metric['real_jet_rate_test'])
                fjr_test.append(perf_metric['fake_jet_rate_test'])
                sjr_test.append(perf_metric['squish_jet_rate_test'])
                rjr_train.append(perf_metric['real_jet_rate_train'])
                fjr_train.append(perf_metric['fake_jet_rate_train'])
                sjr_train.append(perf_metric['squish_jet_rate_train'])
                fo_test.append(perf_metric['fall_out_test'])
                fi_test.append(perf_metric['fall_in_test'])
                fo_train.append(perf_metric['fall_out_train'])
                fi_train.append(perf_metric['fall_in_train'])

        
        plt.plot([10,20,30,40], rjr_test, 'b.')
        plt.ylim(0.98,1)
        plt.title("Real Jet Rate for %s, test"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Real Jets Classified Real")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/RealJetRate_test_.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fjr_test, 'b.')
        plt.ylim(0.6,1)
        plt.title("Fake Jet Rate for %s, test"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Fake Jets Classified Fake")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/FakeJetRate_test_.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], sjr_test, 'b.')
        plt.ylim(0,1)
        plt.title("Squish Jet Rate for %s, test"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Squish Jets Classified Squish")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/SquishJetRate_test_.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], rjr_train, 'b.')
        plt.ylim(0.98,1)
        plt.title("Real Jet Rate for %s, train"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Real Jets Classified Real")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/RealJetRate_train.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fjr_train, 'b.')
        plt.ylim(0.6,1)
        plt.title("Fake Jet Rate for %s, train"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Fake Jets Classified Fake")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/FakeJetRate_train.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], sjr_train, 'b.')
        plt.ylim(0,1)
        plt.title("Squish Jet Rate for %s, train"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Squish Jets Classified Squish")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/SquishJetRate_train.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fo_test, 'b.')
        plt.ylim(0,1)
        plt.title("Percent Fake predicted Real for %s, test"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Fake Jets Classified Real")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/FakePredReal_test.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fi_test, 'b.')
        plt.ylim(0,1)
        plt.title("Percent Real predicted Fake for %s, test"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Real Jets Classified Fake")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/RealPredFake_test.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fo_train, 'b.')
        plt.ylim(0,1)
        plt.title("Percent Fake predicted Real for %s, train"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Fake Jets Classified Real")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/FakePredReal_train.png"%fil)
        plt.close()
        
        

        plt.plot([10,20,30,40], fo_train, 'b.')
        plt.ylim(0,1)
        plt.title("Percent Real predicted Fake for %s, train"%fil)
        plt.xlabel("p_T hard min[GeV]")
        plt.ylabel("% Real Jets Classified Fake")
        plt.grid()
        plt.tight_layout()
        plt.savefig("Plots/%s/RealPredFake_train.png"%fil)
        plt.close()
        

             


                