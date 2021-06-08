import numpy as np
import os
import pickle

def compute_performance_metrics(model, X, Y, X_test, Y_test, rad, ptm):
    
    encode_jet_labels = lambda label: 3 if label=="real" else (2 if label=="squish" else (1 if label=="fake" else None)) 

    get_jet_indices = lambda data, ind: np.array(data==encode_jet_labels(ind)).nonzero()[0] 

    perf_metrics = {}

    predictions_train = np.array(model.predict(X)) # Get model's predictions on training data
    predictions_test = np.array(model.predict(X_test)) # Get the model's prediction on the test data

    print(f"Number of Real Jets train: {len(np.array(Y==3).nonzero()[0])}")
    print(f"Number of Fake Jets train: {len(np.array(Y==1).nonzero()[0])}")
    print(f"Number of Squishy Jets train: {len(np.array(Y==2).nonzero()[0])}")

    real_jets_train_mask = np.array(Y==3).nonzero()[0] # Store locations of real/fake/squish jets in dataset for comparison to predictions
    fake_jets_train_mask = np.array(Y==1).nonzero()[0] # Basically class masks
    squish_jets_train_mask = np.asarray(Y==2).nonzero()[0] # These are lists of indices

    pred_for_real_jets_train = predictions_train[real_jets_train_mask] # Apply masks to the predictions
    pred_for_fake_jets_train = predictions_train[fake_jets_train_mask] # These are lists of predictions
    pred_for_squish_jets_train = predictions_train[squish_jets_train_mask] 

    real_jet_rate_train = float(len(np.asarray(pred_for_real_jets_train==3).nonzero()[0]))/float(len(real_jets_train_mask))
    fake_jet_rate_train = float(len(np.asarray(pred_for_fake_jets_train==1).nonzero()[0]))/float(len(fake_jets_train_mask))
    squish_jet_rate_train = float(len(np.asarray(pred_for_squish_jets_train==2).nonzero()[0]))/float(len(squish_jets_train_mask))

    fall_out_train = float(len(np.asarray(pred_for_fake_jets_train==3).nonzero()[0]))/float(len(fake_jets_train_mask))
    fall_in_train = float(len(np.asarray(pred_for_real_jets_train==1).nonzero()[0]))/float(len(real_jets_train_mask))

    perf_metrics['real_jet_rate_train'] = real_jet_rate_train
    perf_metrics['fake_jet_rate_train'] = fake_jet_rate_train
    perf_metrics['squish_jet_rate_train'] = squish_jet_rate_train
    perf_metrics['fall_out_train'] = fall_out_train
    perf_metrics['fall_in_train'] = fall_in_train

    print("Real Jet Rate train: "+str(real_jet_rate_train)+"\nFake Jet Rate train: "+str(fake_jet_rate_train)+"\nSquish Jet Rate train: "+str(squish_jet_rate_train)+"\nFake predicted real train: "+str(fall_out_train)+"\nReal predicted fake train: "+str(fall_in_train))
    log("Real Jet Rate train: "+str(real_jet_rate_train)+"\nFake Jet Rate train: "+str(fake_jet_rate_train)+"\nSquish Jet Rate train: "+str(squish_jet_rate_train)+"\nFake predicted real train: "+str(fall_out_train)+"\nReal predicted fake train: "+str(fall_in_train))

    #Real Jet rate = real jets pred/real jets
    
    print(f"Number of Real Jets test: {len(np.array(Y_test==3).nonzero()[0])}")
    print(f"Number of Fake Jets test: {len(np.array(Y_test==1).nonzero()[0])}")
    real_jet_ind_test = np.asarray(Y_test==3).nonzero()[0]
    fake_jet_ind_test = np.asarray(Y_test==1).nonzero()[0]
    squish_jet_ind_test = np.asarray(Y_test==2).nonzero()[0]
    pred_for_real_jets_test = predictions_test[real_jet_ind_test]
    pred_for_fake_jets_test = predictions_test[fake_jet_ind_test]
    pred_for_squish_jets_test = predictions_test[squish_jet_ind_test] 

    print(str(len(real_jet_ind_test))+' '+str(len(np.asarray(pred_for_real_jets_test==3).nonzero()[0])))

    real_jet_rate_test = float(len(np.asarray(pred_for_real_jets_test==3).nonzero()[0]))/float(len(real_jet_ind_test))
    fake_jet_rate_test = float(len(np.asarray(pred_for_fake_jets_test==1).nonzero()[0]))/float(len(fake_jet_ind_test))
    squish_jet_rate_test = float(len(np.asarray(pred_for_squish_jets_test==2).nonzero()[0]))/float(len(squish_jet_ind_test))

    fall_out_test = float(len(np.asarray(pred_for_fake_jets_test==3).nonzero()[0]))/float(len(fake_jet_ind_test))
    fall_in_test = float(len(np.asarray(pred_for_real_jets_test==1).nonzero()[0]))/float(len(real_jet_ind_test))


    perf_metrics['real_jet_rate_test'] = real_jet_rate_test
    perf_metrics['fake_jet_rate_test'] = fake_jet_rate_test
    perf_metrics['squish_jet_rate_test'] = squish_jet_rate_test
    perf_metrics['fall_out_test'] = fall_out_test
    perf_metrics['fall_in_test'] = fall_in_test


    print("Real Jet Rate test: "+str(real_jet_rate_test)+"\nFake Jet Rate test: "+str(fake_jet_rate_test)+"\nSquish Jet Rate test: "+str(squish_jet_rate_test)+"\nFake predicted real test: "+str(fall_out_test)+"\nReal predicted fake test: "+str(fall_in_test))
    log("Real Jet Rate test: "+str(real_jet_rate_test)+"\nFake Jet Rate test: "+str(fake_jet_rate_test)+"\nSquish Jet Rate test: "+str(squish_jet_rate_test)+"\nFake predicted real test: "+str(fall_out_test)+"\nReal predicted fake test: "+str(fall_in_test))


    if(not os.path.isdir("Objects/R=%1.1f"%rad)):
        os.mkdir("Objects/R=%1.1f"%rad)
    with open("Objects/R=%1.1f/perf_metric_pthard=%d"%(rad, ptm), 'wb') as perffile:
        pickle.dump(perf_metrics, perffile)