//This file creates the interface b/w root and pythia libraries for macros you want to run in root that use
//pythia libraries

void maindriver_newton(Int_t nEvents = 100000000, Int_t jobID =39034, Int_t tune = 100, Double_t Jet_Radius = 0.4, Int_t HF = 0 , Double_t DCA = 1000 , Int_t cent_bin = 0 , Int_t corrjet_bin = 1, Double_t constit_cut = 0.,  Bool_t Data = kTRUE, Bool_t RT_Stats = kFALSE ){

//________________load Pythia 6 and link to your library (libPythia6.s0)_________________________________________________//
gSystem->Load("libEGPythia6");
//gSystem->Load("/home/alidock/ML/pythia6/libPythia6.so"); //change to your setup
gSystem->Load("/home/charles/Documents/alice_software/pythia6/libPythia6.so"); //change to your setup
  gSystem->Load("$PYTHIA6428/libPythia6"); //change to your setup

//__________________________________________________________________________________________________//

//____//below lines needed on ACF !!!!__________________________________________________________//

  //gSystem->Load("libpythia6_4_28");
  //gSystem->Load("libEGPythia6");
  //gSystem->Load("libEG");
  //gSystem->Load("libAliPythia6");

//______________________________________________________________________________________________//

  gRandom->SetSeed(jobID); //needed for random number generatrion DO NOT FORGET

//__________________________________________________________________________________________________//

//____//below lines needed for FASTJET !!!!__________________________________________________________//

  //gSystem->Load("/nics/a/proj/UTK0019/alicesw/CGAL/lib/libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libfastjettools");
  gSystem->AddIncludePath("-I$FASTJET/include");
//  gSystem->AddIncludePath("-I/nics/c/home/jneuhau1/code/include");
//
  
//__________________________________________________________________________________________________//

//compile macro

//YOU MUST EXPORT THESE ENVIRONMENTAL VARIABLES TO RUN MY JET CODE

  gSystem->Load("libPhysics");
  gROOT->ProcessLine(".L BkgrLoad.cxx+");
  gROOT->ProcessLine(".L Toy_Model_ML_Study.C++");

  Toy_Model_ML_Study( nEvents, jobID , tune, Jet_Radius, HF , DCA, cent_bin , corrjet_bin, constit_cut, Data, RT_Stats) ;
}	
