// HEADER FILES WILL NEED TO BE FIRST

#ifndef __CINT__
#include "TFile.h"
#include "TClonesArray.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TRandom3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include "THnSparse.h"
#include "TObject.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TTimeStamp.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include <stdio.h>
#include <algorithm>
#include <unistd.h>
#include <iostream>
#include <numeric>
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair
#include <stdexcept> // std::runtime_error
#include <sstream> // std::stringstream
#include <iomanip> // std::setprecision
#endif



using namespace std;

// declare the function and define it
void Christine_Plot_Macro( const char * input_file_name="Background_And_Pythia_Test_3_HF=0_PJ_10_30_FF_test2_DATA.root" , const char * output_file_name="Christine.root" ){

  TFile * infile1 = TFile::Open(input_file_name);

  TTree *treePythia;// = (TTree*) file->Get("Pythia");
  TTree *treePythiaTenn;
  treePythiaTenn = (TTree*) infile1->Get("Pythia-and-TennGen-HF-0-CB-0");
  treePythia = (TTree*) infile1->Get("Pythia");


  treePythiaTenn->Draw("pTcorr","MeanpT>0.68");
  TCanvas *can = new TCanvas("can", "can");
  treePythia->Draw("pTcorr","MeanpT>0.68");


}
