// Author: Christian Holm Christensen <cholm@hilux15.nbi.dk>
// Update: 2002-08-16 16:40:27+0200
// Copyright: 2002 (C) Christian Holm Christensen
// Copyright (C) 2006, Rene Brun and Fons Rademakers.
// Copyright (c) 2017; Charles Hughes
// For the licensing terms see $ROOTSYS/LICENSE.
//
// 2nd Author: Richard Corke, Copyright 2013
// Roughly follows analysis of:
// * T. Aaltonen et al. [CDF Collaboration],
//
// 3rd Author: Joel Mazer <jmazer@utk.edu> (UTK collaboration)
//
// MODIFIED BY CHARLES HUGHES <chughe26@vols.utk.edu>
// University of Tennessee ALICE collaboration
// LAST UPDATE: 2019-10-04 18:45 EST
//


//This code is currently intended to be the framework for a jet finding code that provides
//fragmentation functions for final state particles in p-p (and eventually Pb-Pb, Au-Au, or others)
//collisions. This code will eventually
// 1) Find fragmentation functions for different partons (final state)
// 2) make pT threshold cuts
// 3) add detector effects (track reconstruction efficiency, pT resolution), add switches
// 4) calculate resolution matrix and unfold 
// 5) calculate raw spectrum of particles using dihadron correlations
// 6) apply unfolding and show result is the same as dihadron correlations 
// 7) add a realistic background and calculated corrected jets and frag functions
// 8) add reaction plane dependence and compare 
// 9) experiment with different energy loss models 

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
#include "TMCParticle.h"
#include "TObject.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/AreaDefinition.hh"
//#include "Pythia8/Pythia.h"
#include "TGraph.h"
#include "fastjet/Selector.hh"
//#include "BackgroundGenerator.h"
#include "BkgrLoad.h"
//#include "JetMatch.h"
#include "TLegend.h"
#include "TTimeStamp.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TF2.h"
#include <stdio.h>
#include <algorithm>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <numeric>

using namespace fastjet;

#include "TPythia6.h"
#endif
class TPythia6;
		

//class BackgroundGenerator;

//#include "TPythia8.h"
//class TPythia8;
int HEADER__=1;

// nEvents is how many events we want. This is dynamic user input to be run from the command prompt and user specified
// jobID is random number for seeding
// Jet_Radius is desired resolution paramater for jets
// HF is harmonic flag of background generator (only use this when generating background events, you don't use it when pre-loading)
// DCA - distance to closest approach (can use to discriminate pythia primaries from other particles)
// cent_bin - again for generating background events, don't need it if already pre-loading
// corrjet_bin - setting the bias for pythia
// constit_cut - for pythia + bkgd jets
// Data - file naming
// RT_Stats - use to get a feel for how long events take
void Toy_Model_ML_Study(Int_t nEvents, Int_t jobID , Int_t tune, Double_t Jet_Radius, Int_t HF , Double_t DCA , Int_t cent_bin , Int_t corrjet_bin , Double_t constit_cut, Bool_t Data = kTRUE, Bool_t RT_Stats = kFALSE ){

  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

// corrjet_bin: (1 : 10 - 30 GeV , pT_hard_min_pythia = 9 GeV) , (2 : 30 - 50 GeV , pT_hard_min_pythia = 29 GeV) , (3: 50 - 70 GeV , pT_hard_min_pythia = 49 GeV)

  char expression1[256];

  Double_t jet_pT_cut_low;
  Double_t jet_pT_cut_high;

  if( corrjet_bin == 1 ){

    if(Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_10_30_FF_test2_DATA.root" , HF);
    }
    else if(!Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_10_30_FF_test2_RESPONSE.root" , HF);
    }

    jet_pT_cut_low = 10.;
    jet_pT_cut_high = 30.;
    cout<<"\n\nlower pT bound for cut is "<< jet_pT_cut_low <<" GeV/c\n\n" << endl;
    cout<<"\n\nupper pT bound for cut is "<< jet_pT_cut_high <<" GeV/c\n\n" << endl;
  }

  else if( corrjet_bin == 2){

    if(Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_30_50_FF_test2_DATA.root" , HF);
    }
    else if(!Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_30_50_FF_test2_RESPONSE.root" , HF);
    }

    jet_pT_cut_low = 30.;
    jet_pT_cut_high = 50.;
    cout<<"\n\nlower pT bound for cut is "<< jet_pT_cut_low <<" GeV/c\n\n" << endl;
    cout<<"\n\nupper pT bound for cut is "<< jet_pT_cut_high <<" GeV/c\n\n" << endl;

  }

  else if( corrjet_bin == 3){

    if(Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_50_70_FF_test2_DATA.root" , HF);
    }
    else if(!Data){
      sprintf(expression1 , "Background_And_Pythia_Test_3_HF=%d_PJ_50_70_FF_test2_RESPONSE.root" , HF);
    }

    jet_pT_cut_low = 50.;
    jet_pT_cut_high = 70.;
    cout<<"\n\nlower pT bound for cut is "<< jet_pT_cut_low <<" GeV/c\n\n" << endl;
    cout<<"\n\nupper pT bound for cut is "<< jet_pT_cut_high <<" GeV/c\n\n" << endl;

  }


  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(0);

  //_____________________________________________________Creating Plots and Titles___________________________________________//
  //_____________________________________________________PYHTIA ONLY PLOTS FIRST___________________________________________//

  TH1D *histpythia_events = new TH1D( "histpythia_events" , "Number of Pythia Events" , 1 , 0.5 , 1.5 );
  TH1D *histpythia_jets = new TH1D( "histpythia_jets" , "Number of Pythia Jets" , 1 , 0.5 , 1.5 );

  char expression2[256];
  char expression3[256];
  char expression4[256];
  char expression5[256];
  char expression6[256];
  char expression7[256];
  char expression8[256];
  char expression9[256];
  char expression10[256];
  char expression11[256];
  char expression12[256];
  char expression13[256];
  char expression14[256];
  char expression15[256];

  sprintf(expression2, "p_{T} distribution of all anti-k_{T} Pythia Particles with |#eta| < 0.9" );
  sprintf(expression3, "#eta distribution of all anti-k_{T} Pythia Particles with |#eta| < 0.9" );
  sprintf(expression4, "#phi distribution of all anti-k_{T} Pythia Particles with |#eta| < 0.9" );

  sprintf(expression5, "p_{T} distribution of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius  );
  sprintf(expression6, "#eta distribution of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius);
  sprintf(expression7, "#phi distribution of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );

  sprintf(expression8, "Number of tracks for all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius  );
  sprintf(expression9, "Mean track p_{T} (<p_{T}>) of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius);
  sprintf(expression10, "Jet Angularity of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );

  sprintf(expression11, "Leading Track p_{T} of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );
  sprintf(expression12, "Sub-Leading Track p_{T} of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );
  sprintf(expression13, "3rd Leading Track p_{T} of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );

  sprintf(expression14, "Raw Fragmentation Function of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );
  sprintf(expression15, "LeSub of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );


  TH1D *histpT_pyth_part = new TH1D("histpT_pyth_part", expression2,100,0.,100.);
  histpT_pyth_part -> Sumw2();
  histpT_pyth_part->SetXTitle("p_{T}^{track} (GeV/c)");
  histpT_pyth_part->SetYTitle("dN/dp_{T}");

  TH1D *histeta_pyth_part = new TH1D("histeta_pyth_part", expression3,100,-0.9,0.9);
  histeta_pyth_part -> Sumw2();
  histeta_pyth_part->SetXTitle("#eta^{track} (pseudo-rap units)");
  histeta_pyth_part->SetYTitle("dN/d#eta");

  TH1D *histphi_pyth_part = new TH1D("histphi_pyth_part", expression4,100,0.,2.*TMath::Pi());
  histphi_pyth_part -> Sumw2();
  histphi_pyth_part->SetXTitle("#phi^{track} (radians)");
  histphi_pyth_part->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histpT_pyth_jet = new TH1D("histpT_pyth_jet", expression5,100,0.,100.);
  histpT_pyth_jet -> Sumw2();
  histpT_pyth_jet->SetXTitle("p_{T}^{jet} (GeV/c)");
  histpT_pyth_jet->SetYTitle("dN/dp_{T} ");

  TH1D *histeta_pyth_jet = new TH1D("histeta_pyth_jet", expression6,100,-0.9,0.9);
  histeta_pyth_jet -> Sumw2();
  histeta_pyth_jet->SetXTitle("#eta^{jet} (pseudo-rap units)");
  histeta_pyth_jet->SetYTitle("dN/d#eta");

  TH1D *histphi_pyth_jet = new TH1D("histphi_pyth_jet", expression7,100,0.,2.*TMath::Pi());
  histphi_pyth_jet -> Sumw2();
  histphi_pyth_jet->SetXTitle("#phi^{jet} (radians)");
  histphi_pyth_jet->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histnumtr_pyth_jet = new TH1D("histnumtr_pyth_jet", expression8,500,0,500);
  histnumtr_pyth_jet -> Sumw2();
  histnumtr_pyth_jet->SetXTitle("N_{tracks}");
  histnumtr_pyth_jet->SetYTitle("dN/dN_{tracks}");

  TH1D *histmeanpTtr_pyth_jet = new TH1D("histmeanpTtr_pyth_jet", expression9,200,0,100);
  histmeanpTtr_pyth_jet -> Sumw2();
  histmeanpTtr_pyth_jet->SetXTitle("<p_{T}^{track}> (GeV/c)");
  histmeanpTtr_pyth_jet->SetYTitle("dN/d<p_{T}^{track}>");

  TH1D *histangularity_pyth_jet = new TH1D("histangularity_pyth_jet", expression10,100,0.,10.);
  histangularity_pyth_jet -> Sumw2();
  histangularity_pyth_jet->SetXTitle("g (#phi - #eta distance units)");
  histangularity_pyth_jet->SetYTitle("dN/dg");

  TH1D *histldtr_pyth_jet = new TH1D("histldtr_pyth_jet", expression11,200,0.,100.);
  histldtr_pyth_jet -> Sumw2();
  histldtr_pyth_jet->SetXTitle("p_{T}^{ld. tr.} (GeV/c)");
  histldtr_pyth_jet->SetYTitle("dN/dp_{T}^{ld. tr.} ");

  TH1D *histsubldtr_pyth_jet = new TH1D("histsubldtr_pyth_jet", expression12,200,0.,100.);
  histsubldtr_pyth_jet -> Sumw2();
  histsubldtr_pyth_jet->SetXTitle("p_{T}^{sub-ld. tr.} (GeV/c)");
  histsubldtr_pyth_jet->SetYTitle("dN/dp_{T}^{sub-ld. tr.}");

  TH1D *histthrdldtr_pyth_jet = new TH1D("histthrdldtr_pyth_jet", expression13,200,0.,100.);
  histthrdldtr_pyth_jet -> Sumw2();
  histthrdldtr_pyth_jet->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histthrdldtr_pyth_jet->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histFF_pyth_jet = new TH1D("histFF_pyth_jet", expression14,100,0.,1.);
  histFF_pyth_jet -> Sumw2();
  histFF_pyth_jet->SetXTitle("z");
  histFF_pyth_jet->SetYTitle("dN/dz");

  TH1D *histlesub_pyth_jet = new TH1D("histlesub_pyth_jet", expression15,200,0.,100.);
  histlesub_pyth_jet -> Sumw2();
  histlesub_pyth_jet->SetXTitle("LeSub = p_{T}^{ld. tr.} - p_{T}^{sub-ld. tr.} (GeV/c)");
  histlesub_pyth_jet->SetYTitle("dN/dLeSub");

//___________________________________BACKGROUND ONLY PARTICLES_________________________________//

  TH1D *histbackground_events = new TH1D( "histbackground_events" , "Number of background Events" , 1 , 0.5 , 1.5 );
  TH1D *histbackground_jets = new TH1D( "histbackground_jets" , "Number of background Jets" , 1 , 0.5 , 1.5 );

  char expression16[256];
  char expression17[256];
  char expression18[256];
  char expression19[256];
  char expression20[256];
  char expression21[256];
  char expression22[256];
  char expression23[256];
  char expression24[256];
  char expression25[256];
  char expression26[256];
  char expression27[256];
  char expression28[256];
  char expression29[256];

  char expression30[256];
  char expression31[256];
  char expression32[256];
  char expression33[256];
  char expression34[256];

  char expression35[256];

  sprintf(expression16, "p_{T} distribution of all k_{T} Background Particles with |#eta| < 0.9" );
  sprintf(expression17, "#eta distribution of all k_{T} Background Particles with |#eta| < 0.9" );
  sprintf(expression18, "#phi distribution of all k_{T} Background Particles with |#eta| < 0.9" );

  sprintf(expression19, "p_{T} distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression20, "#eta distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression21, "#phi distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression22, "Number of tracks for all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression23, "Mean track p_{T} (<p_{T}>) of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression24, "Jet Angularity of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression25, "Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression26, "Sub-Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression27, "3rd Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression28, "Raw Fragmentation Function of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29, "LeSub of all k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression30, "Reconstructed #Psi_{1}^{EP} distribution for background particles");
  sprintf(expression31, "Reconstructed #Psi_{2}^{EP} distribution for background particles");
  sprintf(expression32, "Reconstructed #Psi_{3}^{EP} distribution for background particles");
  sprintf(expression33, "Reconstructed #Psi_{4}^{EP} distribution for background particles");
  sprintf(expression34, "Reconstructed #Psi_{5}^{EP} distribution for background particles");

  sprintf(expression35, "Median Background only k_{T} jets, #rho event by event");

//_______________________________________________________________________________________//

  TH1D *histpT_bkgd_part = new TH1D("histpT_bkgd_part", expression16,200,0.,10.);
  histpT_bkgd_part -> Sumw2();
  histpT_bkgd_part->SetXTitle("p_{T}^{tracks} (GeV/c)");
  histpT_bkgd_part->SetYTitle("dN/dp_{T} ");

  TH1D *histeta_bkgd_part = new TH1D("histeta_bkgd_part", expression17,100,-0.9,0.9);
  histeta_bkgd_part -> Sumw2();
  histeta_bkgd_part->SetXTitle("#eta^{tracks} (pseudo-rap units)");
  histeta_bkgd_part->SetYTitle("dN/d#eta");

  TH1D *histphi_bkgd_part = new TH1D("histphi_bkgd_part", expression18,100,0.,2.*TMath::Pi());
  histphi_bkgd_part -> Sumw2();
  histphi_bkgd_part->SetXTitle("#phi^{tracks} (radians)");
  histphi_bkgd_part->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histpT_bkgd_jet = new TH1D("histpT_bkgd_jet", expression19,100,0.,100.);
  histpT_bkgd_jet -> Sumw2();
  histpT_bkgd_jet->SetXTitle("p_{T} (GeV/c)");
  histpT_bkgd_jet->SetYTitle("dN/dp_{T} (GeV/c)");

  TH1D *histeta_bkgd_jet = new TH1D("histeta_bkgd_jet", expression20,100,-0.9,0.9);
  histeta_bkgd_jet -> Sumw2();
  histeta_bkgd_jet->SetXTitle("#eta^{jet} (pseudo-rap units)");
  histeta_bkgd_jet->SetYTitle("dN/d#eta");

  TH1D *histphi_bkgd_jet = new TH1D("histphi_bkgd_jet", expression21,100,0.,2.*TMath::Pi());
  histphi_bkgd_jet -> Sumw2();
  histphi_bkgd_jet->SetXTitle("#phi)^{jet} (radians)");
  histphi_bkgd_jet->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histnumtr_bkgd_jet = new TH1D("histnumtr_bkgd_jet", expression22,500,0,500);
  histnumtr_bkgd_jet -> Sumw2();
  histnumtr_bkgd_jet->SetXTitle("N_{tracks}");
  histnumtr_bkgd_jet->SetYTitle("dN/dN_{tracks}");

  TH1D *histmeanpTtr_bkgd_jet = new TH1D("histmeanpTtr_bkgd_jet", expression23,100,0,10);
  histmeanpTtr_bkgd_jet -> Sumw2();
  histmeanpTtr_bkgd_jet->SetXTitle("<p_{T}^{track}> (GeV/c)");
  histmeanpTtr_bkgd_jet->SetYTitle("dN/d<p_{T}^{track}>");

  TH1D *histangularity_bkgd_jet = new TH1D("histangularity_bkgd_jet", expression24,100,0.,10.);
  histangularity_bkgd_jet -> Sumw2();
  histangularity_bkgd_jet->SetXTitle("g (#phi - #eta distance units)");
  histangularity_bkgd_jet->SetYTitle("dN/dg");

  TH1D *histldtr_bkgd_jet = new TH1D("histldtr_bkgd_jet", expression25,100,0.,10.);
  histldtr_bkgd_jet -> Sumw2();
  histldtr_bkgd_jet->SetXTitle("p_{T}^{ld. tr.} (GeV/c)");
  histldtr_bkgd_jet->SetYTitle("dN/dp_{T}^{ld. tr.} ");

  TH1D *histsubldtr_bkgd_jet = new TH1D("histsubldtr_bkgd_jet", expression26,100,0.,10.);
  histsubldtr_bkgd_jet -> Sumw2();
  histsubldtr_bkgd_jet->SetXTitle("p_{T}^{sub-ld. tr.} (GeV/c)");
  histsubldtr_bkgd_jet->SetYTitle("dN/dp_{T}^{sub-ld. tr.}");

  TH1D *histthrdldtr_bkgd_jet = new TH1D("histthrdldtr_bkgd_jet", expression27,100,0.,10.);
  histthrdldtr_bkgd_jet -> Sumw2();
  histthrdldtr_bkgd_jet->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histthrdldtr_bkgd_jet->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histFF_bkgd_jet = new TH1D("histFF_bkgd_jet", expression28,100,0.,1.);
  histFF_bkgd_jet -> Sumw2();
  histFF_bkgd_jet->SetXTitle("z");
  histFF_bkgd_jet->SetYTitle("dN/dz");

  TH1D *histlesub_bkgd_jet = new TH1D("histlesub_bkgd_jet", expression29,100,0.,10.);
  histlesub_bkgd_jet -> Sumw2();
  histlesub_bkgd_jet->SetXTitle("LeSub = p_{T}^{ld. tr.} - p_{T}^{sub-ld. tr.} (GeV/c)");
  histlesub_bkgd_jet->SetYTitle("dN/dLeSub");

//____________________________EVENT PLANES___________________________________________________________//

  TH1D *histPsi_1 = new TH1D("histPsi_1", expression30,200,-0.5*TMath::Pi(),2.5*TMath::Pi());
  histPsi_1 -> Sumw2();
  histPsi_1->SetXTitle("#Psi_{1}^{reco}");
  histPsi_1->SetYTitle("dN/d#Psi_{1}^{reco}");

  TH1D *histPsi_2 = new TH1D("histPsi_2", expression31,200,-0.5*TMath::Pi(),2.5*TMath::Pi());
  histPsi_2 -> Sumw2();
  histPsi_2->SetXTitle("#Psi_{2}^{reco}");
  histPsi_2->SetYTitle("dN/d#Psi_{2}^{reco}");

  TH1D *histPsi_3 = new TH1D("histPsi_3", expression32,200,-0.5*TMath::Pi(),2.5*TMath::Pi());
  histPsi_3 -> Sumw2();
  histPsi_3->SetXTitle("#Psi_{3}^{reco}");
  histPsi_3->SetYTitle("dN/d#Psi_{3}^{reco}");

  TH1D *histPsi_4 = new TH1D("histPsi_4", expression33,200,-0.5*TMath::Pi(),2.5*TMath::Pi());
  histPsi_4 -> Sumw2();
  histPsi_4->SetXTitle("#Psi_{4}^{reco}");
  histPsi_4->SetYTitle("dN/d#Psi_{4}^{reco}");

  TH1D *histPsi_5 = new TH1D("histPsi_5", expression34,200,-0.5*TMath::Pi(),2.5*TMath::Pi());
  histPsi_5 -> Sumw2();
  histPsi_5->SetXTitle("#Psi_{5}^{reco}");
  histPsi_5->SetYTitle("dN/d#Psi_{5}^{reco}");

//______________________________________________________________________________________________//

  TH1D *histmedianrhoebe_bgkd = new TH1D("histmedianrhoebe_bgkd",expression35,600,0.,600.);
  histmedianrhoebe_bgkd -> Sumw2();
  histmedianrhoebe_bgkd->SetXTitle("#rho_{median}^{e.b.e.}");
  histmedianrhoebe_bgkd->SetYTitle("dN/d#rho_{median}^{e.b.e.}");


//___________________________________BACKGROUND + PYTHIA PARTICLES_________________________________//

  TH1D *histpythia_AND_bkgd_jets = new TH1D( "histpythia_AND_bkgd_jets" , "Number of Pythia & Jets" , 1 , 0.5 , 1.5 );

  char expression36[256];
  char expression37[256];
  char expression38[256];
  char expression39[256];
  char expression40[256];
  char expression41[256];
  char expression42[256];
  char expression43[256];
  char expression44[256];
  char expression45[256];
  char expression46[256];
  char expression47[256];
  char expression48[256];
  char expression49[256];

  char expression50[256];
  char expression51[256];

  char expression52[256];
  char expression53[256];

  sprintf(expression36, "p_{T} distribution of all anti-k_{T} Background + Pythia Particles with |#eta| < 0.9" );
  sprintf(expression37, "#eta distribution of all anti-k_{T} Background + Pythia Particles with |#eta| < 0.9" );
  sprintf(expression38, "#phi distribution of all anti-k_{T} Background + Pythia Particles with |#eta| < 0.9" );

  sprintf(expression39, "p_{T} distribution of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius  );
  sprintf(expression40, "#eta distribution of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius);
  sprintf(expression41, "#phi distribution of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );


  sprintf(expression42, "Jet p_{T} (area-based corr) for all anti-k_{T} Backgrouns + Pytha Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression43, "Number of tracks for all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius  );
  sprintf(expression44, "Mean track p_{T} (<p_{T}>) of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100 , 0.9 - Jet_Radius);
  sprintf(expression45, "Jet Angularity of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );

  sprintf(expression46, "Leading Track p_{T} of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression47, "Sub-Leading Track p_{T} of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression48, "3rd Leading Track p_{T} of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );

  sprintf(expression49, "Raw Fragmentation Function of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression50, "LeSub of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );

  sprintf(expression51, "Median Background + Particles (exclude 2 leading anti-k_{T} jets) only #rho event by event");

  sprintf(expression52, "Fraction of Pythia (true) tracks for all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius);

  sprintf(expression53, "Fraction of Background (fake) tracks for all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius);

//___________________________________________________________________________________________//

  TH1D *histpT_pytha_AND_bkgd_part = new TH1D("histpT_pytha_AND_bkgd_part", expression36,100,0.,100.);
  histpT_pytha_AND_bkgd_part -> Sumw2();
  histpT_pytha_AND_bkgd_part->SetXTitle("p_{T}^{tracks} (GeV/c)");
  histpT_pytha_AND_bkgd_part->SetYTitle("dN/dp_{T}");

  TH1D *histeta_pytha_AND_bkgd_part = new TH1D("histeta_pytha_AND_bkgd_part", expression37,100,-0.9,0.9);
  histeta_pytha_AND_bkgd_part -> Sumw2();
  histeta_pytha_AND_bkgd_part->SetXTitle("#eta^{track} (pseudo-rap units)");
  histeta_pytha_AND_bkgd_part->SetYTitle("dN/d#eta");

  TH1D *histphi_pytha_AND_bkgd_part = new TH1D("histphi_pytha_AND_bkgd_part", expression38,100,0.,2.*TMath::Pi());
  histphi_pytha_AND_bkgd_part -> Sumw2();
  histphi_pytha_AND_bkgd_part->SetXTitle("#phi^{track} (radians)");
  histphi_pytha_AND_bkgd_part->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histpT_pytha_AND_bkgd_jet = new TH1D("histpT_pytha_AND_bkgd_jet", expression39,100,0.,100.);
  histpT_pytha_AND_bkgd_jet -> Sumw2();
  histpT_pytha_AND_bkgd_jet->SetXTitle("p_{T}^{jet} (GeV/c)");
  histpT_pytha_AND_bkgd_jet->SetYTitle("dN/dp_{T} ");

  TH1D *histeta_pytha_AND_bkgd_jet = new TH1D("histeta_pytha_AND_bkgd_jet", expression40,100,-0.9,0.9);
  histeta_pytha_AND_bkgd_jet -> Sumw2();
  histeta_pytha_AND_bkgd_jet->SetXTitle("#eta^{jet} (pseudo-rap units)");
  histeta_pytha_AND_bkgd_jet->SetYTitle("dN/d#eta");

  TH1D *histphi_pytha_AND_bkgd_jet = new TH1D("histphi_pytha_AND_bkgd_jet", expression41,100,0.,2.*TMath::Pi());
  histphi_pytha_AND_bkgd_jet -> Sumw2();
  histphi_pytha_AND_bkgd_jet->SetXTitle("#phi^{jet} (radians)");
  histphi_pytha_AND_bkgd_jet->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histpT_area_corr_pytha_AND_bkgd_jet = new TH1D("histpT_area_corr_pytha_AND_bkgd_jet", expression42,300,-100,100);
  histpT_area_corr_pytha_AND_bkgd_jet -> Sumw2();
  histpT_area_corr_pytha_AND_bkgd_jet->SetXTitle("p_{T}^{are based correc. jet}");
  histpT_area_corr_pytha_AND_bkgd_jet->SetYTitle("dN/dp_{T}^{are based correc. jet}");

  TH1D *histnumtr_pytha_AND_bkgd_jet = new TH1D("histnumtr_pytha_AND_bkgd_jet", expression43,500,0,500);
  histnumtr_pytha_AND_bkgd_jet -> Sumw2();
  histnumtr_pytha_AND_bkgd_jet->SetXTitle("N_{tracks}");
  histnumtr_pytha_AND_bkgd_jet->SetYTitle("dN/dN_{tracks}");

  TH1D *histmeanpTtr_pytha_AND_bkgd_jet = new TH1D("histmeanpTtr_pytha_AND_bkgd_jet", expression44,200,0,100);
  histmeanpTtr_pytha_AND_bkgd_jet -> Sumw2();
  histmeanpTtr_pytha_AND_bkgd_jet->SetXTitle("<p_{T}^{track}> (GeV/c)");
  histmeanpTtr_pytha_AND_bkgd_jet->SetYTitle("dN/d<p_{T}^{track}>");

  TH1D *histangularity_pytha_AND_bkgd_jet = new TH1D("histangularity_pytha_AND_bkgd_jet", expression45,100,0.,10.);
  histangularity_pytha_AND_bkgd_jet -> Sumw2();
  histangularity_pytha_AND_bkgd_jet->SetXTitle("g (#phi - #eta distance units)");
  histangularity_pytha_AND_bkgd_jet->SetYTitle("dN/dg");

  TH1D *histldtr_pytha_AND_bkgd_jet = new TH1D("histldtr_pytha_AND_bkgd_jet", expression46,200,0.,100.);
  histldtr_pytha_AND_bkgd_jet -> Sumw2();
  histldtr_pytha_AND_bkgd_jet->SetXTitle("p_{T}^{ld. tr.} (GeV/c)");
  histldtr_pytha_AND_bkgd_jet->SetYTitle("dN/dp_{T}^{ld. tr.} ");

  TH1D *histsubldtr_pytha_AND_bkgd_jet = new TH1D("histsubldtr_pytha_AND_bkgd_jet", expression47,200,0.,100.);
  histsubldtr_pytha_AND_bkgd_jet -> Sumw2();
  histsubldtr_pytha_AND_bkgd_jet->SetXTitle("p_{T}^{sub-ld. tr.} (GeV/c)");
  histsubldtr_pytha_AND_bkgd_jet->SetYTitle("dN/dp_{T}^{sub-ld. tr.}");

  TH1D *histthrdldtr_pytha_AND_bkgd_jet = new TH1D("histthrdldtr_pytha_AND_bkgd_jet", expression48,200,0.,100.);
  histthrdldtr_pytha_AND_bkgd_jet -> Sumw2();
  histthrdldtr_pytha_AND_bkgd_jet->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histthrdldtr_pytha_AND_bkgd_jet->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histFF_pytha_AND_bkgd_jet = new TH1D("histFF_pytha_AND_bkgd_jet", expression49,100,0.,1.);
  histFF_pytha_AND_bkgd_jet -> Sumw2();
  histFF_pytha_AND_bkgd_jet->SetXTitle("z");
  histFF_pytha_AND_bkgd_jet->SetYTitle("dN/dz");

  TH1D *histlesub_pytha_AND_bkgd_jet = new TH1D("histlesub_pytha_AND_bkgd_jet", expression50,200,0.,100.);
  histlesub_pytha_AND_bkgd_jet -> Sumw2();
  histlesub_pytha_AND_bkgd_jet->SetXTitle("LeSub = p_{T}^{ld. tr.} - p_{T}^{sub-ld. tr.} (GeV/c)");
  histlesub_pytha_AND_bkgd_jet->SetYTitle("dN/dLeSub");

//______________________________________________________

  TH1D *histmedianrhoebe_pythia_AND_bkgd = new TH1D("histmedianrhoebe_pythia_AND_bkgd", expression51,600,0.,600.);
  histmedianrhoebe_pythia_AND_bkgd -> Sumw2();
  histmedianrhoebe_pythia_AND_bkgd->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histmedianrhoebe_pythia_AND_bkgd->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histX_pythia_tru_pythia_AND_bkgd = new TH1D("histX_pythia_tru_pythia_AND_bkgd", expression52,100,0.,1.);
  histX_pythia_tru_pythia_AND_bkgd -> Sumw2();
  histX_pythia_tru_pythia_AND_bkgd->SetXTitle("X^{tru}");
  histX_pythia_tru_pythia_AND_bkgd->SetYTitle("dN/dX^{tru}");

  TH1D *histX_bkgd_fake_pythia_AND_bkgd = new TH1D("histX_bkgd_fake_pythia_AND_bkgd", expression53,100,0.,1.);
  histX_bkgd_fake_pythia_AND_bkgd -> Sumw2();
  histX_bkgd_fake_pythia_AND_bkgd->SetXTitle("X^{fake}");
  histX_bkgd_fake_pythia_AND_bkgd->SetYTitle("dN/dX^{fake}");

  TH2D *histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw = new TH2D("histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw", expression52,100,0.,1.,100,0.,100.);
  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw -> Sumw2();
  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->SetXTitle("X^{tru}");
  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->SetYTitle("p_{T}^{jet, uncorr.}");
  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->SetZTitle("dN/dX^{tru}dp_{T}^{jet}");


//___________________________________________________________________________________________________//
//____________________________In-Plane_Angle_Cuts___________________________________//

  Double_t ip_low_1 = 0;
  Double_t ip_high_1 = 0.524;

  Double_t ip_low_2 = 5.760;
  Double_t ip_high_2 = 6.283;

  Double_t ip_low_3 = 2.618;
  Double_t ip_high_3 = 3.665;

//____________________________Mid-Plane_Angle_Cuts___________________________________//

  Double_t mp_low_1 = 0.524;
  Double_t mp_high_1 = 1.047;

  Double_t mp_low_2 = 5.236;
  Double_t mp_high_2 = 5.760;

  Double_t mp_low_3 = 2.094;
  Double_t mp_high_3 = 2.618;

  Double_t mp_low_4 = 3.665;
  Double_t mp_high_4 = 4.189;

//____________________________Out-Of-Plane_Angle_Cuts___________________________________//

  Double_t oop_low_1 = 1.047;
  Double_t oop_high_1 = 2.094;

  Double_t oop_low_2 = 4.189;
  Double_t oop_high_2 = 5.236;

//______________________________________________________________________________________//

  Int_t num_pythia_events = 1;
  Int_t num_background_events = 1;

//___________________________________________________________________________//

  //initializing run time stats

  Int_t* time_array;
  Int_t* time_interval_array;
  Int_t* event_array;

  time_array = new Int_t[nEvents + 1];
  time_interval_array = new Int_t[nEvents];
  event_array = new Int_t[nEvents];

  ///////////////////////////////////////////////////////////////////////
  //Creating a Jet Definition, choosing a Jet Algorithm
  ///////////////////////////////////////////////////////////////////////
  double Rparam = Jet_Radius; // <- Jet Radius  
  double ghost_maxrap = 2.0; //ghost particles go up to y=2
  int    repeat = 1; //only repeats once
  double ghost_area = 0.001; //area of ghost particles
  double grid_scatter = 1.0; //grid scatter 
  double pt_scatter = 0.1 ; //pt scatter
  double mean_ghost_pt = 1e-100; //mean pt of ghosts
  fastjet::Strategy               strategy = fastjet::Best;
  fastjet::RecombinationScheme    recombScheme = fastjet::WTA_pt_scheme;
  //fastjet::RecombinationScheme    recombScheme2 = fastjet::BIpt_scheme;
  fastjet::JetDefinition         *jetDef = NULL;
  fastjet::JetDefinition         *jetDef2 = NULL;
  fastjet::JetAlgorithm           algorithm=fastjet::antikt_algorithm;
  fastjet::JetAlgorithm           algorithm2=fastjet::kt_algorithm;
  fastjet::AreaType               area_type=fastjet::active_area_explicit_ghosts;
  fastjet::AreaDefinition        *AreaDef = NULL;
  fastjet::GhostedAreaSpec       *ghost_area_spec = NULL ;
  
  ghost_area_spec = new fastjet::GhostedAreaSpec(ghost_maxrap , repeat , ghost_area , grid_scatter, pt_scatter, mean_ghost_pt );
  AreaDef = new fastjet::AreaDefinition(area_type, GhostedAreaSpec(ghost_maxrap , repeat , ghost_area , grid_scatter, pt_scatter, mean_ghost_pt ));
  jetDef = new fastjet::JetDefinition(algorithm,Rparam,recombScheme,strategy);
  jetDef2 = new fastjet::JetDefinition(algorithm2,Rparam,recombScheme,strategy);
  std::vector <fastjet::PseudoJet> fjInputs_TOTAL;
  std::vector <fastjet::PseudoJet> fjInputs_TOTAL_kT;
  std::vector <fastjet::PseudoJet> fjInputs_Background;
  std::vector <fastjet::PseudoJet> fjInputs_Pythia_Particle;
  ///////////////////////////////////////////////////////////////////////

  //the speed of light in a vaccum is 2.998E8 m/s
  Double_t sol_vacuu = 299800000. ; // in m/s

  //initialize counting variables 
  Int_t N_counter = 0 ; //counts events
  TRandom3 *rando = new TRandom3((jobID)*17); //use this for the seed
  Int_t seed_b =  TMath::CeilNint((rando->Uniform(1.,1000.)));    

  // Create an instance of the Pythia event generator ...
  TPythia6* pythia = new TPythia6;
  //cout<<"I made it here = load pythia"<<endl;
  
  pythia -> SetMSTP(5, tune); //Tune A
  pythia -> SetCKIN(3,(jet_pT_cut_low - 1.)); // pt hard min (jet_pT_cut_low - 1.0 GeV/c)
  //pythia -> SetCKIN(3,49.);
   //pythia -> SetCKIN(4,40.); // pt hard max 40 GeV

  
  //The following block needs to be used on Newton Only
  //____________________________________________________________________
  TTimeStamp *timestamp = new TTimeStamp();
  UInt_t seed = ((jobID+1)*17) + (timestamp->GetSec() - static_cast<Int_t> (1446000000));
  if( (seed>=0) && (seed<=900000000) ) {
    //pythia->SetMRPY(1, ((rando->Uniform(1.,1000.))*(N_counter + 5 + seed)) );                   // set seed
    pythia->SetMRPY(1, seed);
    pythia->SetMRPY(2, 0);                      // use new seed
   // TRandom *Rand = new TRandom(seed); //this I will use for the pT efficiency
    cout<<"Random Seed : "<<seed<<endl;
  } else {cout << "error: time " << seed << " is not valid" << endl; exit(2);}
  
  //________Initiliaze Pythia and Background Generator________________________________________________//
  
  pythia -> Initialize("cms", "p", "p", 2760); //2760 for 2.76 TeV	  				    
  TClonesArray* particles = (TClonesArray*) pythia ->GetListOfParticles();
  Int_t N_reco_det_level = 0;
  TObjArray *particles_eff = new TObjArray(10000);
  Int_t N_reco_det_level_bkgd = 0;
  TObjArray *particles_eff_bkgd = new TObjArray(10000);

  TDatime *time = new TDatime(); //creating new time stamp
  cout<< "\n\nStart Time is " << time->GetTime() <<"\n\n\n" <<endl;
  time_array[0] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;

  Int_t user_index = 0;

  char run_title_str[256];
  sprintf( run_title_str, "BG_Gen_Cent%d_Harm%d" , cent_bin , HF);
 
  Float_t px; 
  Float_t py;
  Float_t pz;
  //Float_t m;
  //Float_t E;
  Short_t KF;
	  
  UInt_t centrality = cent_bin;
  UInt_t harmonics = HF;
  UInt_t version = 0;
  UInt_t events = nEvents;
                                                       
  BkgrLoad *bkgd2 = new BkgrLoad();

  UInt_t seed6 = ((jobID+834)*33) + (timestamp->GetSec() - static_cast<Int_t> (1446000000)); // for file shuffling

//____________________________________________________________________________________________________________________//

  char filepathstr[512];

  sprintf(filepathstr,"/home/alidock/ML/BKGD_ROOT_FILES");

  //getcwd(filepathstr, sizeof(filepathstr));

  bkgd2->PassInSettings(filepathstr,  0 , 0 , 0 , seed6 );
  bkgd2->Load();
  bkgd2->PrintStatus();
  TClonesArray *BKGD = (TClonesArray *)bkgd2->GetEventBG(); //(outside particle loop)

  
  std::vector <Double_t> rho_vec; //vector of rhos to get median from background only;
  std::vector <Double_t> rho_vec_TOTAL; //vector of rhos to get median from TOTAL while excluding 2 leading jets;

  //_____________________________________________Starting Event Loop_______________________________________________________________//

  for(Int_t event = 0; event < nEvents; event++) {
    //if( (bkgd2->GetActiveFile()) == (bkgd2->GetNumFiles()) ) //GetActiveFile() starts at 0 

    if( (bkgd2->GetActiveFile()) == 5 ){ //GetActiveFile() starts at 0  
      break; //break out of the loop if you have used up all the background events
    }
    if (event % 10 == 0) cout << "========================Event # " << event << endl;  
    N_counter ++;
    pythia->GenerateEvent();		
    Int_t npart = particles->GetEntries(); //This returns the number of particles in particles

    num_pythia_events++;
    histpythia_events->Fill(1.00);
    
      if( event != 0){
        free(time);
        TDatime *time = new TDatime(); //creating new time stamp
        //cout << "\n\n\n\n\nTime at Event # "<< event << " is " << time->GetTime() <<"\n\n\n\n" << endl;
        time_array[event] = ((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        Int_t user_index = 0; //must reset this each event iteration loop
        Int_t N_reco_det_level = 0; //must reset this each event iteration loop
        Int_t N_reco_det_level_bkgd = 0;
      }//end get new time stamp for dianostic purposes

//___________________________________Starting Pythia Particle Loop___________________________________________//
    for (Int_t part=0; part<npart; part++) {
    
      TMCParticle *MPart = (TMCParticle *) particles->At(part); //indexing the particles array for the current member in loop
     
      Double_t pt =  TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy() );
      Double_t p = TMath::Sqrt(MPart->GetPx() * MPart->GetPx() + MPart->GetPy() * MPart->GetPy() + MPart->GetPz() * MPart->GetPz() );
      Double_t E = MPart -> GetEnergy();
      Float_t Vx = MPart -> GetVx(); //x coordinate of particle vertex in mm
      Float_t Vy = MPart -> GetVy(); //y coordinate of particle vertex in mm 
      Float_t Vz = MPart -> GetVz(); //z coordinate of particle vertex in mm
      Int_t KF = MPart -> GetKF(); // getting flavour code
      Float_t DCA_calc = (TMath::Sqrt( (Vx*Vx) + (Vy*Vy) + (Vz*Vz) )) * (1000000/1000) ; //answer in microns
      Double_t p_eta =  0.5*(TMath::Log((p+MPart->GetPz())/(p-MPart->GetPz()))) ;
      Double_t phi = 1.0*(TMath::Pi()+TMath::ATan2(-MPart->GetPy(),-MPart->GetPx()));
     
      if(pt > constit_cut && pt < jet_pT_cut_high && (TMath::Abs(p_eta)) < 0.9  ){ //looking for particles in the EMCAL with pT = 0.15 GeV/c (c = 1)   
        if( DCA_calc < DCA ){ //you want primary particles only

            fastjet::PseudoJet p1( MPart->GetPx(), MPart->GetPy() , MPart->GetPz() , MPart->GetEnergy() );//creating temp pseudojet p1
            p1.set_user_index(user_index);
            fjInputs_Pythia_Particle.push_back(p1); //store in pseudo-jets to be clustered for PYTHIA only
            fjInputs_TOTAL.push_back(p1); //store in pseudo-jets to be clustered for PYTHIA + bkgd 
            fjInputs_TOTAL_kT.push_back(p1); //meant for a realistic assessment of mixed signal + bkgd

            user_index++;

            histpT_pyth_part->Fill(pt);
            histeta_pyth_part->Fill(p_eta);
            histphi_pyth_part->Fill(phi);

            histpT_pytha_AND_bkgd_part->Fill(pt);
            histeta_pytha_AND_bkgd_part->Fill(p_eta);
            histphi_pytha_AND_bkgd_part->Fill(phi);

        } //end if DCA cut
      } //end if pT soft cut and eta acceptance cut


    }//end Pythia Particle Loop


   //apply jet definitions 
    vector <fastjet::PseudoJet> inclusiveJetsPythia_part, sortedJetsPythia_part; // 
    fastjet::ClusterSequenceArea clustSeqPythia_Particle(fjInputs_Pythia_Particle, *jetDef , *AreaDef ); // 
    //fastjet::ClusterSequence clustSeqPythia_Particle(fjInputs_Pythia_Particle, *jetDef ); // 

   //fjInputs_Pythia_Particle.clear();

    //MCTRUTH (all particles)
    ////////////////////////////////////////////
    inclusiveJetsPythia_part = clustSeqPythia_Particle.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
    ////////////////////////////////////////////  

   //sorting
    ////////////////////////////////////////////
    sortedJetsPythia_part   = sorted_by_pt(inclusiveJetsPythia_part); //sort by decreasing transverse momentum
    ////////////////////////////////////////////


    /////////////////////////making jet cuts on eta and pT 
    fastjet::Selector select_pt = fastjet::SelectorPtRange(jet_pT_cut_low,jet_pT_cut_high); // Selects Jets with transverse momentum between jet_pT_cut_low and jet_pT_cut_high 
    /////////////////////////making jet cuts on eta and pT 
    fastjet::Selector select_pt2 = fastjet::SelectorPtRange(constit_cut,1000); // Selects Jets with transverse momentum between constit_cut and 1000 GeV
    /////////////////////////making jet cuts on eta and pT 
    fastjet::Selector select_pt3 = fastjet::SelectorPtRange(jet_pT_cut_low,100); // Selects Jets with transverse momentum between jet_pT_cut_low and 100 GeV
    fastjet::Selector select_rapidity = fastjet::SelectorEtaRange(-0.9 + Rparam , 0.9 - Rparam); // Selects Jets in the desired eta range
    fastjet::Selector select_both = select_pt && select_rapidity ; //combining the two jet cuts
    fastjet::Selector select_both2 = select_pt2 && select_rapidity ; //combining the two jet cuts for the second pT selector
    fastjet::Selector select_both3 = select_pt3 && select_rapidity ; //combining the two jet cuts for the third pT selector
    ////////////////////////////////////////

    vector<fastjet::PseudoJet>selected_jetsPythia_part = select_both(sortedJetsPythia_part); 
    vector<fastjet::PseudoJet>selected_jetsPythia_sorted_part = sorted_by_pt(selected_jetsPythia_part);
    Int_t selected_jetsPythia_sorted_part_size = selected_jetsPythia_sorted_part.size();

/////////////______________FILLING THE PYTHIA JET HISTOGRAMS_________________________//////////////////////

    if( selected_jetsPythia_sorted_part_size > 0){ 

      for(Int_t py_jet_ind2 = 0; py_jet_ind2 < selected_jetsPythia_sorted_part_size ; py_jet_ind2++ ){
        histpythia_jets->Fill(1.00);
        
        histpT_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].pt() ); //basic jet properties
        histeta_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].eta() );
        histphi_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].phi() );

        Double_t jet_phi = selected_jetsPythia_sorted_part[py_jet_ind2].phi();
        Double_t jet_eta = selected_jetsPythia_sorted_part[py_jet_ind2].eta();

        Double_t jet_pt = selected_jetsPythia_sorted_part[py_jet_ind2].pt();

        vector <fastjet::PseudoJet> constituents_part = selected_jetsPythia_sorted_part[py_jet_ind2].constituents(); //grab the constituents
        vector <Double_t> pythia_constit_pT_vec;

        Int_t nPart_part = constituents_part.size();

        Double_t angularity_sum =0;
        Int_t num_const_gt_cut =0;

        
        for (Int_t np_part = 0; np_part < nPart_part; np_part++){ //loop over the constituents
          if(constituents_part[np_part].perp() > 1e-50){ //MAKE SURE you have no ghosts
            num_const_gt_cut++;
            pythia_constit_pT_vec.push_back(constituents_part[np_part].perp()); 

            Double_t z = constituents_part[np_part].perp()/jet_pt;

            Double_t pphi = constituents_part[np_part].phi();
            Double_t peta = constituents_part[np_part].eta();

            Double_t delri = TMath::Sqrt( ((TMath::Abs(jet_phi - pphi))*(TMath::Abs(jet_phi - pphi))) + ((TMath::Abs(jet_eta - peta))*(TMath::Abs(jet_eta - peta))) );
   
            histFF_pyth_jet->Fill(z); //frag func

            angularity_sum += (z*delri); //summing for the angularity
          }
	}

        histangularity_pyth_jet->Fill(angularity_sum); //angularity
        histnumtr_pyth_jet->Fill(num_const_gt_cut); //number of consituents
        constituents_part.clear();

        if(pythia_constit_pT_vec.size() > 0){
          
         // histmeanpTtr_pyth_jet->Fill( std::accumulate(pythia_constit_pT_vec.begin(), pythia_constit_pT_vec.end(), 0.0) / pythia_constit_pT_vec.size() ); //mean track pT in the jet

          std::sort(pythia_constit_pT_vec.begin(), pythia_constit_pT_vec.end()); //sorts least to greatest
          histldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size()-1]); // grab leading track in jet

          if(pythia_constit_pT_vec.size() > 1){
            histsubldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 2] ); // grab sub-leading track in jet
            histlesub_pyth_jet->Fill( pythia_constit_pT_vec[pythia_constit_pT_vec.size()-1] - pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 2] );
            if(pythia_constit_pT_vec.size() > 2){
              histthrdldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 3]); //grad sub-sub leading track in jet
            } //for sub sub leading
          } //for sub leading
        } //make sure there are constits


      pythia_constit_pT_vec.clear();
      } //end looping over pythia jets
    } //end if there are jets that pass your cuts

//_______________________________________________________________LETS GET INTO THE BACKGROUND NOW___________________________________//

    if( selected_jetsPythia_sorted_part.size()==0){
      cout<<"\n\nNO jets in Pythia jets group AFTER CUTS  !\n\n"<<endl;
    }

    else if( selected_jetsPythia_sorted_part.size() > 0) { //you had a PYTHIA jet that passed your cuts

      if( event != 0){
        delete BKGD;
        TClonesArray *BKGD = (TClonesArray *)bkgd2->GetEventBG(); //(outside particle loop)
        num_background_events++;
        cout<<"you got the new background event\n\n"<<endl;
      } 

      histbackground_events->Fill(1.00);

      Double_t mass;
      Double_t nParticles= BKGD ->GetEntries();

      Double_t Sum_Psi_1_numer = 0;
      Double_t Sum_Psi_1_denom = 0;

      Double_t Sum_Psi_2_numer = 0;
      Double_t Sum_Psi_2_denom = 0;

      Double_t Sum_Psi_3_numer = 0;
      Double_t Sum_Psi_3_denom = 0;

      Double_t Sum_Psi_4_numer = 0;
      Double_t Sum_Psi_4_denom = 0;

      Double_t Sum_Psi_5_numer = 0;
      Double_t Sum_Psi_5_denom = 0;

    //____________________________Starting Background Particle Loop___________________________________________________//
      for(Int_t N = 0 ; N < nParticles ; N++) {
        TMCParticle* bob = (TMCParticle*)BKGD->At(N); 
        Int_t K_F = bob -> GetKF();

        if (K_F == 211 || K_F == -211) mass = 0.13957; //pi+ & pi- (GeV/c^2) (inside particle loop)
	else if (K_F  == 321 || K_F == -321) mass = 0.493677; //k+ && k- (GeV/c^2)
	else if (K_F == 2212 || K_F == -2212) mass = 0.938272; //p & pbar (GeV/c^2)
	else if (K_F == 111) mass = 0.134977; //pi0 (GeV/c^2)

        Double_t p = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) + ((bob->GetPz())*(bob->GetPz())) );
        Double_t pT = TMath::Sqrt( ((bob->GetPx())*(bob->GetPx())) + ((bob->GetPy())*(bob->GetPy())) );
        Double_t eta = 0.5*(TMath::Log((p+bob->GetPz())/(p-bob->GetPz())));
        Double_t Energy = TMath::Sqrt( (p*p) + (mass*mass) );
        Double_t phi = 1.0*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx()));
      
        fastjet::PseudoJet pa( bob->GetPx(), bob->GetPy() , bob->GetPz() , Energy );//creating temp pseudojet pa
        pa.set_user_index(-1); //we want all the *Fake* particles to have negative 1 as 

        if( pT > constit_cut ){

          histpT_bkgd_part->Fill(pT);
          histpT_pytha_AND_bkgd_part->Fill(pT);

          histeta_bkgd_part->Fill(eta);
          histeta_pytha_AND_bkgd_part->Fill(eta);

          histphi_bkgd_part->Fill( phi ); //azimuthal angle
          histphi_pytha_AND_bkgd_part->Fill( phi ); //azimuthal angle

          //FOR FILLING THE EVENT PLANE ANGLES (USING IT AS A CHECK, EVEN ANGLES SHOULD MOSTLY BE 0, ODD ANGLES RANDOM)
          Sum_Psi_1_numer = Sum_Psi_1_numer + pT*TMath::Sin( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_1_denom = Sum_Psi_1_denom + pT*TMath::Cos( (1.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_2_numer = Sum_Psi_2_numer + pT*TMath::Sin( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_2_denom = Sum_Psi_2_denom + pT*TMath::Cos( (2.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_3_numer = Sum_Psi_3_numer + pT*TMath::Sin( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_3_denom = Sum_Psi_3_denom + pT*TMath::Cos( (3.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_4_numer = Sum_Psi_4_numer + pT*TMath::Sin( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_4_denom = Sum_Psi_4_denom + pT*TMath::Cos( (4.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          Sum_Psi_5_numer = Sum_Psi_5_numer + pT*TMath::Sin( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;
          Sum_Psi_5_denom = Sum_Psi_5_denom + pT*TMath::Cos( (5.0)*(TMath::Pi()+TMath::ATan2(-bob->GetPy(),-bob->GetPx())) ) ;

          fjInputs_Background.push_back(pa); //store in pseudo-jets to be clustered for bkgd only
          fjInputs_TOTAL.push_back(pa); //store in pseudo-jets to be clustered for PYTHIA + bkgd 
          fjInputs_TOTAL_kT.push_back(pa); //meant for a realistic assessment of mixed signal + bkgd
         }
     
      }//END PARTICLE LOOP BACKGROUND

      histPsi_1->Fill( (1.0/1.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_1_numer,-Sum_Psi_1_denom)));
      histPsi_2->Fill( (1.0/2.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_2_numer,-Sum_Psi_2_denom)));
      histPsi_3->Fill( (1.0/3.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_3_numer,-Sum_Psi_3_denom)));
      histPsi_4->Fill( (1.0/4.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_4_numer,-Sum_Psi_4_denom)));
      histPsi_5->Fill( (1.0/5.0)*(TMath::Pi()+TMath::ATan2(-Sum_Psi_5_numer,-Sum_Psi_5_denom)));

      if (fjInputs_TOTAL.size() == 0  ) {
        cout << "Error: event with no final state particles" << endl;
        continue;
      } //endif looking for events with final state particles 

      //apply jet definitions

      vector <fastjet::PseudoJet> inclusiveJetsTOTAL, sortedJetsTOTAL ; // 
      fastjet::ClusterSequenceArea clustSeqTOTAL(fjInputs_TOTAL, *jetDef , *AreaDef ); 

      vector <fastjet::PseudoJet> inclusiveJetsBackground, sortedJetsBackground; // 
      fastjet::ClusterSequenceArea clustSeqBackground(fjInputs_Background, *jetDef2 , *AreaDef ); //
      //fastjet::ClusterSequence clustSeqBackground(fjInputs_Background, *jetDef2 ); // 
      fjInputs_Background.clear();

      vector <fastjet::PseudoJet> inclusiveJetsTOTAL_kT, sortedJets_TOTAL_kT; // 
      fastjet::ClusterSequenceArea clustSeqTOTAL_kT(fjInputs_TOTAL_kT, *jetDef2 , *AreaDef ); // 

      ////////////////////////////////////////////
      inclusiveJetsTOTAL = clustSeqTOTAL.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      ////////////////////////////////////////////

      ////////////////////////////////////////////
      inclusiveJetsBackground = clustSeqBackground.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      inclusiveJetsTOTAL_kT = clustSeqTOTAL_kT.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      ////////////////////////////////////////////


      //sorting
      ////////////////////////////////////////////
      sortedJetsTOTAL   = sorted_by_pt(inclusiveJetsTOTAL); //sort by decreasing transverse momentum
      sortedJetsBackground   = sorted_by_pt(inclusiveJetsBackground); //sort by decreasing transverse momentum
      sortedJets_TOTAL_kT = sorted_by_pt(inclusiveJetsTOTAL_kT); //sort by decreasing transverse momentum
      ////////////////////////////////////////////

      vector<fastjet::PseudoJet>selected_jetsTOTAL = select_both3(sortedJetsTOTAL); 
      vector<fastjet::PseudoJet>selected_jetsTOTAL_sorted = sorted_by_pt(selected_jetsTOTAL);

      vector<fastjet::PseudoJet>selected_jetsBackground = select_both2(sortedJetsBackground); 
      vector<fastjet::PseudoJet>selected_jetsBackground_sorted = sorted_by_pt(selected_jetsBackground);

      vector<fastjet::PseudoJet>selected_jetsTOTAL_kT = select_both2(sortedJets_TOTAL_kT); 
      vector<fastjet::PseudoJet>selected_jetsTOTAL_kT_sorted = sorted_by_pt(selected_jetsTOTAL_kT);


      Double_t running_rho_average = 0;

      cout<<"\n\n\n|||||||||||The total background jet population is = "<<selected_jetsBackground_sorted.size()<<"|||||||||||||||||"<<endl;

      if(selected_jetsBackground_sorted.size() > 0){
        for(unsigned i_jet = 0 ; i_jet < selected_jetsBackground_sorted.size() ; i_jet++ ){ 

          histbackground_jets->Fill(1.00);
    
          histpT_bkgd_jet -> Fill(selected_jetsBackground_sorted[i_jet].pt());
          histeta_bkgd_jet -> Fill(selected_jetsBackground_sorted[i_jet].eta());
          histphi_bkgd_jet -> Fill(selected_jetsBackground_sorted[i_jet].phi());

          Double_t jet_pT_back = selected_jetsBackground_sorted[i_jet].pt();
          Double_t jet_eta_back = selected_jetsBackground_sorted[i_jet].eta();
          Double_t jet_phi_back = selected_jetsBackground_sorted[i_jet].phi();

          rho_vec.push_back(selected_jetsBackground_sorted[i_jet].pt()/selected_jetsBackground_sorted[i_jet].area());
      
          vector <fastjet::PseudoJet> constituents_bkgd = selected_jetsBackground_sorted[i_jet].constituents();
          vector <Double_t> bkgd_constit_pT_vec;

          Int_t nPart_bkgd = constituents_bkgd.size();
          histnumtr_bkgd_jet->Fill(constituents_bkgd.size()); //number of consituents

          Double_t angularity_sum_back =0;
        
          for(Int_t np_bkgd = 0; np_bkgd < nPart_bkgd; np_bkgd++){
            if(constituents_bkgd[np_bkgd].perp() > 1e-50){ //don't use ghosts
              bkgd_constit_pT_vec.push_back(constituents_bkgd[np_bkgd].perp());    

              Double_t z = constituents_bkgd[np_bkgd].perp()/jet_pT_back;
              histFF_bkgd_jet->Fill(z);

              Double_t pphi_back = constituents_bkgd[np_bkgd].phi();
              Double_t peta_back = constituents_bkgd[np_bkgd].eta();

              Double_t delri_back = TMath::Sqrt( ((TMath::Abs(jet_phi_back - pphi_back))*(TMath::Abs(jet_phi_back - pphi_back))) + ((TMath::Abs(jet_eta_back - peta_back))*(TMath::Abs(jet_eta_back - peta_back))) );

              angularity_sum_back += (delri_back * z);
            }
          }

          constituents_bkgd.clear(); 
          histangularity_bkgd_jet->Fill(angularity_sum_back);

          if(bkgd_constit_pT_vec.size() > 0){
          
           // histmeanpTtr_bkgd_jet->Fill( std::accumulate(bkgd_constit_pT_vec.begin(), bkgd_constit_pT_vec.end(), 0.0) / bkgd_constit_pT_vec.size() ); //mean track pT in the jet

            std::sort(bkgd_constit_pT_vec.begin(), bkgd_constit_pT_vec.end()); //sorts least to greatest
            histldtr_bkgd_jet->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size()-1]); // grab leading track in jet

            if(bkgd_constit_pT_vec.size() > 1){
              histsubldtr_bkgd_jet->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 2]); // grab sub-leading track in jet
              histlesub_bkgd_jet->Fill( bkgd_constit_pT_vec[bkgd_constit_pT_vec.size()-1] - bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 2] );
              if(bkgd_constit_pT_vec.size() > 2){
                histthrdldtr_bkgd_jet->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 3]); //grad sub-sub leading track in jet
              }
            }


          } //end if check for constits
  
          bkgd_constit_pT_vec.clear();
      
        } //end looping over the background jets

        Double_t event_median;
        Int_t index_rho_vec;
        
        if(rho_vec.size() > 0){
          std::sort(rho_vec.begin(), rho_vec.end());
          index_rho_vec = TMath::Ceil(rho_vec.size() / 2);
          if( rho_vec.size() % 2 == 0){ //median if even number of detector level jets
            event_median = 0.5*(rho_vec[(rho_vec.size()/2)-1] + rho_vec[(rho_vec.size()/2)])*TMath::Pi()*Jet_Radius*Jet_Radius;
          }
          else {
            event_median = rho_vec[index_rho_vec]*TMath::Pi()*Jet_Radius*Jet_Radius; //median if odd number of detector level jets
          }
        }
        else {
          event_median = 0.;
        }
        
        cout<<"\n\n\n|||||||||||||The median rho*area is = "<<event_median<<"||||||||||||||||||"<<endl;
        histmedianrhoebe_bgkd->Fill(event_median/(TMath::Pi()*Jet_Radius*Jet_Radius));
        rho_vec.clear();

      } //end if check for background jets that pass your cuts

       
      Double_t event_median_realistic;

      if(sortedJets_TOTAL_kT.size()>0){

//__________________________MOVE INTO BACKGROUND + PYTHIA JUST TO GET THE MEDIAN RHO (WITH REALISITIC LEADING JET SUPPRESSION)_________//

        if(selected_jetsTOTAL_kT_sorted.size() > 0){ //you must have at least 3 kT jets in TOTAL
          for(unsigned i = 0 ; i < selected_jetsTOTAL_kT_sorted.size() ; i++ ){ //exclude no leading jet !!!!!!!
            rho_vec_TOTAL.push_back(selected_jetsTOTAL_kT_sorted[i].pt()/selected_jetsTOTAL_kT_sorted[i].area());
          }
        }
        else if(selected_jetsTOTAL_kT_sorted.size() > 1){
          rho_vec_TOTAL.clear();
          for(unsigned i = 1 ; i < selected_jetsTOTAL_kT_sorted.size() ; i++ ){ //exclude the leading jet !!!!!!!
            rho_vec_TOTAL.push_back(selected_jetsTOTAL_kT_sorted[i].pt()/selected_jetsTOTAL_kT_sorted[i].area());
          }
        }
        else if(selected_jetsTOTAL_kT_sorted.size() > 2){ //most restricitive
          rho_vec_TOTAL.clear();
          for(unsigned i = 2 ; i < selected_jetsTOTAL_kT_sorted.size() ; i++ ){ //exclude the two leading jets !!!!!!!
            rho_vec_TOTAL.push_back(selected_jetsTOTAL_kT_sorted[i].pt()/selected_jetsTOTAL_kT_sorted[i].area());
          }
        }

        Int_t index_rho_vec_realistic;

        if(rho_vec_TOTAL.size() > 0){
          std::sort(rho_vec_TOTAL.begin(), rho_vec_TOTAL.end());
          index_rho_vec_realistic = TMath::Ceil(rho_vec_TOTAL.size() / 2);
          if( rho_vec_TOTAL.size() % 2 == 0){ //median if even number of detector level jets
            event_median_realistic = 0.5*(rho_vec_TOTAL[(rho_vec_TOTAL.size()/2)-1] + rho_vec_TOTAL[(rho_vec_TOTAL.size()/2)])*TMath::Pi()*Jet_Radius*Jet_Radius;
          }
          else {
            event_median_realistic = rho_vec_TOTAL[index_rho_vec_realistic]*TMath::Pi()*Jet_Radius*Jet_Radius; //median if odd number of detector level jets
          }
        }
        else {
          event_median_realistic = 0.;
        }
      }
      cout<<"\n\n\n|||||||||||||The median rho*area (drop 2) is = "<<event_median_realistic<<"||||||||||||||||||"<<endl;
      histmedianrhoebe_pythia_AND_bkgd->Fill(event_median_realistic/(TMath::Pi()*Jet_Radius*Jet_Radius));
      rho_vec_TOTAL.clear();

//__________________________NOW GET INTO THE BACKGROUND + PYTHIA _______________________________________________________________//

      if( selected_jetsTOTAL_sorted.size()==0){
        cout<<"\n\nNO jets in re-clustered Total (background + pythia) jets group !\n\n"<<endl;
        rho_vec.clear();
      }
      else if(selected_jetsTOTAL_sorted.size()>0){
        ofstream pbgOut;
        pbgOut.open("pbgOut.csv", fstream::app);
        if(HEADER__){
          pbgOut<<"p_T, Eta, Phi, p_T-corr, N-Trk, Angularity, Mean-p_T, p_T_1, p_T_2, p_T_3, p_T_4, p_T_5, X_tru"<<endl;
          HEADER__--;
        }

        for(unsigned t_jet = 0 ; t_jet < selected_jetsTOTAL_sorted.size() ; t_jet++ ){ 

          histpythia_AND_bkgd_jets->Fill(1.00);
    
          histpT_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].pt());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].pt()<<", ";

          histeta_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].eta());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].eta()<<", ";

          histphi_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].phi());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].phi()<<", ";

          histpT_area_corr_pytha_AND_bkgd_jet->Fill(selected_jetsTOTAL_sorted[t_jet].pt() - event_median_realistic );
          pbgOut<<selected_jetsTOTAL_sorted[t_jet].pt()-event_median_realistic<<", ";
      
          vector <fastjet::PseudoJet> constituents_total = selected_jetsTOTAL_sorted[t_jet].constituents();
          vector <Double_t> total_constit_pT_vec;

          Int_t nPart_total = constituents_total.size();
          histnumtr_pytha_AND_bkgd_jet->Fill(constituents_total.size()); //number of consituents
          pbgOut<<constituents_total.size()<<", ";

          Double_t angularity_sum =0;

          Double_t tru_sum=0;
          Double_t fake_sum=0;
        
          for(Int_t np_total = 0; np_total < nPart_total; np_total++){
            if(constituents_total[np_total].perp() > 1e-50){ //don't grab the ghosts !!!!

              total_constit_pT_vec.push_back(constituents_total[np_total].perp());    
              histFF_pytha_AND_bkgd_jet->Fill(constituents_total[np_total].perp()/selected_jetsTOTAL_sorted[t_jet].pt() ); //frag func

              

              angularity_sum += ( (constituents_total[np_total].perp()/selected_jetsTOTAL_sorted[t_jet].pt())*TMath::Sqrt( (selected_jetsTOTAL_sorted[t_jet].phi() - constituents_total[np_total].phi())*(selected_jetsTOTAL_sorted[t_jet].phi() - constituents_total[np_total].phi()) + (selected_jetsTOTAL_sorted[t_jet].eta() - constituents_total[np_total].eta())*(selected_jetsTOTAL_sorted[t_jet].eta() - constituents_total[np_total].eta()) )); //summing for the angularity

              if(constituents_total[np_total].user_index() < 0){
                fake_sum = fake_sum+constituents_total[np_total].perp();
              }
              else{
                tru_sum = tru_sum+constituents_total[np_total].perp();
              }
            }
          }

          Double_t X_tru = tru_sum/total_constit_pT_vec.size();
          Double_t X_fake = fake_sum/total_constit_pT_vec.size();

          //cout<<"\ntrusum = "<<tru_sum<<" jet size (no ghosts) = "<<total_constit_pT_vec.size()<<" Xtru = "<<X_tru<<endl;
          //cout<<"fakesum = "<<fake_sum<<" jet size (no ghosts) = "<<total_constit_pT_vec.size()<<" Xfake = "<<X_fake<<endl;

          histX_pythia_tru_pythia_AND_bkgd->Fill(X_tru);

          Int_t label = X_tru<0.1?0:X_tru<0.5?1:X_tru<0.9?2:3;//0-fake, 1-mostly fake, 2-mostly real, 3-real

          

          histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->Fill(X_tru, selected_jetsTOTAL_sorted[t_jet].pt());

          histX_bkgd_fake_pythia_AND_bkgd->Fill(X_fake);

          histangularity_pytha_AND_bkgd_jet->Fill(angularity_sum);

          pbgOut<<angularity_sum<<", ";

          if(total_constit_pT_vec.size() > 0){
          
            histmeanpTtr_pytha_AND_bkgd_jet->Fill( std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size() ); //mean track pT in the jet

            pbgOut<<std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size()<<", ";

            std::sort(total_constit_pT_vec.begin(), total_constit_pT_vec.end()); //sorts least to greatest
            histldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 1]); // grab leading track in jet

            pbgOut<<total_constit_pT_vec[total_constit_pT_vec.size() - 1]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 2]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 3]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 4]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 5]<<", ";


            if(total_constit_pT_vec.size() > 1){
              histsubldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 2]); // grab sub-leading track in jet
              histlesub_pytha_AND_bkgd_jet->Fill( total_constit_pT_vec[total_constit_pT_vec.size() - 1] - total_constit_pT_vec[total_constit_pT_vec.size() - 2] );
              if(total_constit_pT_vec.size() > 2){
                histthrdldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 3]); //grab sub-sub leading track in jet
              }
            }
          } //end if check for total jet constit
          pbgOut<<X_tru<<endl;
           total_constit_pT_vec.clear();
          constituents_total.clear();

        } //end looping over the background + pythiajets
        pbgOut.close();
      } //end if you have jets in the total




//________________________________________________________________________________________________________________________________//

    } //endif you have a pythia jet that passes your cuts

    fjInputs_TOTAL.clear();
    fjInputs_TOTAL_kT.clear();
    fjInputs_Pythia_Particle.clear();

    rho_vec.clear();
    rho_vec_TOTAL.clear();
  }//________________________________________________END EVENT LOOP (line 1749)________________________________________

  bkgd2->CloseFiles(); //close all the damn files
  //bkgd->~BkgrLoad();

  TFile *ff = new TFile(expression1,"RECREATE"); //make output files and directories
  TDirectory *pythia_only = ff->mkdir("pythia_only");
  TDirectory *background_only = ff->mkdir("background_only");
  TDirectory *background_and_pythia = ff->mkdir("background_and_pythia");

  Int_t Sum_Time = 0;
  Int_t Total_Time = 0;

  for( Int_t i = 0 ; i < nEvents; i++){
    time_interval_array[i] = time_array[i + 1] - time_array[i];
    Sum_Time = Sum_Time + time_interval_array[i];
    event_array[i] = i + 1;
  }

  if(RT_Stats){
    TCanvas *c1 = new TCanvas("c1" , "Run Time For 1 Background Event vs Event Number (Events Run in Sequence)", 200, 10,500,300); 

    TH1D *histruntime = new TH1D("histruntime", "Run Times for Background Generator",8,((Sum_Time/nEvents) - 1),((Sum_Time/nEvents) + 1)); //for non-random distribution
    histruntime -> Sumw2();
    histruntime->SetXTitle("Run Time For 1 Event (Seconds)");
    histruntime->SetYTitle("Counts");

    TGraph *rtgraph = new TGraph(nEvents , event_array , time_interval_array);
    rtgraph->SetTitle("Run Time For 1 Background Event vs Event Number (Events Run in Sequence)");
    rtgraph->GetXaxis()->SetTitle("Event Iteration Number");
    rtgraph->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
    rtgraph->Write();
    c1->cd(1); 
    rtgraph->Draw("AC");
    c1->Draw();
    c1->Write("");
    for( Int_t j = 0 ; j < nEvents; j++){
      histruntime->Fill(time_interval_array[j]);
    }
    histruntime->Write();
  }


//_______________________________________PYTHIA ONLY STUFF FIRST______________________________________________________

  pythia_only->cd(); 

  histpythia_events->Write();
  histpythia_jets->Write();

  histpT_pyth_part->Write();
  histeta_pyth_part->Write();
  histphi_pyth_part->Write();

  histpT_pyth_jet->Write();
  histeta_pyth_jet->Write();
  histphi_pyth_jet->Write();

  histnumtr_pyth_jet->Write();
  histmeanpTtr_pyth_jet->Write();
  histangularity_pyth_jet->Write();

  histldtr_pyth_jet->Write();
  histsubldtr_pyth_jet->Write();

  histthrdldtr_pyth_jet->Write();
  histFF_pyth_jet->Write();
  histlesub_pyth_jet->Write();
  

//______________________________________END PYTHIA ONLY STUFF______________________________________________________________//

//______________________________________START OF BACKGROUND ONLY STUFF______________________________________________________________//

  background_only->cd();

  histbackground_events->Write();
  histbackground_jets->Write();

  histpT_bkgd_part->Write();
  histeta_bkgd_part->Write();
  histphi_bkgd_part->Write();

  histpT_bkgd_jet->Write();
  histeta_bkgd_jet->Write();
  histphi_bkgd_jet->Write();

  histnumtr_bkgd_jet->Write();
  histmeanpTtr_bkgd_jet->Write();
  histangularity_bkgd_jet->Write();

  histldtr_bkgd_jet->Write();
  histsubldtr_bkgd_jet->Write();
  histthrdldtr_bkgd_jet->Write();

  histFF_bkgd_jet->Write();
  histlesub_bkgd_jet->Write();

  histPsi_1->Write();
  histPsi_2->Write();
  histPsi_3->Write();
  histPsi_4->Write();
  histPsi_5->Write();

  histmedianrhoebe_bgkd->Write();


//______________________________________END BACKGROUND ONLY STUFF______________________________________________________________// 

//_______________________________________START OF BACKGROUND AND PYTHIA_________________________________________________________//

  background_and_pythia->cd();

  histpythia_AND_bkgd_jets->Write();

  histpT_pytha_AND_bkgd_part->Write();
  histeta_pytha_AND_bkgd_part->Write();
  histphi_pytha_AND_bkgd_part->Write();

  histpT_pytha_AND_bkgd_jet->Write();
  histeta_pytha_AND_bkgd_jet->Write();
  histphi_pytha_AND_bkgd_jet->Write();

  histpT_area_corr_pytha_AND_bkgd_jet->Write();
  histnumtr_pytha_AND_bkgd_jet->Write();
  histmeanpTtr_pytha_AND_bkgd_jet->Write();

  histangularity_pytha_AND_bkgd_jet->Write();

  histldtr_pytha_AND_bkgd_jet->Write();
  histsubldtr_pytha_AND_bkgd_jet->Write();
  histthrdldtr_pytha_AND_bkgd_jet->Write();

  histFF_pytha_AND_bkgd_jet->Write();
  histlesub_pytha_AND_bkgd_jet->Write();

  histmedianrhoebe_pythia_AND_bkgd->Write();

  histX_pythia_tru_pythia_AND_bkgd->Write();
  histX_bkgd_fake_pythia_AND_bkgd->Write();

  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->Write();



//_______________________________________END OF BACKGROUND AND PYTHIA_________________________________________________________//

  ff->Write();
  ff->Close();

} //end function

