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
void Toy_Model_ML_Study(Int_t nEvents, Int_t jobID , Int_t tune, Double_t Jet_Radius, Int_t HF , Double_t DCA , Int_t cent_bin , Int_t corrjet_bin , Double_t constit_cut, Bool_t Data = kTRUE, Bool_t RT_Stats = kFALSE , Bool_t GRID = kFALSE ){

  //Int_t HF = 0; // set harmonics flag (0 : v1 - v5) , (1 : v2 - v5) , (2: v3 - v5) , (3: v1 - v4) , (4: v1 - v3) , (5: no vn, uniform phi) , (6: v1 - v2 , v4 - v5) , (7: v1 -v3 , v5) , (8: v1 , v3 - v5) , (9: v1 only) , (10: v2 only) , (11: v3 only) , (12: v4 only) , (13: v5 only)

 // cent_bin: (0 : 0 - 5 %) , (1 : 5 - 10 % ) , (2 : 10 - 20 % ) , (3 : 20 - 30 % ) , (4 : 30 - 40 %) , (5 : 40 - 50 %) , (6 : 50 - 60 %) , (7: 60 - 70 %)

// corrjet_bin: (1 : 10 - 30 GeV , pT_hard_min_pythia = 9 GeV) , (2 : 30 - 50 GeV , pT_hard_min_pythia = 29 GeV) , (3: 50 - 70 GeV , pT_hard_min_pythia = 49 GeV)

  auto jet_constit_loop = []( fastjet::PseudoJet parent_jet ,  std::vector<fastjet::PseudoJet> all_constits , TH2D *hist_diagnostic_cone , TH2D *hist_diagnostic_cone_pT , TH2D *hist_diagnostic_cone_z , TH1D *histFF )
  { //parent jet , //vector of constits, // vector to fill with non-ghost pT , histo, histo, histo, first histo, pythia particle index

    Int_t num_all_constits = all_constits.size(); //get all the constituents (including ghosts) to loop over
    //Intializing variables, need to retrun first 6 (z_sum is for diagnostic purposes)
    Double_t angularity_sum =0; //return
    Double_t angularity_sum_nw = 0; //return
    Int_t num_non_ghosts =0; //return
    Double_t first_pythia_particle_index = 0; //return
    Double_t tru_sum=0; //return
    Double_t fake_sum=0; //return
    Bool_t got_pythia_particle;
    Double_t z_sum = 0;

    

    got_pythia_particle = kFALSE;
    for(Int_t i = 0 ; i < num_all_constits ; i++ ){
      if( all_constits[i].perp() > 1e-50 ){ //DO NOT GRAB GHOSTS !!!

        if( !got_pythia_particle ){
          if( all_constits[i].user_index() != -1 ){
            first_pythia_particle_index = all_constits[i].user_index();
            got_pythia_particle = kTRUE; //you found the pythia particle
          }
        }
      
        num_non_ghosts++;

        Double_t z = all_constits[i].perp() / parent_jet.pt();
        z_sum += z;
        histFF->Fill( z); //frag func

        Double_t relphi = all_constits[i].delta_phi_to(parent_jet);
        Double_t releta = parent_jet.eta() - all_constits[i].eta();
        Double_t rely = parent_jet.rap() - all_constits[i].rap();

        hist_diagnostic_cone->Fill(relphi,rely); //checking jet sub-structure with these diagnostics
        hist_diagnostic_cone_pT->Fill(relphi,rely,all_constits[i].perp());
        hist_diagnostic_cone_z->Fill(relphi,rely,z);

        Double_t delri = TMath::Sqrt( (relphi*relphi) + (releta*releta) );

        angularity_sum += (z*delri); //summing for the traditional angularity
        angularity_sum_nw += (delri); //summing the delri's to divide by the number of consituents (that are not ghosts)       

        if( all_constits[i].user_index() < 0){
          fake_sum = fake_sum+all_constits[i].perp();
        }
        else{
          tru_sum = tru_sum+all_constits[i].perp();
        } 

      } 
    }
    
    //cout<<"\n\nz_sum = "<<z_sum<<endl; //diganostic purposes , should always be 1 (or extremely close when excluding ghosts)

    static Double_t return_vals[4];
    return_vals[0] = angularity_sum; //angularity (traditional)
    return_vals[1] = angularity_sum_nw; //angularity (number-weighted)
    return_vals[2] = num_non_ghosts; //number of non-ghosts
    return_vals[3] = first_pythia_particle_index; //index of first pythia particle
    return_vals[4] = tru_sum; //sum of pythia particle momentum in jet
    return_vals[5] = fake_sum; //sum of bkgd particle momentum in jet
    return return_vals;
  };

 // exit(0);

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
  char expression15x[256];
  char expression15y[256];
  char expression15z[256];

  sprintf(expression2, "p_{T} distribution of all Pythia Particles with |#eta| < 0.9" );
  sprintf(expression3, "#eta distribution of all Pythia Particles with |#eta| < 0.9" );
  sprintf(expression4, "#phi distribution of all Pythia Particles with |#eta| < 0.9" );

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

  sprintf(expression15x, "Area of anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );
  sprintf(expression15y, "#rho = p_{T}^{jet}/Area of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );
  sprintf(expression15z, "Number-weighted Angularity of all anti-k_{T} Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, jet_pT_cut_high , 0.9 - Jet_Radius );


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

  TH1D *histangularity_pyth_jet = new TH1D("histangularity_pyth_jet", expression10,100,-0.5,Jet_Radius + 0.5);
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

  TH1D *histarea_pyth_jet = new TH1D("histarea_pyth_jet", expression15x,200,0.,2.);
  histarea_pyth_jet -> Sumw2();
  histarea_pyth_jet->SetXTitle("Area^{jet}");
  histarea_pyth_jet->SetYTitle("dN/dArea^{jet}");

  TH1D *histrho_pyth_jet = new TH1D("histrho_pyth_jet", expression15y,600,0.,600.);
  histrho_pyth_jet -> Sumw2();
  histrho_pyth_jet->SetXTitle("#rho^{jet} = p_{T}^{jet} / Area^{jet} (GeV/c)");
  histrho_pyth_jet->SetYTitle("dN/d#rho^{jet}");

  TH1D *histangularity_nw_pyth_jet = new TH1D("histangularity_nw_pyth_jet", expression15z,100,-0.5,Jet_Radius + 0.5);
  histangularity_nw_pyth_jet -> Sumw2();
  histangularity_nw_pyth_jet->SetXTitle("Angularity^{num. weighted}");
  histangularity_nw_pyth_jet->SetYTitle("dN/dAngularity");


  TH2D *hist_diagnostic_jet_cone_pyth_jet = new TH2D("hist_diagnostic_jet_cone_pyth_jet" , "Diagnostic Test of un-weighted distribution of particles in pythia only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pyth_jet -> Sumw2();
  hist_diagnostic_jet_cone_pyth_jet->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pyth_jet->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_pyth_jet_pT = new TH2D("hist_diagnostic_jet_cone_pyth_jet_pT" , "Diagnostic Test of p_{T}^{constit.}-weighted distribution of particles in pythia only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pyth_jet_pT -> Sumw2();
  hist_diagnostic_jet_cone_pyth_jet_pT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pyth_jet_pT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_pyth_jet_z = new TH2D("hist_diagnostic_jet_cone_pyth_jet_z" , "Diagnostic Test of z^{constit.}-weighted distribution of particles in pythia only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pyth_jet_z -> Sumw2();
  hist_diagnostic_jet_cone_pyth_jet_z->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pyth_jet_z->SetYTitle("#Deltay = y^{jet} - y^{constit.}");


//___________________________________BACKGROUND ONLY PARTICLES_________________________________//

  TH1D *histbackground_events = new TH1D( "histbackground_events" , "Number of background Events" , 1 , 0.5 , 1.5 );
  TH1D *histbackground_kT_jets = new TH1D( "histbackground_kT_jets" , "Number of k_{T} background Jets" , 1 , 0.5 , 1.5 );
  TH1D *histbackground_antikT_jets = new TH1D( "histbackground_antikT_jets" , "Number of anti-k_{T} background Jets" , 1 , 0.5 , 1.5 );

  char expression16[256];
  char expression17[256];
  char expression18[256];
  char expression19[256];
  char expression19a[256];
  char expression20[256];
  char expression20a[256];
  char expression21[256];
  char expression21a[256];
  char expression22[256];
  char expression22a[256];
  char expression23[256];
  char expression23a[256];
  char expression24[256];
  char expression24a[256];
  char expression25[256];
  char expression25a[256];
  char expression26[256];
  char expression26a[256];
  char expression27[256];
  char expression27a[256];
  char expression28[256];
  char expression28a[256];
  char expression29[256];
  char expression29a[256];

  char expression29ax[256];
  char expression29bx[256];
  char expression29ay[256];
  char expression29by[256];
  char expression29az[256];
  char expression29bz[256];

  char expression30[256];
  char expression31[256];
  char expression32[256];
  char expression33[256];
  char expression34[256];

  char expression35[256];

  sprintf(expression16, "p_{T} distribution of all Background Particles with |#eta| < 0.9" );
  sprintf(expression17, "#eta distribution of all Background Particles with |#eta| < 0.9" );
  sprintf(expression18, "#phi distribution of all Background Particles with |#eta| < 0.9" );

  sprintf(expression19, "p_{T} distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression19a, "p_{T} distribution of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression20, "#eta distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression20a, "#eta distribution of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression21, "#phi distribution of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression21a, "#phi distribution of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression22, "Number of tracks for all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression22a, "Number of tracks for all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius  );
  sprintf(expression23, "Mean track p_{T} (<p_{T}>) of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression23a, "Mean track p_{T} (<p_{T}>) of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius);
  sprintf(expression24, "Jet Angularity of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression24a, "Jet Angularity of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression25, "Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression25a, "Leading Track p_{T} of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression26, "Sub-Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression26a, "Sub-Leading Track p_{T} of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression27, "3rd Leading Track p_{T} of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression27a, "3rd Leading Track p_{T} of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression28, "Raw Fragmentation Function of all k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression28a, "Raw Fragmentation Function of all anti-k_{T} Background Jets %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29, "LeSub of all k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29a, "LeSub of all anti-k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

  sprintf(expression29ax, "Area of all k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29bx, "Area of all anti-k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29ay, "#rho of all k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29by, "#rho of all anti-k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29az, "Number-weighted angularity of all k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );
  sprintf(expression29bz, "Number-weighted angularity of all anti-k_{T} Background Jets, %lf - %lf GeV/c with |#eta| < %lf", constit_cut, 1000. , 0.9 - Jet_Radius );

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

  TH1D *histpT_bkgd_jet_kT = new TH1D("histpT_bkgd_jet_kT", expression19,100,0.,100.);
  histpT_bkgd_jet_kT -> Sumw2();
  histpT_bkgd_jet_kT->SetXTitle("p_{T} (GeV/c)");
  histpT_bkgd_jet_kT->SetYTitle("dN/dp_{T} (GeV/c)");

  TH1D *histpT_bkgd_jet_antikT = new TH1D("histpT_bkgd_jet_antikT", expression19a,100,0.,100.);
  histpT_bkgd_jet_antikT -> Sumw2();
  histpT_bkgd_jet_antikT->SetXTitle("p_{T} (GeV/c)");
  histpT_bkgd_jet_antikT->SetYTitle("dN/dp_{T} (GeV/c)");

  TH1D *histeta_bkgd_jet_kT = new TH1D("histeta_bkgd_jet_kT", expression20,100,-0.9,0.9);
  histeta_bkgd_jet_kT -> Sumw2();
  histeta_bkgd_jet_kT->SetXTitle("#eta^{jet} (pseudo-rap units)");
  histeta_bkgd_jet_kT->SetYTitle("dN/d#eta");

  TH1D *histeta_bkgd_jet_antikT = new TH1D("histeta_bkgd_jet_antikT", expression20a,100,-0.9,0.9);
  histeta_bkgd_jet_antikT -> Sumw2();
  histeta_bkgd_jet_antikT->SetXTitle("#eta^{jet} (pseudo-rap units)");
  histeta_bkgd_jet_antikT->SetYTitle("dN/d#eta");

  TH1D *histphi_bkgd_jet_kT = new TH1D("histphi_bkgd_jet_kT", expression21,100,0.,2.*TMath::Pi());
  histphi_bkgd_jet_kT -> Sumw2();
  histphi_bkgd_jet_kT->SetXTitle("#phi)^{jet} (radians)");
  histphi_bkgd_jet_kT->SetYTitle("dN/d#phi");

  TH1D *histphi_bkgd_jet_antikT = new TH1D("histphi_bkgd_jet_antikT", expression21a,100,0.,2.*TMath::Pi());
  histphi_bkgd_jet_antikT -> Sumw2();
  histphi_bkgd_jet_antikT->SetXTitle("#phi)^{jet} (radians)");
  histphi_bkgd_jet_antikT->SetYTitle("dN/d#phi");

//______________________________________________________________//

  TH1D *histnumtr_bkgd_jet_kT = new TH1D("histnumtr_bkgd_jet_kT", expression22,500,0,500);
  histnumtr_bkgd_jet_kT -> Sumw2();
  histnumtr_bkgd_jet_kT->SetXTitle("N_{tracks}");
  histnumtr_bkgd_jet_kT->SetYTitle("dN/dN_{tracks}");

  TH1D *histnumtr_bkgd_jet_antikT = new TH1D("histnumtr_bkgd_jet_antikT", expression22a,500,0,500);
  histnumtr_bkgd_jet_antikT -> Sumw2();
  histnumtr_bkgd_jet_antikT->SetXTitle("N_{tracks}");
  histnumtr_bkgd_jet_antikT->SetYTitle("dN/dN_{tracks}");

  TH1D *histmeanpTtr_bkgd_jet_kT = new TH1D("histmeanpTtr_bkgd_jet_kT", expression23,100,0,10);
  histmeanpTtr_bkgd_jet_kT -> Sumw2();
  histmeanpTtr_bkgd_jet_kT->SetXTitle("<p_{T}^{track}> (GeV/c)");
  histmeanpTtr_bkgd_jet_kT->SetYTitle("dN/d<p_{T}^{track}>");

  TH1D *histmeanpTtr_bkgd_jet_antikT = new TH1D("histmeanpTtr_bkgd_jet_antikT", expression23a,100,0,10);
  histmeanpTtr_bkgd_jet_antikT -> Sumw2();
  histmeanpTtr_bkgd_jet_antikT->SetXTitle("<p_{T}^{track}> (GeV/c)");
  histmeanpTtr_bkgd_jet_antikT->SetYTitle("dN/d<p_{T}^{track}>");

  TH1D *histangularity_bkgd_jet_kT = new TH1D("histangularity_bkgd_jet_kT", expression24,100,-0.5,Jet_Radius + 0.5);
  histangularity_bkgd_jet_kT -> Sumw2();
  histangularity_bkgd_jet_kT->SetXTitle("g (#phi - #eta distance units)");
  histangularity_bkgd_jet_kT->SetYTitle("dN/dg");

  TH1D *histangularity_bkgd_jet_antikT = new TH1D("histangularity_bkgd_jet_antikT", expression24a,100,-0.5,Jet_Radius + 0.5);
  histangularity_bkgd_jet_antikT -> Sumw2();
  histangularity_bkgd_jet_antikT->SetXTitle("g (#phi - #eta distance units)");
  histangularity_bkgd_jet_antikT->SetYTitle("dN/dg");

  TH1D *histldtr_bkgd_jet_kT = new TH1D("histldtr_bkgd_jet_kT", expression25,100,0.,10.);
  histldtr_bkgd_jet_kT -> Sumw2();
  histldtr_bkgd_jet_kT->SetXTitle("p_{T}^{ld. tr.} (GeV/c)");
  histldtr_bkgd_jet_kT->SetYTitle("dN/dp_{T}^{ld. tr.} ");

  TH1D *histldtr_bkgd_jet_antikT = new TH1D("histldtr_bkgd_jet_antikT", expression25a,100,0.,10.);
  histldtr_bkgd_jet_antikT -> Sumw2();
  histldtr_bkgd_jet_antikT->SetXTitle("p_{T}^{ld. tr.} (GeV/c)");
  histldtr_bkgd_jet_antikT->SetYTitle("dN/dp_{T}^{ld. tr.} ");

  TH1D *histsubldtr_bkgd_jet_kT = new TH1D("histsubldtr_bkgd_jet_kT", expression26,100,0.,10.);
  histsubldtr_bkgd_jet_kT -> Sumw2();
  histsubldtr_bkgd_jet_kT->SetXTitle("p_{T}^{sub-ld. tr.} (GeV/c)");
  histsubldtr_bkgd_jet_kT->SetYTitle("dN/dp_{T}^{sub-ld. tr.}");

  TH1D *histsubldtr_bkgd_jet_antikT = new TH1D("histsubldtr_bkgd_jet_antikT", expression26a,100,0.,10.);
  histsubldtr_bkgd_jet_antikT -> Sumw2();
  histsubldtr_bkgd_jet_antikT->SetXTitle("p_{T}^{sub-ld. tr.} (GeV/c)");
  histsubldtr_bkgd_jet_antikT->SetYTitle("dN/dp_{T}^{sub-ld. tr.}");

  TH1D *histthrdldtr_bkgd_jet_kT = new TH1D("histthrdldtr_bkgd_jet_kT", expression27,100,0.,10.);
  histthrdldtr_bkgd_jet_kT -> Sumw2();
  histthrdldtr_bkgd_jet_kT->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histthrdldtr_bkgd_jet_kT->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histthrdldtr_bkgd_jet_antikT = new TH1D("histthrdldtr_bkgd_jet_antikT", expression27a,100,0.,10.);
  histthrdldtr_bkgd_jet_antikT -> Sumw2();
  histthrdldtr_bkgd_jet_antikT->SetXTitle("p_{T}^{3rd-ld. tr.} (GeV/c");
  histthrdldtr_bkgd_jet_antikT->SetYTitle("dN/dp_{T}^{3rd-ld. tr.}");

  TH1D *histFF_bkgd_jet_kT = new TH1D("histFF_bkgd_jet_kT", expression28,100,0.,1.);
  histFF_bkgd_jet_kT -> Sumw2();
  histFF_bkgd_jet_kT->SetXTitle("z");
  histFF_bkgd_jet_kT->SetYTitle("dN/dz");

  TH1D *histFF_bkgd_jet_antikT = new TH1D("histFF_bkgd_jet_antikT", expression28a,100,0.,1.);
  histFF_bkgd_jet_antikT -> Sumw2();
  histFF_bkgd_jet_antikT->SetXTitle("z");
  histFF_bkgd_jet_antikT->SetYTitle("dN/dz");

  TH1D *histlesub_bkgd_jet_kT = new TH1D("histlesub_bkgd_jet_kT", expression29,100,0.,10.);
  histlesub_bkgd_jet_kT -> Sumw2();
  histlesub_bkgd_jet_kT->SetXTitle("LeSub = p_{T}^{ld. tr.} - p_{T}^{sub-ld. tr.} (GeV/c)");
  histlesub_bkgd_jet_kT->SetYTitle("dN/dLeSub");

  TH1D *histlesub_bkgd_jet_antikT = new TH1D("histlesub_bkgd_jet_antikT", expression29a,100,0.,10.);
  histlesub_bkgd_jet_antikT -> Sumw2();
  histlesub_bkgd_jet_antikT->SetXTitle("LeSub = p_{T}^{ld. tr.} - p_{T}^{sub-ld. tr.} (GeV/c)");
  histlesub_bkgd_jet_antikT->SetYTitle("dN/dLeSub");



  TH1D *histarea_bkgd_jet_kT = new TH1D("histarea_bkgd_jet_kT", expression29ax,200,0.,2.);
  histarea_bkgd_jet_kT -> Sumw2();
  histarea_bkgd_jet_kT->SetXTitle("Area^{jet}");
  histarea_bkgd_jet_kT->SetYTitle("dN/dArea^{jet}");

  TH1D *histarea_bkgd_jet_antikT = new TH1D("histarea_bkgd_jet_antikT", expression29bx,200,0.,2.);
  histarea_bkgd_jet_antikT -> Sumw2();
  histarea_bkgd_jet_antikT->SetXTitle("Area^{jet}");
  histarea_bkgd_jet_antikT->SetYTitle("dN/dArea^{jet}");

  TH1D *histrho_bkgd_jet_kT = new TH1D("histrho_bkgd_jet_kT", expression29ay,600,0.,600.);
  histrho_bkgd_jet_kT -> Sumw2();
  histrho_bkgd_jet_kT->SetXTitle("#rho^{jet} = p_{T}^{jet} / Area^{jet} (GeV/c)");
  histrho_bkgd_jet_kT->SetYTitle("dN/d#rho^{jet}");

  TH1D *histrho_bkgd_jet_antikT = new TH1D("histrho_bkgd_jet_antikT", expression29by,600,0.,600.);
  histrho_bkgd_jet_antikT -> Sumw2();
  histrho_bkgd_jet_antikT->SetXTitle("#rho^{jet} = p_{T}^{jet} / Area^{jet} (GeV/c)");
  histrho_bkgd_jet_antikT->SetYTitle("dN/d#rho^{jet}");

  TH1D *histangularity_nw_bkgd_jet_kT = new TH1D("histangularity_nw_bkgd_jet_kT", expression29az,100,-0.5,Jet_Radius + 0.5);
  histangularity_nw_bkgd_jet_kT -> Sumw2();
  histangularity_nw_bkgd_jet_kT->SetXTitle("Angularity^{num. weighted}");
  histangularity_nw_bkgd_jet_kT->SetYTitle("dN/dAngularity");

  TH1D *histangularity_nw_bkgd_jet_antikT = new TH1D("histangularity_nw_bkgd_jet_antikT", expression29bz,100,-0.5,Jet_Radius + 0.5);
  histangularity_nw_bkgd_jet_antikT -> Sumw2();
  histangularity_nw_bkgd_jet_antikT->SetXTitle("Angularity^{num. weighted}");
  histangularity_nw_bkgd_jet_antikT->SetYTitle("dN/dAngularity");




  TH2D *hist_diagnostic_jet_cone_bkgd_jet_kT = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_kT" , "Diagnostic Test of distribution of un-weighted particles in background only k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_kT -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_kT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_kT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_bkgd_jet_kT_pT = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_kT_pT" , "Diagnostic Test of p_{T}^{constit.}-weighted distribution of particles in background only k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_kT_pT -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_kT_pT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_kT_pT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_bkgd_jet_kT_z = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_kT_z" , "Diagnostic Test of z^{constit.}-weighted distribution of particles in background only k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_kT_z -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_kT_z->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_kT_z->SetYTitle("#Deltay = y^{jet} - y^{constit.}");




  TH2D *hist_diagnostic_jet_cone_bkgd_jet_antikT = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_antikT" , "Diagnostic Test of distribution of particles in background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_antikT -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_antikT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_antikT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_bkgd_jet_antikT_pT = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_antikT_pT" , "Diagnostic Test of p_{T}^{constit.}-weighted distribution of particles in background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_antikT_pT -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_antikT_pT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_antikT_pT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_bkgd_jet_antikT_z = new TH2D("hist_diagnostic_jet_cone_bkgd_jet_antikT_z" , "Diagnostic Test of z^{constit.}-weighted distribution of particles in background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_bkgd_jet_antikT_z -> Sumw2();
  hist_diagnostic_jet_cone_bkgd_jet_antikT_z->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_bkgd_jet_antikT_z->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

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
  char expression50x[256];
  char expression50y[256];
  char expression50z[256];
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
  sprintf(expression50x, "Area of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression50y, "#rho of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );
  sprintf(expression50z, "Number-weighted Angularity of all anti-k_{T} Background + Pythia Jets %lf - %lf GeV/c with |#eta| < %lf", jet_pT_cut_low, 100. , 0.9 - Jet_Radius );

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

  TH1D *histangularity_pytha_AND_bkgd_jet = new TH1D("histangularity_pytha_AND_bkgd_jet", expression45,100,-0.5,Jet_Radius + 0.5);
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


  TH1D *histarea_pytha_AND_bkgd_jet = new TH1D("histarea_pytha_AND_bkgd_jet", expression50x,200,0.,2.);
  histarea_pytha_AND_bkgd_jet -> Sumw2();
  histarea_pytha_AND_bkgd_jet->SetXTitle("Area^{jet}");
  histarea_pytha_AND_bkgd_jet->SetYTitle("dN/dArea^{jet}");

  TH1D *histrho_pytha_AND_bkgd_jet = new TH1D("histrho_pytha_AND_bkgd_jet", expression50y,600,0.,600.);
  histrho_pytha_AND_bkgd_jet -> Sumw2();
  histrho_pytha_AND_bkgd_jet->SetXTitle("#rho^{jet} = p_{T}^{jet} / Area^{jet} (GeV/c)");
  histrho_pytha_AND_bkgd_jet->SetYTitle("dN/d#rho^{jet}");

  TH1D *histangularity_nw_pytha_AND_bkgd_jet = new TH1D("histangularity_nw_pytha_AND_bkgd_jet", expression50z,100,-0.5,Jet_Radius + 0.5);
  histangularity_nw_pytha_AND_bkgd_jet -> Sumw2();
  histangularity_nw_pytha_AND_bkgd_jet->SetXTitle("Angularity^{num. weighted}");
  histangularity_nw_pytha_AND_bkgd_jet->SetYTitle("dN/dAngularity");


  TH2D *hist_diagnostic_jet_cone_pythia_AND_bkgd_jet = new TH2D("hist_diagnostic_jet_cone_pythia_AND_bkgd_jet" , "Diagnostic Test of distribution of particles in pythia + background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet -> Sumw2();
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT = new TH2D("hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT" , "Diagnostic Test of p_{T}^{constit.}-weighted distribution of particles in pythia + background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT -> Sumw2();
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

  TH2D *hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z = new TH2D("hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z" , "Diagnostic Test of z^{constit.}-weighted distribution of particles in pythia + background only anti-k_{T} jets",200,-1.0,1.0,200,-1.0,1.0);
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z -> Sumw2();
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z->SetXTitle("#Delta#phi = #phi^{jet} - #phi^{constit.}");
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z->SetYTitle("#Deltay = y^{jet} - y^{constit.}");

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

  //Int_t* time_array;
  //Int_t* time_interval_array;
  //Int_t* event_array;

  std::vector<Int_t> time_array;
  //time_array = new Int_t[nEvents + 1];


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
  jetDef = new fastjet::JetDefinition(algorithm,Rparam,recombScheme,strategy); //anti-kT
  jetDef2 = new fastjet::JetDefinition(algorithm2,Rparam,recombScheme,strategy); //kT
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
  time_array.push_back((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;

  Int_t user_index = 1; //must start @ 1 (not 0 !)

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

  TFile *ff = new TFile(expression1,"RECREATE"); //make output files and directories
  TDirectory *pythia_only = ff->mkdir("pythia_only");
  TDirectory *background_only = ff->mkdir("background_only");
  TDirectory *background_and_pythia = ff->mkdir("background_and_pythia");
  TDirectory *Run_Time_Stats = ff->mkdir("Run_Time_Stats");


//____________________________________________________________________________________________________________________//

//Patrick these are the lines for the TTree Christine wanted

  char treestr_pyth[128];
  sprintf( treestr_pyth , "Pythia");
  TTree *tree_pyth = new TTree(treestr_pyth,"tree_pyth");

  //making variables for the tree branches (pythia jet properties)
  Float_t p_T_tree_pyth; //0th
  Float_t Eta_tree_pyth; //1st
  Float_t Phi_tree_pyth; //2nd
  Float_t Area_tree_pyth; //3rd
  Float_t Rho_tree_pyth; //4th
  Float_t p_T_corr_tree_pyth; //5th
  Float_t N_Trk_tree_pyth; //6th
  Float_t Angularity_tree_pyth; //7th
  Float_t Angularity_NW_tree_pyth; //8th
  Float_t Mean_p_T_tree_pyth; //9th
  Float_t p_T_1_tree_pyth; //10th
  Float_t p_T_2_tree_pyth; //11th
  Float_t p_T_3_tree_pyth; //12th
  Float_t p_T_4_tree_pyth; //13th
  Float_t p_T_5_tree_pyth; //14th
  Float_t X_tru_tree_pyth; //15th


  //Now add these branches to the trees

   tree_pyth->Branch("pT",&p_T_tree_pyth,"pT/F");
   tree_pyth->Branch("Eta",&Eta_tree_pyth,"Eta/F");
   tree_pyth->Branch("Phi",&Phi_tree_pyth,"Phi/F");
   tree_pyth->Branch("Area",&Area_tree_pyth,"Area/F");
   tree_pyth->Branch("Rho",&Rho_tree_pyth,"Rho/F");
   tree_pyth->Branch("pTcorr",&p_T_corr_tree_pyth,"pTcorr/F");
   tree_pyth->Branch("NTrk",&N_Trk_tree_pyth,"NTrk/F");
   tree_pyth->Branch("Angularity",&Angularity_tree_pyth,"Angularity/F");
   tree_pyth->Branch("Angularity_NW",&Angularity_NW_tree_pyth,"Angularity_NW/F");
   tree_pyth->Branch("MeanpT",&Mean_p_T_tree_pyth,"MeanpT/F");
   tree_pyth->Branch("pT1",&p_T_1_tree_pyth,"pT1/F");
   tree_pyth->Branch("pT2",&p_T_2_tree_pyth,"pT2/F");
   tree_pyth->Branch("pT3",&p_T_3_tree_pyth,"pT3/F");
   tree_pyth->Branch("pT4",&p_T_4_tree_pyth,"pT4/F");
   tree_pyth->Branch("pT5",&p_T_5_tree_pyth,"pT5/F");
   tree_pyth->Branch("XTru",&X_tru_tree_pyth,"XTru/F");




//_______________________________________________________//


  //___________________kT below______________________//
  char bkgd_kt_treestr[128];
  sprintf( bkgd_kt_treestr , "TennGen-kT-HF-%d-CB-%d" , HF, cent_bin);
  TTree *tree_bkgd_kt = new TTree(bkgd_kt_treestr,"tree_bkgd_kt");

  //making variables for the tree branches (TennGen kT jet properties)
  Float_t p_T_tree_bkgd; //0th
  Float_t Eta_tree_bkgd; //1st
  Float_t Phi_tree_bkgd; //2nd
  Float_t Area_tree_bkgd; //3rd
  Float_t Rho_tree_bkgd;  //4th
  Float_t p_T_corr_tree_bkgd; //5th
  Float_t N_Trk_tree_bkgd; //6th
  Float_t Angularity_tree_bkgd; //7th
  Float_t Angularity_NW_tree_bkgd; //8th
  Float_t Mean_p_T_tree_bkgd; //9th
  Float_t p_T_1_tree_bkgd; //10th
  Float_t p_T_2_tree_bkgd; //11th
  Float_t p_T_3_tree_bkgd; //12th
  Float_t p_T_4_tree_bkgd; //13th
  Float_t p_T_5_tree_bkgd; //14th
  Float_t X_tru_tree_bkgd; //15th


  //Now add these branches to the trees

   tree_bkgd_kt->Branch("pT",&p_T_tree_bkgd,"pT/F");
   tree_bkgd_kt->Branch("Eta",&Eta_tree_bkgd,"Eta/F");
   tree_bkgd_kt->Branch("Phi",&Phi_tree_bkgd,"Phi/F");
   tree_bkgd_kt->Branch("Area",&Area_tree_bkgd,"Area/F");
   tree_bkgd_kt->Branch("Rho",&Rho_tree_bkgd,"Rho/F");
   tree_bkgd_kt->Branch("pTcorr",&p_T_corr_tree_bkgd,"pTcorr/F");
   tree_bkgd_kt->Branch("NTrk",&N_Trk_tree_bkgd,"NTrk/F");
   tree_bkgd_kt->Branch("Angularity",&Angularity_tree_bkgd,"Angularity/F");
   tree_bkgd_kt->Branch("Angularity_NW",&Angularity_NW_tree_bkgd,"Angularity_NW/F");
   tree_bkgd_kt->Branch("MeanpT",&Mean_p_T_tree_bkgd,"MeanpT/F");
   tree_bkgd_kt->Branch("pT1",&p_T_1_tree_bkgd,"pT1/F");
   tree_bkgd_kt->Branch("pT2",&p_T_2_tree_bkgd,"pT2/F");
   tree_bkgd_kt->Branch("pT3",&p_T_3_tree_bkgd,"pT3/F");
   tree_bkgd_kt->Branch("pT4",&p_T_4_tree_bkgd,"pT4/F");
   tree_bkgd_kt->Branch("pT5",&p_T_5_tree_bkgd,"pT5/F");
   tree_bkgd_kt->Branch("XTru",&X_tru_tree_bkgd,"XTru/F");

   //_______________antikT below______________________//
  char bkgd_antikt_treestr[128];
  sprintf( bkgd_antikt_treestr , "TennGen-antikT-HF-%d-CB-%d" , HF, cent_bin);
  TTree *tree_bkgd_antikt = new TTree(bkgd_antikt_treestr,"tree_bkgd_antikt");

  //making variables for the tree branches (TennGen kT jet properties)
  Float_t p_T_tree_antikt_bkgd; //0th
  Float_t Eta_tree_antikt_bkgd; //1st
  Float_t Phi_tree_antikt_bkgd; //2nd
  Float_t Area_tree_antikt_bkgd; //3rd
  Float_t Rho_tree_antikt_bkgd;  //4th
  Float_t p_T_corr_tree_antikt_bkgd; //5th
  Float_t N_Trk_tree_antikt_bkgd; //6th
  Float_t Angularity_tree_antikt_bkgd; //7th
  Float_t Angularity_NW_tree_antikt_bkgd; //8th
  Float_t Mean_p_T_tree_antikt_bkgd; //9th
  Float_t p_T_1_tree_antikt_bkgd; //10th
  Float_t p_T_2_tree_antikt_bkgd; //11th
  Float_t p_T_3_tree_antikt_bkgd; //12th
  Float_t p_T_4_tree_antikt_bkgd; //13th
  Float_t p_T_5_tree_antikt_bkgd; //14th
  Float_t X_tru_tree_antikt_bkgd; //15th


  //Now add these branches to the trees

   tree_bkgd_antikt->Branch("pT",&p_T_tree_antikt_bkgd,"pT/F");
   tree_bkgd_antikt->Branch("Eta",&Eta_tree_antikt_bkgd,"Eta/F");
   tree_bkgd_antikt->Branch("Phi",&Phi_tree_antikt_bkgd,"Phi/F");
   tree_bkgd_antikt->Branch("Area",&Area_tree_antikt_bkgd,"Area/F");
   tree_bkgd_antikt->Branch("Rho",&Rho_tree_antikt_bkgd,"Rho/F");
   tree_bkgd_antikt->Branch("pTcorr",&p_T_corr_tree_antikt_bkgd,"pTcorr/F");
   tree_bkgd_antikt->Branch("NTrk",&N_Trk_tree_antikt_bkgd,"NTrk/F");
   tree_bkgd_antikt->Branch("Angularity",&Angularity_tree_antikt_bkgd,"Angularity/F");
   tree_bkgd_antikt->Branch("Angularity_NW",&Angularity_NW_tree_antikt_bkgd,"Angularity_NW/F");
   tree_bkgd_antikt->Branch("MeanpT",&Mean_p_T_tree_antikt_bkgd,"MeanpT/F");
   tree_bkgd_antikt->Branch("pT1",&p_T_1_tree_antikt_bkgd,"pT1/F");
   tree_bkgd_antikt->Branch("pT2",&p_T_2_tree_antikt_bkgd,"pT2/F");
   tree_bkgd_antikt->Branch("pT3",&p_T_3_tree_antikt_bkgd,"pT3/F");
   tree_bkgd_antikt->Branch("pT4",&p_T_4_tree_antikt_bkgd,"pT4/F");
   tree_bkgd_antikt->Branch("pT5",&p_T_5_tree_antikt_bkgd,"pT5/F");
   tree_bkgd_antikt->Branch("XTru",&X_tru_tree_antikt_bkgd,"XTru/F");




//___________________________________________________________________________________//
  char treestr[128];
  sprintf( treestr , "Pythia-and-TennGen-HF-%d-CB-%d" , HF, cent_bin);
  TTree *tree = new TTree(treestr,"TreeID");

  //making variables for the tree branches (jet properties)
  Float_t p_T_tree; //0th
  Float_t Eta_tree; //1st
  Float_t Phi_tree; //2nd
  Float_t Area_tree; //3rd
  Float_t Rho_tree;  //4th
  Float_t p_T_corr_tree; //5th
  Float_t N_Trk_tree; //6th
  Float_t Angularity_tree; //7th
  Float_t Angularity_NW_tree; //8th
  Float_t Mean_p_T_tree; //9th
  Float_t p_T_1_tree; //10th
  Float_t p_T_2_tree; //11th
  Float_t p_T_3_tree; //12th
  Float_t p_T_4_tree; //13th
  Float_t p_T_5_tree; //14th
  Float_t X_tru_tree; //15th


  //Now add these branches to the trees

   tree->Branch("pT",&p_T_tree,"pT/F");
   tree->Branch("Eta",&Eta_tree,"Eta/F");
   tree->Branch("Phi",&Phi_tree,"Phi/F");
   tree->Branch("Area",&Area_tree,"Area/F");
   tree->Branch("Rho",&Rho_tree,"Rho/F");
   tree->Branch("pTcorr",&p_T_corr_tree,"pTcorr/F");
   tree->Branch("NTrk",&N_Trk_tree,"NTrk/F");
   tree->Branch("Angularity",&Angularity_tree,"Angularity/F");
   tree->Branch("Angularity_NW",&Angularity_NW_tree,"Angularity_NW/F");
   tree->Branch("MeanpT",&Mean_p_T_tree,"MeanpT/F");
   tree->Branch("pT1",&p_T_1_tree,"pT1/F");
   tree->Branch("pT2",&p_T_2_tree,"pT2/F");
   tree->Branch("pT3",&p_T_3_tree,"pT3/F");
   tree->Branch("pT4",&p_T_4_tree,"pT4/F");
   tree->Branch("pT5",&p_T_5_tree,"pT5/F");
   tree->Branch("XTru",&X_tru_tree,"XTru/F");


//___________________________________________________________________________________//
  char treestr_pat[128];
  sprintf( treestr_pat , "PATS-ALT-Pythia-and-TennGen-HF-%d-CB-%d" , HF, cent_bin);
  TTree *tree_pat = new TTree(treestr_pat,"TreeID");

  //making variables for the tree branches (jet properties)
  Float_t p_T_tree_pat; //0th
  Float_t Eta_tree_pat; //1st
  Float_t Phi_tree_pat; //2nd
  Float_t Area_tree_pat; //3rd
  Float_t Rho_tree_pat;  //4th
  Float_t p_T_corr_tree_pat; //5th
  Float_t N_Trk_tree_pat; //6th
  Float_t Angularity_tree_pat; //7th
  Float_t Angularity_NW_tree_pat; //8th
  Float_t Mean_p_T_tree_pat; //9th
  Float_t p_T_1_tree_pat; //10th
  Float_t p_T_2_tree_pat; //11th
  Float_t p_T_3_tree_pat; //12th
  Float_t p_T_4_tree_pat; //13th
  Float_t p_T_5_tree_pat; //14th
  Float_t geom_match_tree_pat; //15th
  Float_t mom_frac_match_tree_pat; //16th
  Float_t X_tru_tree_pat; //17th


  //Now add these branches to the trees

   tree_pat->Branch("pT",&p_T_tree,"pT/F");
   tree_pat->Branch("Eta",&Eta_tree,"Eta/F");
   tree_pat->Branch("Phi",&Phi_tree,"Phi/F");
   tree_pat->Branch("Area",&Area_tree,"Area/F");
   tree_pat->Branch("Rho",&Rho_tree,"Rho/F");
   tree_pat->Branch("pTcorr",&p_T_corr_tree,"pTcorr/F");
   tree_pat->Branch("NTrk",&N_Trk_tree,"NTrk/F");
   tree_pat->Branch("Angularity",&Angularity_tree,"Angularity/F");
   tree_pat->Branch("Angularity_NW",&Angularity_NW_tree,"Angularity_NW/F");
   tree_pat->Branch("MeanpT",&Mean_p_T_tree,"MeanpT/F");
   tree_pat->Branch("pT1",&p_T_1_tree,"pT1/F");
   tree_pat->Branch("pT2",&p_T_2_tree,"pT2/F");
   tree_pat->Branch("pT3",&p_T_3_tree,"pT3/F");
   tree_pat->Branch("pT4",&p_T_4_tree,"pT4/F");
   tree_pat->Branch("pT5",&p_T_5_tree,"pT5/F");
   tree_pat->Branch("distmatch",&geom_match_tree_pat,"distmatch");
   tree_pat->Branch("XMatch",&mom_frac_match_tree_pat,"XMatch");
   tree_pat->Branch("XTru",&X_tru_tree,"XTru/F");


//____________________________________________________________________________________________________________________//

  char filepathstr[512];

  if( !GRID ){

    //PUT THE PATH TO THE BKGD_ROOT_FILES directory here !!!!!!!//
    sprintf(filepathstr,"/home/alidock/ML/BKGD_ROOT_FILES");
    //sprintf(filepathstr,"/home/charles/Documents/research/Background_Research/forcharles/Newest_Background_Code_05_15_2018/Harmonic_Code_for_copying/DiJet_Asymmetry/Updated_Code/Latest_most_up_to_date_Background_Code_08_31_2018/Latest_Version_meant_for_anti_kT/Frag_Func_Code/Heavy_Ion_BGLoad/Patrick_Studies/ML-Jet-BG-Subtraction/BKGD_ROOT_FILES");

  }
  else if( GRID ){

   getcwd(filepathstr, sizeof(filepathstr));

  }

  bkgd2->PassInSettings(filepathstr,  0 , 0 , 0 , seed6 );
  bkgd2->Load();
  bkgd2->PrintStatus();
  TClonesArray *BKGD = (TClonesArray *)bkgd2->GetEventBG(); //(outside particle loop)

  
  std::vector <Double_t> rho_vec; //vector of rhos to get median from background only;
  std::vector <Double_t> rho_vec_TOTAL; //vector of rhos to get median from TOTAL while excluding 2 leading jets;

  remove( "pbgOut.csv" ); //if the csv file alread exists, delete it.

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
        time_array.push_back((time->GetYear()*31540000) + (time->GetMonth()*2628000) + (time->GetDay()*86400) + (time->GetHour()*3600) + (time->GetMinute()*60) + time->GetSecond()) ;
        Int_t user_index = 1; //must reset this each event iteration loop
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
     
//      if(pt > constit_cut && pt < jet_pT_cut_high && (TMath::Abs(p_eta)) < 0.9  ){ //looking for particles in the EMCAL with pT = 0.15 GeV/c (c = 1)   
      if(pt > constit_cut && (TMath::Abs(p_eta)) < 0.9  ){ //looking for particles in the EMCAL with pT = 0.15 GeV/c (c = 1) 
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
    fastjet::Selector select_pt3 = fastjet::SelectorPtRange(jet_pT_cut_low,1000); // Selects Jets with transverse momentum between jet_pT_cut_low and 100 GeV
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
        p_T_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].pt();
        histeta_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].eta() );
        Eta_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].eta();
        histphi_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].phi() );        
        Phi_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].phi();

        histarea_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].area() );
        histrho_pyth_jet->Fill( selected_jetsPythia_sorted_part[py_jet_ind2].pt()/selected_jetsPythia_sorted_part[py_jet_ind2].area() );

        Area_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].area();
        Rho_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].pt()/selected_jetsPythia_sorted_part[py_jet_ind2].area() ;

        Double_t jet_phi = selected_jetsPythia_sorted_part[py_jet_ind2].phi();
        Double_t jet_eta = selected_jetsPythia_sorted_part[py_jet_ind2].eta();

        Double_t jet_pt = selected_jetsPythia_sorted_part[py_jet_ind2].pt();
        p_T_corr_tree_pyth = selected_jetsPythia_sorted_part[py_jet_ind2].pt(); //ALL PYTHIA so no correction (correction = 0)

        std::vector <fastjet::PseudoJet> constituents_part = selected_jetsPythia_sorted_part[py_jet_ind2].constituents(); //grab the constituents
        std::vector <Double_t> pythia_constit_pT_vec;

        Int_t nPart_part = constituents_part.size();
        for( Int_t i_pyth = 0 ; i_pyth < nPart_part ; i_pyth++ ){
          if( constituents_part[i_pyth].perp() > 1e-50 ){
            pythia_constit_pT_vec.push_back( constituents_part[i_pyth].perp() );
          }
        }
       
        Double_t *return_arr_pyth;

        return_arr_pyth = jet_constit_loop( selected_jetsPythia_sorted_part[py_jet_ind2] , constituents_part , hist_diagnostic_jet_cone_pyth_jet , hist_diagnostic_jet_cone_pyth_jet_pT , hist_diagnostic_jet_cone_pyth_jet_z , histFF_pyth_jet ); 

        Double_t angularity_sum = *(return_arr_pyth + 0);
        Double_t angularity_sum_nw = *(return_arr_pyth + 1);
        Double_t num_const_gt_cut = *(return_arr_pyth + 2);
        Int_t first_pythia_particle_index = *(return_arr_pyth + 3);
        
/* Debug lines below
       if(angularity_sum > Jet_Radius){
         cout<<"\n\n______________________________________________________________"<<endl;
         cout<<"\n\nangularity was too big for this jet, angularity = "<<angularity_sum<<endl;     
         cout<<"\n\n______________________________________________________________"<<endl;
         exit(0);
       }
*/
        histangularity_pyth_jet->Fill(angularity_sum); //angularity
        histangularity_nw_pyth_jet->Fill(angularity_sum_nw/num_const_gt_cut); //number weighted angularity

        Angularity_tree_pyth = angularity_sum;
        Angularity_NW_tree_pyth = angularity_sum_nw/num_const_gt_cut ;

        histnumtr_pyth_jet->Fill(num_const_gt_cut); //number of consituents
        N_Trk_tree_pyth = num_const_gt_cut;

        constituents_part.clear();

        if(pythia_constit_pT_vec.size() > 0){
          
          histmeanpTtr_pyth_jet->Fill( std::accumulate(pythia_constit_pT_vec.begin(), pythia_constit_pT_vec.end(), 0.0) / pythia_constit_pT_vec.size() ); //mean track pT in the jet
          Mean_p_T_tree_pyth = std::accumulate(pythia_constit_pT_vec.begin(), pythia_constit_pT_vec.end(), 0.0) / pythia_constit_pT_vec.size() ;

          std::sort(pythia_constit_pT_vec.begin(), pythia_constit_pT_vec.end()); //sorts least to greatest
          histldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size()-1]); // grab leading track in jet

          if(pythia_constit_pT_vec.size() > 1){
            histsubldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 2] ); // grab sub-leading track in jet
            histlesub_pyth_jet->Fill( pythia_constit_pT_vec[pythia_constit_pT_vec.size()-1] - pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 2] );
            if(pythia_constit_pT_vec.size() > 2){
              histthrdldtr_pyth_jet->Fill(pythia_constit_pT_vec[pythia_constit_pT_vec.size() - 3]); //grad sub-sub leading track in jet
            } //for sub sub leading
          } //for sub leading

        p_T_1_tree_pyth = pythia_constit_pT_vec[pythia_constit_pT_vec.size()-1];
        p_T_2_tree_pyth = pythia_constit_pT_vec[pythia_constit_pT_vec.size()-2];
        p_T_3_tree_pyth = pythia_constit_pT_vec[pythia_constit_pT_vec.size()-3];
        p_T_4_tree_pyth = pythia_constit_pT_vec[pythia_constit_pT_vec.size()-4];
        p_T_5_tree_pyth = pythia_constit_pT_vec[pythia_constit_pT_vec.size()-5];

        } //make sure there are constits

      X_tru_tree_pyth = 1;
      pythia_constit_pT_vec.clear();
      tree_pyth->Fill();
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
/*
      if( num_background_events == 1){ // you have not yet called the background, therefore not filled the XL file !!!
        //cout<<"\n\nyou want to remove the pre-existingXL file!!!!"<<endl;
        remove( "pbgOut.csv" );
        //exit(0);
      }
*/
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
      fastjet::ClusterSequenceArea clustSeqTOTAL(fjInputs_TOTAL, *jetDef , *AreaDef );  //anti-kT

      vector <fastjet::PseudoJet> inclusiveJetsBackground, sortedJetsBackground; // 
      fastjet::ClusterSequenceArea clustSeqBackground(fjInputs_Background, *jetDef2 , *AreaDef ); // kT

      vector <fastjet::PseudoJet> inclusiveJetsBackground_antikT, sortedJetsBackground_antikT; // 
      fastjet::ClusterSequenceArea clustSeqBackground_antikT(fjInputs_Background, *jetDef , *AreaDef ); // antikT
      //fastjet::ClusterSequence clustSeqBackground(fjInputs_Background, *jetDef2 ); // 
      fjInputs_Background.clear();

      vector <fastjet::PseudoJet> inclusiveJetsTOTAL_kT, sortedJets_TOTAL_kT; // 
      fastjet::ClusterSequenceArea clustSeqTOTAL_kT(fjInputs_TOTAL_kT, *jetDef2 , *AreaDef ); // kT

      ////////////////////////////////////////////
      inclusiveJetsTOTAL = clustSeqTOTAL.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      ////////////////////////////////////////////

      ////////////////////////////////////////////
      inclusiveJetsBackground = clustSeqBackground.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      inclusiveJetsTOTAL_kT = clustSeqTOTAL_kT.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      inclusiveJetsBackground_antikT = clustSeqBackground_antikT.inclusive_jets(0.0); //only jets with pt greater than 0.0 GeV
      ////////////////////////////////////////////


      //sorting
      ////////////////////////////////////////////
      sortedJetsTOTAL   = sorted_by_pt(inclusiveJetsTOTAL); //sort by decreasing transverse momentum
      sortedJetsBackground   = sorted_by_pt(inclusiveJetsBackground); //sort by decreasing transverse momentum
      sortedJets_TOTAL_kT = sorted_by_pt(inclusiveJetsTOTAL_kT); //sort by decreasing transverse momentum
      sortedJetsBackground_antikT = sorted_by_pt( inclusiveJetsBackground_antikT ); //sort by decreasing transverse momentum
      ////////////////////////////////////////////

      vector<fastjet::PseudoJet>selected_jetsTOTAL = select_both3(sortedJetsTOTAL); 
      vector<fastjet::PseudoJet>selected_jetsTOTAL_sorted = sorted_by_pt(selected_jetsTOTAL);

      vector<fastjet::PseudoJet>selected_jetsBackground = select_both3(sortedJetsBackground); 
      vector<fastjet::PseudoJet>selected_jetsBackground_sorted = sorted_by_pt(selected_jetsBackground);

      vector<fastjet::PseudoJet>selected_jetsTOTAL_kT = select_both3(sortedJets_TOTAL_kT); 
      vector<fastjet::PseudoJet>selected_jetsTOTAL_kT_sorted = sorted_by_pt(selected_jetsTOTAL_kT);

      vector<fastjet::PseudoJet>selected_jetsBackground_antikT = select_both3(sortedJetsBackground_antikT); 
      vector<fastjet::PseudoJet>selected_jetsBackground_antikT_sorted = sorted_by_pt(selected_jetsBackground_antikT );


      Double_t running_rho_average = 0;

      cout<<"\n\n\n|||||||||||The total kT background jet population is = "<<selected_jetsBackground_sorted.size()<<"|||||||||||||||||"<<endl;

      Double_t event_median_shared;

/////////////////loop over the kT background only jets///////////////////////////////////////////////

      if(selected_jetsBackground_sorted.size() > 0){ 
        for(unsigned i_jet = 0 ; i_jet < selected_jetsBackground_sorted.size() ; i_jet++ ){ 

          histbackground_kT_jets->Fill(1.00);
    
          histpT_bkgd_jet_kT -> Fill(selected_jetsBackground_sorted[i_jet].pt());
          p_T_tree_bkgd = selected_jetsBackground_sorted[i_jet].pt();

          histeta_bkgd_jet_kT -> Fill(selected_jetsBackground_sorted[i_jet].eta());
          Eta_tree_bkgd = selected_jetsBackground_sorted[i_jet].eta();

          histphi_bkgd_jet_kT -> Fill(selected_jetsBackground_sorted[i_jet].phi());
          Phi_tree_bkgd = selected_jetsBackground_sorted[i_jet].phi();

          Area_tree_bkgd = selected_jetsBackground_sorted[i_jet].area();
          Rho_tree_bkgd = selected_jetsBackground_sorted[i_jet].pt()/selected_jetsBackground_sorted[i_jet].area() ;

          Double_t jet_pT_back = selected_jetsBackground_sorted[i_jet].pt();
          Double_t jet_eta_back = selected_jetsBackground_sorted[i_jet].eta();
          Double_t jet_phi_back = selected_jetsBackground_sorted[i_jet].phi();

          p_T_corr_tree_bkgd = 0; //the bkgd jets are all background so we are TEMPORARILY setting these to 0

          histarea_bkgd_jet_kT->Fill( selected_jetsBackground_sorted[i_jet].area() );
          histrho_bkgd_jet_kT->Fill( selected_jetsBackground_sorted[i_jet].pt()/selected_jetsBackground_sorted[i_jet].area() );
          rho_vec.push_back(selected_jetsBackground_sorted[i_jet].pt()/selected_jetsBackground_sorted[i_jet].area());
      
          vector <fastjet::PseudoJet> constituents_bkgd = selected_jetsBackground_sorted[i_jet].constituents();
          vector <Double_t> bkgd_constit_pT_vec;

          Int_t nPart_bkgd = constituents_bkgd.size();
          for( Int_t i_bkgd = 0 ; i_bkgd < nPart_bkgd ; i_bkgd++ ){
            if( constituents_bkgd[i_bkgd].perp() > 1e-50 ){
              bkgd_constit_pT_vec.push_back( constituents_bkgd[i_bkgd].perp() );
            }
          }

          Double_t *return_arr_bkgd;

          return_arr_bkgd = jet_constit_loop( selected_jetsBackground_sorted[i_jet] , constituents_bkgd , hist_diagnostic_jet_cone_bkgd_jet_kT , hist_diagnostic_jet_cone_bkgd_jet_kT_pT , hist_diagnostic_jet_cone_bkgd_jet_kT_z , histFF_bkgd_jet_kT ); 

          Double_t angularity_sum_back = *(return_arr_bkgd + 0);
          Double_t angularity_sum_back_nw = *(return_arr_bkgd + 1);
          Double_t bkgd_const_gt_cut = *(return_arr_bkgd + 2);
          Int_t first_pythia_particle_index_in_bkgd = *(return_arr_bkgd + 3);

          histnumtr_bkgd_jet_kT->Fill(bkgd_const_gt_cut); //number of consituents
          N_Trk_tree_bkgd = bkgd_const_gt_cut;

          histangularity_nw_bkgd_jet_kT->Fill( angularity_sum_back_nw / bkgd_constit_pT_vec.size() );

          constituents_bkgd.clear(); 
          histangularity_bkgd_jet_kT->Fill(angularity_sum_back);
          Angularity_tree_bkgd = angularity_sum_back;
          Angularity_NW_tree_bkgd = angularity_sum_back_nw / bkgd_constit_pT_vec.size() ;

          if(bkgd_constit_pT_vec.size() > 0){
          
            histmeanpTtr_bkgd_jet_kT->Fill( std::accumulate(bkgd_constit_pT_vec.begin(), bkgd_constit_pT_vec.end(), 0.0) / bkgd_constit_pT_vec.size() ); //mean track pT in the jet
            Mean_p_T_tree_bkgd = std::accumulate(bkgd_constit_pT_vec.begin(), bkgd_constit_pT_vec.end(), 0.0) / bkgd_constit_pT_vec.size();

            std::sort(bkgd_constit_pT_vec.begin(), bkgd_constit_pT_vec.end()); //sorts least to greatest
            histldtr_bkgd_jet_kT->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size()-1]); // grab leading track in jet

            if(bkgd_constit_pT_vec.size() > 1){
              histsubldtr_bkgd_jet_kT->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 2]); // grab sub-leading track in jet
              histlesub_bkgd_jet_kT->Fill( bkgd_constit_pT_vec[bkgd_constit_pT_vec.size()-1] - bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 2] );
              if(bkgd_constit_pT_vec.size() > 2){
                histthrdldtr_bkgd_jet_kT->Fill(bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 3]); //grad sub-sub leading track in jet
              }
            }

            p_T_1_tree_bkgd = bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 1]; 
            p_T_2_tree_bkgd = bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 2]; 
            p_T_3_tree_bkgd = bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 3]; 
            p_T_4_tree_bkgd = bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 4]; 
            p_T_5_tree_bkgd = bkgd_constit_pT_vec[bkgd_constit_pT_vec.size() - 5]; 
          } //end if check for constits
  
          bkgd_constit_pT_vec.clear();
          X_tru_tree_bkgd = 0; //all backgorund no true
          tree_bkgd_kt->Fill();
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

        event_median_shared = event_median;
        
        cout<<"\n\n\n|||||||||||||The TRUE median rho*area is = "<<event_median<<"||||||||||||||||||"<<endl;
        histmedianrhoebe_bgkd->Fill(event_median/(TMath::Pi()*Jet_Radius*Jet_Radius));
        rho_vec.clear();

      } //end if check for background jets that pass your cuts


/////////////////loop over the anti-kT background only jets///////////////////////////////////////////////

      if(selected_jetsBackground_antikT_sorted.size() > 0){
        for(unsigned j_jet = 0 ; j_jet < selected_jetsBackground_antikT_sorted.size() ; j_jet++ ){ 

          histbackground_antikT_jets->Fill(1.00);
    
          histpT_bkgd_jet_antikT -> Fill(selected_jetsBackground_antikT_sorted[j_jet].pt());
          p_T_tree_antikt_bkgd = selected_jetsBackground_antikT_sorted[j_jet].pt();

          histeta_bkgd_jet_antikT -> Fill(selected_jetsBackground_antikT_sorted[j_jet].eta());
          Eta_tree_antikt_bkgd = selected_jetsBackground_antikT_sorted[j_jet].eta();

          histphi_bkgd_jet_antikT -> Fill(selected_jetsBackground_antikT_sorted[j_jet].phi());
          Phi_tree_antikt_bkgd = selected_jetsBackground_antikT_sorted[j_jet].phi();

          Area_tree_antikt_bkgd = selected_jetsBackground_antikT_sorted[j_jet].area();
          Rho_tree_antikt_bkgd = selected_jetsBackground_antikT_sorted[j_jet].pt()/selected_jetsBackground_antikT_sorted[j_jet].area();

          Double_t jet_pT_antikT_back = selected_jetsBackground_antikT_sorted[j_jet].pt();
          Double_t jet_eta_antikT_back = selected_jetsBackground_antikT_sorted[j_jet].eta();
          Double_t jet_phi_antikT_back = selected_jetsBackground_antikT_sorted[j_jet].phi();

          p_T_corr_tree_antikt_bkgd = jet_pT_antikT_back - event_median_shared; //subtract off the median from kT

          histarea_bkgd_jet_antikT->Fill(selected_jetsBackground_antikT_sorted[j_jet].area());
          histrho_bkgd_jet_antikT->Fill(selected_jetsBackground_antikT_sorted[j_jet].pt()/selected_jetsBackground_antikT_sorted[j_jet].area());

          //rho_vec.push_back(selected_jetsBackground_sorted[i_jet].pt()/selected_jetsBackground_sorted[i_jet].area());
      
          vector <fastjet::PseudoJet> constituents_antikT_bkgd = selected_jetsBackground_antikT_sorted[j_jet].constituents();
          vector <Double_t> antikT_bkgd_constit_pT_vec;

          Int_t nPart_antikT_bkgd = constituents_antikT_bkgd.size();
          for( Int_t i_bkgd_antikT = 0 ; i_bkgd_antikT < nPart_antikT_bkgd ; i_bkgd_antikT++ ){
            if( constituents_antikT_bkgd[i_bkgd_antikT].perp() > 1e-50 ){
              antikT_bkgd_constit_pT_vec.push_back( constituents_antikT_bkgd[i_bkgd_antikT].perp() );
            }
          }

          Double_t *return_arr_bkgd_antikT;

          return_arr_bkgd_antikT = jet_constit_loop( selected_jetsBackground_antikT_sorted[j_jet] , constituents_antikT_bkgd , hist_diagnostic_jet_cone_bkgd_jet_antikT , hist_diagnostic_jet_cone_bkgd_jet_antikT_pT , hist_diagnostic_jet_cone_bkgd_jet_antikT_z , histFF_bkgd_jet_antikT ); 

          Double_t angularity_sum_antikT_back= *(return_arr_bkgd_antikT + 0);
          Double_t angularity_sum_antikT_back_nw = *(return_arr_bkgd_antikT + 1);
          Double_t bkgd_anti_kT_const_gt_cut = *(return_arr_bkgd_antikT + 2);
          Int_t first_pythia_particle_index_in_antikT_bkgd = *(return_arr_bkgd_antikT + 3);
        

          histnumtr_bkgd_jet_antikT->Fill(bkgd_anti_kT_const_gt_cut); //number of consituents
          N_Trk_tree_antikt_bkgd = bkgd_anti_kT_const_gt_cut;

          histangularity_nw_bkgd_jet_antikT->Fill( angularity_sum_antikT_back_nw/antikT_bkgd_constit_pT_vec.size() );
          histangularity_bkgd_jet_antikT->Fill(angularity_sum_antikT_back);

          Angularity_tree_antikt_bkgd = angularity_sum_antikT_back;
          Angularity_NW_tree_antikt_bkgd = angularity_sum_antikT_back_nw/antikT_bkgd_constit_pT_vec.size() ;

          constituents_antikT_bkgd.clear(); 

          if(antikT_bkgd_constit_pT_vec.size() > 0){
          
            histmeanpTtr_bkgd_jet_antikT->Fill( std::accumulate(antikT_bkgd_constit_pT_vec.begin(), antikT_bkgd_constit_pT_vec.end(), 0.0) / antikT_bkgd_constit_pT_vec.size() ); //mean track pT in the jet
           Mean_p_T_tree_antikt_bkgd = std::accumulate(antikT_bkgd_constit_pT_vec.begin(), antikT_bkgd_constit_pT_vec.end(), 0.0) / antikT_bkgd_constit_pT_vec.size();

            std::sort(antikT_bkgd_constit_pT_vec.begin(), antikT_bkgd_constit_pT_vec.end()); //sorts least to greatest
            histldtr_bkgd_jet_antikT->Fill(antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size()-1]); // grab leading track in jet

            if(antikT_bkgd_constit_pT_vec.size() > 1){
              histsubldtr_bkgd_jet_antikT->Fill(antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 2]); // grab sub-leading track in jet
              histlesub_bkgd_jet_antikT->Fill( antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size()-1] - antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 2] );
              if(antikT_bkgd_constit_pT_vec.size() > 2){
                histthrdldtr_bkgd_jet_antikT->Fill(antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 3]); //grad sub-sub leading track in jet
              }
            }

            p_T_1_tree_antikt_bkgd = antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 1]; 
            p_T_2_tree_antikt_bkgd = antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 2]; 
            p_T_3_tree_antikt_bkgd = antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 3]; 
            p_T_4_tree_antikt_bkgd = antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 4]; 
            p_T_5_tree_antikt_bkgd = antikT_bkgd_constit_pT_vec[antikT_bkgd_constit_pT_vec.size() - 5]; 
          } //end if check for constits
  
          antikT_bkgd_constit_pT_vec.clear();
          X_tru_tree_antikt_bkgd = 0; //all backgorund no true
          tree_bkgd_antikt->Fill();
        } //end looping over the background antikT jets

      } //end if check for antikT background jets that pass your cuts

/////////////////////////////////////////////////////////////////////////////////////////////////////////
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
      cout<<"\n\n\n|||||||||||||The estimated median rho*area (drop 2) is = "<<event_median_realistic<<"||||||||||||||||||"<<endl;
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
          pbgOut<<"p_T, Eta, Phi, Area, Rho, p_T-corr, N-Trk, Angularity, Angularity-NW, Mean-p_T, p_T_1, p_T_2, p_T_3, p_T_4, p_T_5, distmatch, XMatch, X_tru"<<endl;
          HEADER__--;
        }

        for(unsigned t_jet = 0 ; t_jet < selected_jetsTOTAL_sorted.size() ; t_jet++ ){ 

          histpythia_AND_bkgd_jets->Fill(1.00);
    
          histpT_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].pt());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].pt()<<", ";
          p_T_tree = selected_jetsTOTAL_sorted[t_jet].pt();
          p_T_tree_pat = selected_jetsTOTAL_sorted[t_jet].pt();

          histeta_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].eta());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].eta()<<", ";
          Eta_tree = selected_jetsTOTAL_sorted[t_jet].eta();
          Eta_tree_pat = selected_jetsTOTAL_sorted[t_jet].eta();

          histphi_pytha_AND_bkgd_jet -> Fill(selected_jetsTOTAL_sorted[t_jet].phi());

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].phi()<<", ";
          Phi_tree = selected_jetsTOTAL_sorted[t_jet].phi();
          Phi_tree_pat = selected_jetsTOTAL_sorted[t_jet].phi();

          pbgOut<<selected_jetsTOTAL_sorted[t_jet].area()<<", ";
          Area_tree = selected_jetsTOTAL_sorted[t_jet].area() ;
          Area_tree_pat = selected_jetsTOTAL_sorted[t_jet].area() ;

          pbgOut<<(selected_jetsTOTAL_sorted[t_jet].pt() / selected_jetsTOTAL_sorted[t_jet].area())<<", ";
          Rho_tree = selected_jetsTOTAL_sorted[t_jet].pt() / selected_jetsTOTAL_sorted[t_jet].area() ;
          Rho_tree_pat = selected_jetsTOTAL_sorted[t_jet].pt() / selected_jetsTOTAL_sorted[t_jet].area() ;

          if( selected_jetsTOTAL_sorted[t_jet].area() > 0.6*TMath::Pi()*Jet_Radius*Jet_Radius ){
            histpT_area_corr_pytha_AND_bkgd_jet->Fill(selected_jetsTOTAL_sorted[t_jet].pt() - event_median_realistic );
          }
          pbgOut<<selected_jetsTOTAL_sorted[t_jet].pt()-event_median_realistic<<", ";
          p_T_corr_tree = selected_jetsTOTAL_sorted[t_jet].pt()-event_median_realistic;
          p_T_corr_tree_pat = selected_jetsTOTAL_sorted[t_jet].pt()-event_median_realistic;

          histarea_pytha_AND_bkgd_jet->Fill(selected_jetsTOTAL_sorted[t_jet].area());
          histrho_pytha_AND_bkgd_jet->Fill( selected_jetsTOTAL_sorted[t_jet].pt() / selected_jetsTOTAL_sorted[t_jet].area() );
      
          vector <fastjet::PseudoJet> constituents_total = selected_jetsTOTAL_sorted[t_jet].constituents();
          vector <Double_t> total_constit_pT_vec;

          Int_t nPart_total = constituents_total.size();
          for( Int_t i_bkgd_total = 0 ; i_bkgd_total < nPart_total ; i_bkgd_total++ ){
            if( constituents_total[i_bkgd_total].perp() > 1e-50 ){
              total_constit_pT_vec.push_back( constituents_total[i_bkgd_total].perp() );
            }
          }

          Double_t *return_arr_tot;
          return_arr_tot = jet_constit_loop( selected_jetsTOTAL_sorted[t_jet] , constituents_total , hist_diagnostic_jet_cone_pythia_AND_bkgd_jet , hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT , hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z , histFF_pytha_AND_bkgd_jet );       

          Double_t angularity_sum= *(return_arr_tot + 0);
          Double_t angularity_sum_nw = *(return_arr_tot + 1);
          Double_t tot_const_gt_cut = *(return_arr_tot + 2);
          Int_t first_pythia_particle_index_in_tot = *(return_arr_tot + 3);
          Double_t tru_sum = *(return_arr_tot + 4);
          Double_t fake_sum = *(return_arr_tot + 5);

          pbgOut<<total_constit_pT_vec.size()<<", ";
          N_Trk_tree = total_constit_pT_vec.size();
          N_Trk_tree_pat = total_constit_pT_vec.size();

          histangularity_nw_pytha_AND_bkgd_jet->Fill( angularity_sum_nw / total_constit_pT_vec.size() );

          histnumtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec.size()); //number of constituents

          Double_t X_tru = tru_sum/total_constit_pT_vec.size();
          Double_t X_fake = fake_sum/total_constit_pT_vec.size();

          //cout<<"\ntrusum = "<<tru_sum<<" jet size (no ghosts) = "<<total_constit_pT_vec.size()<<" Xtru = "<<X_tru<<endl;
          //cout<<"fakesum = "<<fake_sum<<" jet size (no ghosts) = "<<total_constit_pT_vec.size()<<" Xfake = "<<X_fake<<endl;

          histX_pythia_tru_pythia_AND_bkgd->Fill((X_tru*total_constit_pT_vec.size())/selected_jetsTOTAL_sorted[t_jet].pt());

          //___________________DOING PAT AND ANTONIOS MATCHING METIRC_______________________________________________________//
          Bool_t found_match = kFALSE;
          Double_t delr_match = 1000000; // 1000000 default for no match (way too big to be realistic !!)
          Double_t match_momentum_fraction = 0; //no match is the default !!!!!
          for(Int_t i_pj = 0 ; i_pj < selected_jetsPythia_sorted_part.size() ; i_pj++ ){ //you have the index now loop through all the pythia jets
            vector <fastjet::PseudoJet> constituents_pyth = selected_jetsPythia_sorted_part[i_pj].constituents();
            for(Int_t i_pc = 0 ; i_pc < constituents_pyth.size() ; i_pc++ ){
              if( constituents_pyth[i_pc].user_index() == first_pythia_particle_index_in_tot  ){
                Double_t delphi_match = selected_jetsPythia_sorted_part[i_pj].delta_phi_to(selected_jetsTOTAL_sorted[t_jet]) ;
                Double_t deleta_match = selected_jetsTOTAL_sorted[t_jet].eta() - selected_jetsPythia_sorted_part[i_pj].eta() ;
                delr_match = TMath::Sqrt( (delphi_match*delphi_match) + (deleta_match*deleta_match) );
                match_momentum_fraction = tru_sum / selected_jetsPythia_sorted_part[i_pj].pt() ;
                found_match = kTRUE;
                break; //break out of inner loop
              }
            }
            if( found_match ){
              break; //break out of outer loop
            }
          }
         geom_match_tree_pat = delr_match;
         mom_frac_match_tree_pat = match_momentum_fraction;
         //________________________________________________________________________________________________________________//
 
          Int_t label = X_tru<0.1?0:X_tru<0.5?1:X_tru<0.9?2:3;//0-fake, 1-mostly fake, 2-mostly real, 3-real

          histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->Fill(X_tru, selected_jetsTOTAL_sorted[t_jet].pt());

          histX_bkgd_fake_pythia_AND_bkgd->Fill((X_fake*total_constit_pT_vec.size())/selected_jetsTOTAL_sorted[t_jet].pt());

          histangularity_pytha_AND_bkgd_jet->Fill(angularity_sum);

          pbgOut<<angularity_sum<<", ";
          pbgOut<<(angularity_sum_nw / total_constit_pT_vec.size())<<", ";
          Angularity_tree = angularity_sum;
          Angularity_NW_tree = angularity_sum_nw / total_constit_pT_vec.size() ;

          Angularity_tree_pat = angularity_sum;
          Angularity_NW_tree_pat = angularity_sum_nw / total_constit_pT_vec.size() ;

          if(total_constit_pT_vec.size() > 0){
          
            histmeanpTtr_pytha_AND_bkgd_jet->Fill( std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size() ); //mean track pT in the jet

            pbgOut<<std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size()<<", ";
            Mean_p_T_tree = std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size();
            Mean_p_T_tree_pat = std::accumulate(total_constit_pT_vec.begin(), total_constit_pT_vec.end(), 0.0) /  total_constit_pT_vec.size();

            std::sort(total_constit_pT_vec.begin(), total_constit_pT_vec.end()); //sorts least to greatest
            histldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 1]); // grab leading track in jet

            pbgOut<<total_constit_pT_vec[total_constit_pT_vec.size() - 1]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 2]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 3]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 4]<<", "<<total_constit_pT_vec[total_constit_pT_vec.size() - 5]<<", ";
            p_T_1_tree = total_constit_pT_vec[total_constit_pT_vec.size() - 1];
            p_T_2_tree = total_constit_pT_vec[total_constit_pT_vec.size() - 2];
            p_T_3_tree = total_constit_pT_vec[total_constit_pT_vec.size() - 3];
            p_T_4_tree = total_constit_pT_vec[total_constit_pT_vec.size() - 4];
            p_T_5_tree = total_constit_pT_vec[total_constit_pT_vec.size() - 5];


            p_T_1_tree_pat = total_constit_pT_vec[total_constit_pT_vec.size() - 1];
            p_T_2_tree_pat = total_constit_pT_vec[total_constit_pT_vec.size() - 2];
            p_T_3_tree_pat = total_constit_pT_vec[total_constit_pT_vec.size() - 3];
            p_T_4_tree_pat = total_constit_pT_vec[total_constit_pT_vec.size() - 4];
            p_T_5_tree_pat = total_constit_pT_vec[total_constit_pT_vec.size() - 5];

            if(total_constit_pT_vec.size() > 1){
              histsubldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 2]); // grab sub-leading track in jet
              histlesub_pytha_AND_bkgd_jet->Fill( total_constit_pT_vec[total_constit_pT_vec.size() - 1] - total_constit_pT_vec[total_constit_pT_vec.size() - 2] );
              if(total_constit_pT_vec.size() > 2){
                histthrdldtr_pytha_AND_bkgd_jet->Fill(total_constit_pT_vec[total_constit_pT_vec.size() - 3]); //grab sub-sub leading track in jet
              }
            }
          } //end if check for total jet constit
          pbgOut<<delr_match<<", ";
          pbgOut<<match_momentum_fraction<<", ";
          pbgOut<<X_tru<<endl;
          X_tru_tree = (X_tru*total_constit_pT_vec.size())/selected_jetsTOTAL_sorted[t_jet].pt();
          X_tru_tree_pat = (X_tru*total_constit_pT_vec.size())/selected_jetsTOTAL_sorted[t_jet].pt();
          tree->Fill();
          tree_pat->Fill();
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
  Run_Time_Stats->cd();

  Int_t Sum_Time = 0;
  Int_t Total_Time = 0;

  Int_t *time_interval_array;
  Int_t *event_array;

  time_interval_array = new Int_t[time_array.size()-1];
  event_array = new Int_t[time_array.size()-1];

  for( Int_t i = 0 ; i < time_array.size() - 1; i++){
    time_interval_array[i] = time_array[i + 1] - time_array[i];
    Sum_Time = Sum_Time + time_interval_array[i];
    event_array[i] = i + 1;
  }

  if(RT_Stats){
    TCanvas *c1 = new TCanvas("c1" , "Run Time For 1 Background Event vs Event Number (Events Run in Sequence)", 200, 10,500,300); 

    TH1D *histruntime = new TH1D("histruntime", "Run Times for Background Generator",8,((Sum_Time/(time_array.size() - 1)) - 1),((Sum_Time/(time_array.size() - 1)) + 1)); //for non-random distribution
    histruntime -> Sumw2();
    histruntime->SetXTitle("Run Time For 1 Event (Seconds)");
    histruntime->SetYTitle("Counts");

    TGraph *rtgraph = new TGraph(time_array.size()-1 , event_array , time_interval_array);
    rtgraph->SetTitle("Run Time For 1 Background Event vs Event Number (Events Run in Sequence)");
    rtgraph->GetXaxis()->SetTitle("Event Iteration Number");
    rtgraph->GetYaxis()->SetTitle("Run Time for 1 Event (Seconds)");
    rtgraph->Write();
    c1->cd(1); 
    rtgraph->Draw("AC");
    c1->Draw();
    c1->Write("");
    for( Int_t j = 0 ; j < (time_array.size() - 1); j++){
      histruntime->Fill(time_interval_array[j]);
    }
    histruntime->Write();
  }


//_______________________________________PYTHIA ONLY STUFF FIRST______________________________________________________

  pythia_only->cd(); 

  //tree_pyth->Write();

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
  histarea_pyth_jet->Write();
  histrho_pyth_jet->Write();
  histangularity_nw_pyth_jet->Write();

  hist_diagnostic_jet_cone_pyth_jet->Write();
  hist_diagnostic_jet_cone_pyth_jet_pT->Write();
  hist_diagnostic_jet_cone_pyth_jet_z->Write();
  

//______________________________________END PYTHIA ONLY STUFF______________________________________________________________//

//______________________________________START OF BACKGROUND ONLY STUFF______________________________________________________________//

  background_only->cd();

  //tree_bkgd_kt->Write();
  //tree_bkgd_antikt->Write();

  histbackground_events->Write();
  histbackground_kT_jets->Write();
  histbackground_antikT_jets->Write();

  histpT_bkgd_part->Write();
  histeta_bkgd_part->Write();
  histphi_bkgd_part->Write();

  histpT_bkgd_jet_kT->Write();
  histpT_bkgd_jet_antikT->Write();
  histeta_bkgd_jet_kT->Write();
  histeta_bkgd_jet_antikT->Write();
  histphi_bkgd_jet_kT->Write();
  histphi_bkgd_jet_antikT->Write();

  histnumtr_bkgd_jet_kT->Write();
  histnumtr_bkgd_jet_antikT->Write();
  histmeanpTtr_bkgd_jet_kT->Write();
  histmeanpTtr_bkgd_jet_antikT->Write();
  histangularity_bkgd_jet_kT->Write();
  histangularity_bkgd_jet_antikT->Write();

  histldtr_bkgd_jet_kT->Write();
  histldtr_bkgd_jet_antikT->Write();
  histsubldtr_bkgd_jet_kT->Write();
  histsubldtr_bkgd_jet_antikT->Write();
  histthrdldtr_bkgd_jet_kT->Write();
  histthrdldtr_bkgd_jet_antikT->Write();

  histFF_bkgd_jet_kT->Write();
  histFF_bkgd_jet_antikT->Write();
  histlesub_bkgd_jet_kT->Write();
  histlesub_bkgd_jet_antikT->Write();
  histarea_bkgd_jet_kT->Write();
  histarea_bkgd_jet_antikT->Write();
  histrho_bkgd_jet_kT->Write();
  histrho_bkgd_jet_antikT->Write();
  histangularity_nw_bkgd_jet_kT->Write();
  histangularity_nw_bkgd_jet_antikT->Write();

  hist_diagnostic_jet_cone_bkgd_jet_kT->Write();
  hist_diagnostic_jet_cone_bkgd_jet_kT_pT->Write();
  hist_diagnostic_jet_cone_bkgd_jet_kT_z->Write();
  hist_diagnostic_jet_cone_bkgd_jet_antikT->Write();
  hist_diagnostic_jet_cone_bkgd_jet_antikT_pT->Write();
  hist_diagnostic_jet_cone_bkgd_jet_antikT_z->Write();

  histPsi_1->Write();
  histPsi_2->Write();
  histPsi_3->Write();
  histPsi_4->Write();
  histPsi_5->Write();

  histmedianrhoebe_bgkd->Write();


//______________________________________END BACKGROUND ONLY STUFF______________________________________________________________// 

//_______________________________________START OF BACKGROUND AND PYTHIA_________________________________________________________//

  background_and_pythia->cd();

  //tree->Write();

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
  histarea_pytha_AND_bkgd_jet->Write();
  histrho_pytha_AND_bkgd_jet->Write();
  histangularity_nw_pytha_AND_bkgd_jet->Write();

  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet->Write();
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_pT->Write();
  hist_diagnostic_jet_cone_pythia_AND_bkgd_jet_z->Write();

  histmedianrhoebe_pythia_AND_bkgd->Write();

  histX_pythia_tru_pythia_AND_bkgd->Write();
  histX_bkgd_fake_pythia_AND_bkgd->Write();

  histX_pythia_tru_pythia_AND_bkgd_jet_pT_raw->Write();



//_______________________________________END OF BACKGROUND AND PYTHIA_________________________________________________________//

  ff->Write();
  ff->Close();

} //end function


