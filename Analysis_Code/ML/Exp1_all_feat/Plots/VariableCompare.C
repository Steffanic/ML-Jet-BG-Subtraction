void SetHistoStyle(TH1F *histo,int marker, int color,char *xtitle,char *ytitle){
  histo->SetMarkerColor(color);
  histo->SetLineColor(color);
  histo->SetMarkerStyle(marker);
  histo->SetMarkerSize(2.0);
  histo->GetXaxis()->SetTitle(xtitle);
  histo->GetYaxis()->SetTitle(ytitle);
  histo->GetXaxis()->SetTitleSize(0.08);
  histo->GetYaxis()->SetTitleSize(0.08);
  histo->GetXaxis()->SetTitleOffset(1.0);
  histo->GetYaxis()->SetTitleOffset(1.0);
  histo->GetXaxis()->SetLabelSize(0.08);
  histo->GetYaxis()->SetLabelSize(0.08);
}
int pythiaMarker = 29;
int pythiaColor = kRed;
int tennGenMarker = 20;
int tennGenColor = kBlue;
int bothMarker = 21;
int bothColor = 1;
void SetPythiaStyle(TH1F *histo,char *xtitle,char *ytitle){
  SetHistoStyle(histo,pythiaMarker,pythiaColor,xtitle,ytitle);
}
void SetTennGenStyle(TH1F *histo,char *xtitle,char *ytitle){
  SetHistoStyle(histo,tennGenMarker,tennGenColor,xtitle,ytitle);
}
void SetBothStyle(TH1F *histo,char *xtitle,char *ytitle){
  SetHistoStyle(histo,bothMarker,bothColor,xtitle,ytitle);
}
void Normalize(TH1F *histo){
  Int_t entries = histo->GetEntries();
  Float_t binwidth = histo->GetBinWidth(1);//assumes constant bin width
  histo->Scale(1.0/((Float_t)entries)/binwidth);
}
TCanvas *GetCanvas(char *name, Bool_t logy=kTRUE){
   TCanvas *cArea = new TCanvas(name,name,600,400);
   cArea->SetFillColor(0);
   cArea->SetBorderMode(0);
   cArea->SetBorderSize(0);
   if(logy)cArea->SetLogy();
   cArea->SetLeftMargin(0.166107);
   cArea->SetRightMargin(0.03);
   cArea->SetTopMargin(0.026738);
   cArea->SetBottomMargin(0.17);
   cArea->SetFrameBorderMode(0);
   return cArea;
}
TTree *treePythia;// = (TTree*) file->Get("Pythia");
TTree *treeTennGen;// = (TTree*) file->Get("TennGen-antikT-HF-0-CB-0");
TTree *treeBoth;// = (TTree*) file->Get("Pythia-and-TennGen-HF-0-CB-0");
TCanvas *GetVarHisto(char *varname,int low = 0, int high = 1.0,int nbins = 50, Bool_t logy=kTRUE){
  TCanvas *canvas = GetCanvas(varname,logy);
  char name1[200];
  char name2[200];

  sprintf(name1,"h%sPythia",varname);
  sprintf(name2,"%s in Pythia",varname);
  TH1F *hPythia = new TH1F(name1,name2,nbins,low,high);
  sprintf(name2,"%s>>%s",varname,name1);
  canvas->cd();
  treePythia->Draw(name2,"","e");
  sprintf(name1,"dN_{jet}/d%s",varname);
  SetPythiaStyle(hPythia,varname,name1);
  Normalize(hPythia);

  sprintf(name1,"h%sTennGen",varname);
  sprintf(name2,"%s in TennGen",varname);
  TH1F *hTennGen = new TH1F(name1,name2,nbins,low,high);
  sprintf(name2,"%s>>%s",varname,name1);
  treeTennGen->Draw(name2,"","e");
  sprintf(name1,"dN_{jet}/d%s",varname);
  SetTennGenStyle(hTennGen,varname,name1);
  Normalize(hTennGen);

  sprintf(name1,"h%sBoth",varname);
  sprintf(name2,"%s in Both",varname);
  TH1F *hBoth = new TH1F(name1,name2,nbins,low,high);
  sprintf(name2,"%s>>%s",varname,name1);
  treeBoth->Draw(name2,"","e");
  sprintf(name1,"dN_{jet}/d%s",varname);
  SetBothStyle(hBoth,varname,name1);
  Normalize(hBoth);
  if(hTennGen->GetMaximum()> hPythia->GetMaximum()) hPythia->SetMaximum(hTennGen->GetMaximum());
  if(hBoth->GetMaximum()> hPythia->GetMaximum()) hPythia->SetMaximum(hBoth->GetMaximum());
  hPythia->Draw("e");
  hTennGen->Draw("e same");
  hBoth->Draw("e same");
  sprintf(name1,"%s.png",varname);
  TLegend *leg3 = new TLegend(0.666107,0.737968,0.823826,0.933155);
   leg3->SetBorderSize(0);
leg3->SetTextFont(72);
 leg3->SetTextSize(0.0481283);
 leg3->SetLineColor(0);
 leg3->SetLineStyle(0);
 leg3->SetLineWidth(0);
 leg3->SetFillColor(0);
 leg3->SetFillStyle(1);
 leg3->AddEntry(hPythia,"Signal","p");
 leg3->AddEntry(hTennGen,"Background","p");
 leg3->AddEntry(hBoth,"Signal+Background","p");
 leg3->Draw();
  canvas->SaveAs(name1);
  return canvas;
}
void CompareTwo(char *var1, char *var2,int low1 = 0, int high1 = 1.0,int low2 = 0, int high2 = 1.0,int nbins1 = 50,int nbins2 = 50,Bool_t logz=kTRUE){

  char name1[200];
  char name2[200];
  sprintf(name2,"%s vs %s",var1,var2);
   TCanvas *cArea = new TCanvas(name2,name2,800,800);
   cArea->SetFillColor(0);
   cArea->SetBorderMode(0);
   cArea->SetBorderSize(0);
   if(logz)cArea->SetLogz();
   cArea->SetLeftMargin(0.166107);
   cArea->SetRightMargin(0.03);
   cArea->SetTopMargin(0.026738);
   cArea->SetBottomMargin(0.17);
   cArea->SetFrameBorderMode(0);
   cArea->Divide(2,2);
  sprintf(name2,"h%svs%sPythia",var1,var2);
  sprintf(name1,"%s:%s>>%s",var1,var2,name2);
  TH2F *hPythia = new TH2F(name2,name2,nbins1,low1,high1,nbins2,low2,high2);
  cArea->cd(1);
  cout<<name1<<endl;
  treePythia->Draw(name1,"","colz");
  hPythia->GetXaxis()->SetTitle(var1);
  hPythia->GetYaxis()->SetTitle(var2);

  sprintf(name2,"h%svs%sTennGen",var1,var2);
  sprintf(name1,"%s:%s>>%s",var1,var2,name2);
  cout<<name1<<endl;
  TH2F *hTennGen = new TH2F(name2,name2,nbins1,low1,high1,nbins2,low2,high2);
  cArea->cd(2);
  treeTennGen->Draw(name1,"","colz");
  hTennGen->GetXaxis()->SetTitle(var1);
  hTennGen->GetYaxis()->SetTitle(var2);

  sprintf(name2,"h%svs%sBoth",var1,var2);
  sprintf(name1,"%s:%s>>%s",var1,var2,name2);
  TH2F *hBoth = new TH2F(name2,name2,nbins1,low1,high1,nbins2,low2,high2);
  cArea->cd(3);
  cout<<name1<<endl;
  treeBoth->Draw(name1,"","colz");
  hBoth->GetXaxis()->SetTitle(var1);
  hBoth->GetYaxis()->SetTitle(var2);

  sprintf(name2,"h%svs%sRatio",var1,var2);
  TH2F *hRatio = (TH2F*) hPythia->Clone(name2);
  hRatio->Divide(hTennGen);
  cArea->cd(4);
  hRatio->Draw("colz");
  sprintf(name2,"%svs%s.png",var1,var2);
  cArea->SaveAs(name2);


}
void VariableCompare(char *filename = "merged_ML_output_LOW_STATS3.root"){
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
  TFile *file = new TFile(filename);
  treePythia = (TTree*) file->Get("Pythia");
  treeTennGen = (TTree*) file->Get("TennGen-antikT-HF-0-CB-0");
  treeBoth = (TTree*) file->Get("Pythia-and-TennGen-HF-0-CB-0");
  //Vars:
  //pT,Eta,Phi,Area,Rho,pTcorr,NTrk,Angularity,Angularity_NW,MeanpT,pT1,pT2,PT3,pT4,pT5,XTru
  //This line draws a 2D histogram of eta vs phi for all jets with area > 0
  //treePythia->Draw("Phi:Eta","Area>0","surf");

  //TCanvas *cArea = GetCanvas("Area");
  TCanvas *cAngularity = GetVarHisto("Angularity");
  TCanvas *cAngularityNW = GetVarHisto("Angularity_NW");
  TCanvas *cArea = GetVarHisto("Area");
  TCanvas *cmeanpT = GetVarHisto("MeanpT",0.2,2.0);
  TCanvas *cpT1 = GetVarHisto("pT1",0,20);
  TCanvas *cpT2 = GetVarHisto("pT2",0,20);
  CompareTwo("Angularity","Angularity_NW",0,0.5,0,0.5);
  CompareTwo("Area","MeanpT",0,1.0,0,2.0);
  CompareTwo("Area","pT",0,1.0,0,30.0);
  CompareTwo("pT1","pT2",0,10,0,10);


}
