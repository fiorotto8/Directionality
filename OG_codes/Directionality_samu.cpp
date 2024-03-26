#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TArrayF.h>
#include <TH1.h>
#include <TMarker.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TLine.h>
#include <TMath.h>



void ScIndicesElem(int nSc,int* ScNelements, std::vector<int>* B, std::vector<int> *E){
  B->clear();
  E->clear();
  
  int parcount=0;
 
  for(int i=0;i<nSc;i++){
    B->push_back(parcount);
    E->push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}


int main(int argc, char** argv){

  //////////////////////////Remove message of TCanvas saved ///////////////////////////
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0000);
  ////////////////////////////////////Get File //////////////////////////////////////
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");

  ///////////////////////////////////Set Branches and Define Variables////////////////////////////////////

  Int_t NPIP;
  Float_t wFac;
  
  std::ifstream inputFile(argv[2]);
  std::string line;   
  
  while (getline(inputFile, line)){
    std::istringstream ss(line);
    std::string name;
    ss >> NPIP >> wFac;
  }
  
  //////////////////////////////////GetDirParameters////////////////////////////////
  
  int nmax=350000;
  int nscmax=40;
  int npixel=2304;
      
  UInt_t nSc;
  int run;
  int event;
             //Pixels
  int ScNpixels[nscmax];
  Float_t XPix[nmax];
  Float_t YPix[nmax];
  Float_t ZPix[nmax];
  int ScNpixelsall[nscmax];
  float width[nscmax];
  float length[nscmax];
  float skew[nscmax];
  float rms[nscmax];
  float size[nscmax];
  float sc_integral[nscmax];
  float sc_integralCal[nscmax];

  
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);           
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_nintpixels",ScNpixels);
  tree->SetBranchAddress("sc_nallintpixels",ScNpixelsall);
  tree->SetBranchAddress("sc_ypixelcoord",XPix);
  tree->SetBranchAddress("sc_xpixelcoord",YPix);
  tree->SetBranchAddress("sc_zpixel",ZPix);
  tree->SetBranchAddress("sc_width",width);
  tree->SetBranchAddress("sc_length",length);
  tree->SetBranchAddress("sc_size",size);
  tree->SetBranchAddress("sc_integral",sc_integral);

  /////////////////////////////////OutTree////////////////////////////////////////////////////
  TFile* fout = new TFile("Tracks/NtuOut.root","recreate");

  TTree* myTree = new TTree("myTree","");
  myTree->Branch("nSc",&nSc,"nSc/I");
  myTree->Branch("run",&run,"run/I");
  myTree->Branch("event",&event,"event/I");
  myTree->Branch("nSc",&nSc,"nSc/I");
  myTree->Branch("ScNPixel",ScNpixels,"ScNpixels[nSc]/I");
  myTree->Branch("ScNAllPixel",ScNpixelsall,"ScNAllPixel[nSc]/I");
  myTree->Branch("ScWidth",width,"ScWidth[nSc]/F");
  myTree->Branch("ScLength",length,"ScLength[nSc]/F");
  myTree->Branch("ScSize",size,"ScSize[nSc]/F");
  myTree->Branch("ScIntegral",sc_integral,"ScIntegral[nSc]/F");
  myTree->Branch("ScIntegralCal",sc_integralCal,"ScIntegralCal[nSc]/F");  
  myTree->Branch("ScRMSOnMA",rms,"ScRMSOnMA[nSc]/F");
  myTree->Branch("ScSkewOnMA",skew,"ScSkewOnMA[nSc]/F");
  

  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  std::vector<int> BeginScPix;
  std::vector<int> EndScPix;

  int counter=0;
  
  TH1F* HistoDirection = new TH1F("HistoDirection","HistoDirection",100,-3.16,3.16);
    
  for(int k=0;k<tree->GetEntries();k++){

    std::cout << "Event:\t" << k << "/" << tree->GetEntries() << std::endl; 
        
    counter=0;
    tree->GetEntry(k);
    
    ScIndicesElem(nSc,ScNpixels,&BeginScPix,&EndScPix);
    
    for(int i=0;i<nSc;i++){
                  
      if(ScNpixels[i]>14000){continue;}
      if(nSc>1){continue;}
            
      Analyzer* Traccia = new Analyzer(Form("Track%i_run%i_evt%i",event,run,counter),XPix,YPix,ZPix,BeginScPix[i],EndScPix[i]);

      Traccia->SetWScal(wFac);
      Traccia->SetNPIP(NPIP);
      Traccia->ApplyThr();
      Traccia->RemoveNoise();
      Traccia->ImpactPoint(Form("TrackIPRegion%i_run%i_evt%i",event,run,counter));
      Traccia->ScaledTrack(Form("TrackScaled%i_run%i_evt%i",event,run,counter));
      Traccia->Direction();
      Traccia->ImprCorrectAngle();
      Traccia->BuildLineDirection();

      Traccia->SavePicDir(Form("Track%i.png",k));
      if(k%10==0)Traccia->SaveRootFile(Form("Track%i.root",k));
      
      counter++;
      
      if(nSc==1) HistoDirection->Fill(Traccia->GetDir());

    }//Chiudo for i on Superclusters
    
    //myTree->Fill();
        
  }//chiudo for k TreeEntries

  TF1* fitGaus = new TF1("fitGaus","gaus");
  HistoDirection->Fit(fitGaus);

  gStyle->SetOptStat(0000);
  gStyle->SetOptFit(1111);
  TCanvas* c = new TCanvas();
  HistoDirection->Draw();
  c->SaveAs("Resolution.png");
  
  myTree->Write();
  fout->Save();
  fout->Close();


    
  return 0;
  
}
