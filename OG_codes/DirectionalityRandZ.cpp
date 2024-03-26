//USE `root-config --cflags --libs` TO COMPILE
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

#include <TTree.h>
#include <TFile.h>
#include <TArrayF.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TSystem.h>
#include <TArc.h>
#include <TMath.h>

std::vector<int> SortSlices(Int_t nSc,Float_t* nSl, Float_t* X, Float_t* Y, Float_t* E);
Int_t GetIndexMinChiOnLine(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,TGraph* plot,TF1* line,Int_t parcount, Int_t radius);
Float_t DistPit(Float_t x1,Float_t y1,Float_t x2,Float_t y2);
Int_t GetIndexMaxDist(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount);
Int_t HasNeighbourTemp(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t radius);
Int_t HasNeighbour(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount,Int_t radius);
Int_t GetIndexMinDistTemp(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl);
Int_t GetIndexMinDist(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount);
void MakeCanvas(TH2F* Plot,Float_t XBar,Float_t YBar,Float_t Phi,Int_t i,Int_t k,Float_t RMS,Float_t rmin, Float_t rmax,Float_t XIP,Float_t YIP,Float_t AngleDir,TF1* DirGraph);
void MakeCanvasDist(TH2F* Plot,Float_t XBar,Float_t YBar,Float_t Phi,Int_t i,Int_t k,Float_t RMS,Float_t rmin, Float_t rmax,Float_t XIP,Float_t YIP,Float_t AngleDir, Float_t RSel,Float_t XIPPrev,Float_t YIPPrev,Float_t ThTrue);

//ScIndicesElem gives back B, E arrays that contain the beggining and end of every pixel
void ScIndicesElem(Int_t nSc,Int_t* ScNelements, std::vector<Int_t>* B, std::vector<Int_t> *E){
  B->clear();
  E->clear();

  Int_t parcount=0;

  for(int i=0;i<nSc;i++){
    B->push_back(parcount);
    E->push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}

void ScIndicesElem(Int_t nSc,Float_t* ScNelements, std::vector<Int_t>* B, std::vector<Int_t> *E){
  B->clear();
  E->clear();

  Int_t parcount=0;

  for(int i=0;i<nSc;i++){
    B->push_back(parcount);
    E->push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}

void Barycenter(TH2F* Track,Float_t* XBar, Float_t* YBar){

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();

  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();

  Float_t z;
  
  Float_t Xb=0;
  Float_t Yb=0;
  Float_t ChargeTot=0;

  for(int i=XBinMin; i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      z=Track->GetBinContent(i,j);
      if(z>0){
	Xb+=(z*Track->GetXaxis()->GetBinCenter(i));
	Yb+=(z*Track->GetYaxis()->GetBinCenter(j));
	ChargeTot+=z;
      }
    }
  }
    
  Xb/=ChargeTot;
  Yb/=ChargeTot;

  *XBar=Xb;
  *YBar=Yb;
}

Float_t RMSOnLine(TH2F* Track,Float_t XBar, Float_t YBar,Float_t Phi){

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();

  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();


  Float_t RMS=0;
  Float_t ChargeTot=0;

  Float_t X,Y,Q;
  
  for(int i= XBinMin;i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      X=Track->GetXaxis()->GetBinCenter(i);
      Y=Track->GetYaxis()->GetBinCenter(j);
      Q=Track->GetBinContent(i,j);
      if(Q>0){
	RMS+=Q*( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi) )*( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi) );
	ChargeTot+=Q;
      }//chiudo if
    }//chiudo for j
  }//chiudo for i
  
  return RMS/=ChargeTot;
  
}

Float_t SkewOnLine(TH2F* Track,Float_t XBar, Float_t YBar,Float_t Phi){

  Float_t Skew=0;
  Float_t ChargeTot=0;

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();

  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();

  
  Float_t X,Y,Q;
  
  for(int i= XBinMin;i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      X=Track->GetXaxis()->GetBinCenter(i);
      Y=Track->GetYaxis()->GetBinCenter(j);
      Q=Track->GetBinContent(i,j);
      if(Q>0){
	Skew+=Q*( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi) )*( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi))*( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi));
	ChargeTot+=Q;
      }
    }
  }
  
  return Skew/=ChargeTot;
  
}

Float_t AngleLineMaxRMS(TH2F* Track,Float_t XBar, Float_t YBar,Float_t* RMSOnLineVal,Float_t* RMSOnLinePerpVal){

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();
  
  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();
  
  Float_t Sum1=0;
  Float_t Sum2=0;
  Float_t Phi;
  Float_t RmsAng;
  Float_t RmsAngPerp;

  Float_t X,Y,Q; 
 
  for(int i=XBinMin; i< XBinMax; i++){
    for(int j=YBinMin; j< YBinMax; j++){
      X=Track->GetXaxis()->GetBinCenter(i);
      Y=Track->GetYaxis()->GetBinCenter(j);
      Q=Track->GetBinContent(i,j);
      if(Q>0){
	Sum1+=Q*( X-XBar )*( Y-YBar );
	Sum2+=Q*( (Y-YBar)*(Y-YBar) - (X-XBar)*(X-XBar) );
      }
    }
  }
  
  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);
  
  RmsAng=RMSOnLine(Track,XBar,YBar,Phi);
  RmsAngPerp=RMSOnLine(Track,XBar,YBar,Phi+TMath::Pi()/2);
  
  if( RmsAng > RmsAngPerp ){
    *RMSOnLineVal=RmsAng;
    *RMSOnLinePerpVal=RmsAngPerp;
    return Phi;
  } else {
    *RMSOnLineVal=RmsAngPerp;
    *RMSOnLinePerpVal=RmsAng;
    if(Phi+TMath::Pi()/2>TMath::Pi()/2){
      return Phi+TMath::Pi()/2-TMath::Pi();
    } else{
      return Phi+TMath::Pi()/2;
    }
  }
  
}

Float_t GetPointSkew(Float_t X,Float_t Y,Float_t XBar,Float_t YBar,Float_t Phi,Float_t Skew){
  
  return ( (X-XBar)*cos(Phi) + (Y-YBar)*sin(Phi) )/Skew;
  
}

Float_t DistCm(Float_t X,Float_t Y,Float_t XBar,Float_t YBar,Float_t RMSMax){

  return sqrt( ( (X-XBar)*(X-XBar) + (Y-YBar)*(Y-YBar) ) / RMSMax );
  
}

TF1* NaiveDirection(TH2F* PlotSel,Float_t XIP,Float_t YIP){
  
  TF1* direction = new TF1("direction","[0]*(x-[1])+[2]",2304);
  direction->FixParameter(1,XIP);
  direction->FixParameter(2,YIP);
  PlotSel->Fit(direction,"Q");

  return direction;
  
}

TF1* NaiveDirection(TGraphErrors* PlotSel,Float_t XIP,Float_t YIP){
  
  TF1* direction = new TF1("direction","[0]*(x-[1])+[2]",0,2304);
  direction->FixParameter(1,XIP);
  direction->FixParameter(2,YIP);
  PlotSel->Fit(direction,"Q");

  return direction;
  
}



void TGraphFitRebin(TH2F* histo, TGraphErrors* graph, Int_t reb,Float_t w){

  TH2F* historeb = new TH2F("historeb","historeb",2304,0,2304,2304,0,2304);

  for(int i=0; i<histo->GetXaxis()->GetNbins();i++){
    for(int k=0; k<histo->GetYaxis()->GetNbins();k++){
      if(histo->GetBinContent(i,k)!=0){
	historeb->SetBinContent(i,k,histo->GetBinContent(i,k));
      }
    }
  }

  historeb->Rebin2D(reb,reb);
  
  int count=0;
  Float_t z;
  
  for(int i=0; i<historeb->GetXaxis()->GetNbins();i++){
    for(int k=0; k<historeb->GetYaxis()->GetNbins();k++){
      z=historeb->GetBinContent(i,k);
      if(z > 0){
	graph->SetPoint(count,historeb->GetXaxis()->GetBinCenter(i),historeb->GetYaxis()->GetBinCenter(k));
	graph->SetPointError(count,w/z,w/z);
	count++;
      }
      
    }//chiudo for k   
  }//chiudo for i
  
  delete historeb;
}

Float_t SupprFactor(Double_t x,Double_t y,Double_t XBar,Double_t YBar,Double_t w){
  
  Float_t dist;
  dist=sqrt( (x-XBar)*(x-XBar) + (y-YBar)*(y-YBar) );

  return exp(-(dist/w) );
}


Float_t CorrAngle(Float_t Angle,Float_t XIP,Float_t YIP,Float_t XBC,Float_t YBC){
  if(XIP<XBC){
    return Angle;

  }else {

    if( (Angle+TMath::Pi()) > TMath::Pi() ){

      return (Angle-TMath::Pi());

    } else {

      return (Angle+TMath::Pi());

    }
  }
}


Float_t ImprCorrAngle(Float_t Angle,Float_t XIP,Float_t YIP,Float_t XPIP,Float_t YPIP){
  Float_t AnglePerp = Angle+TMath::Pi()/2;

  Float_t qIP= YIP-TMath::Tan(AnglePerp)*XIP;
  Float_t qPIP= YPIP-TMath::Tan(AnglePerp)*XPIP;

  if(Angle>0){
    if(qIP<qPIP) {
      return Angle;
    } else {
      return Angle-TMath::Pi();
    }//chiudo else qIP
  } else {
    if(qIP<qPIP){
      return Angle+TMath::Pi();
    } else {
      return Angle;
    }//chiudo else qIP
  }//chiudo esle angolo

}//chiudo Funzione



void ApplyThr(TH2F* Track,Float_t Thr){

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();

  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();

  Float_t z;

  for(int i=XBinMin; i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      z=Track->GetBinContent(i,j);
      if(z>0 && z<Thr){
	Track->SetBinContent(i,j,0);
      }
    }
  }
}

void RemoveNoise(TH2F* Track){

  Int_t XBinMin=Track->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+Track->GetXaxis()->GetNbins();

  Int_t YBinMin=Track->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+Track->GetYaxis()->GetNbins();

  std::vector<Int_t> ToRemoveXBin;
  std::vector<Int_t> ToRemoveYBin;

  Float_t ETemp;

  for(int i=XBinMin; i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      if(Track->GetBinContent(i,j)>0){
	ETemp=0;
	if(Track->GetBinContent(i+1,j)>=0) ETemp+=Track->GetBinContent(i+1,j);
	if(Track->GetBinContent(i-1,j)>=0) ETemp+=Track->GetBinContent(i-1,j);
	if(Track->GetBinContent(i,j+1)>=0) ETemp+=Track->GetBinContent(i,j+1);
	if(Track->GetBinContent(i,j-1)>=0) ETemp+=Track->GetBinContent(i,j-1);
	//std::cout <<"EAround  "<< ETemp << std::endl;
	if(ETemp<=0){
	  ToRemoveXBin.push_back(i);
	  ToRemoveYBin.push_back(j);
	}//chiudo for E
      }//Chiudo if Z>0
    }//chiudo for j
  }//chiudo for j

  // std::cout << ToRemoveXBin.size() << std::endl;
  //std::cout << ToRemoveXBin.size()<< std::endl;
  for(int i=0; i<ToRemoveXBin.size(); i++){
    Track->SetBinContent(ToRemoveXBin[i],ToRemoveYBin[i],0);
  }
}


double GetThetaTrue(double py, double pz){
  
  double RowAngle = TMath::ATan( (py)/(pz) );

  if(pz>0) {
    return RowAngle;
  } else{
    if (RowAngle>0){
      return RowAngle-TMath::Pi();
    } else {
      return RowAngle+TMath::Pi();
    }
  }
  
}

double GetPhiTrue(double px,double py,double pz){
  
  return TMath::ATan(px/(sqrt(py*py+pz*pz)));
  
}

Int_t FlagIP(TH2F* Track){

  Int_t flag=0;

  for(int i=-1;i<2;i++){
    for(int j=-1;j<2;j++){
      if(Track->GetBinContent(Track->GetXaxis()->FindBin(1152+i),Track->GetYaxis()->FindBin(1152+j))>0){
	flag=1;
      }//chiudo if
    }//chiudo for j
  }//chiudo for i

  return flag;
}

int main(int argc, char** argv){
  //Read Variables from File
  std::ifstream inputFile(argv[2]);
  std::string line;

  Float_t rmin;
  Float_t rmin_init;
  Float_t rmax;
  Float_t RSel;
  Float_t NPointSelMin;
  Float_t NPointSelMin_init;
  Float_t ThresholdPoints;
  Float_t Wfac;

  while (getline(inputFile, line)){
    std::istringstream ss(line);

    std::string name;
    ss >> rmin >> rmax >> RSel>> NPointSelMin >> ThresholdPoints >> Wfac;

  }

  std::cout <<"Param: " <<"rmin" << "\t" << "rmax" << "\t" << "RSel" << "\t" << "NPIP" << "\t"<< "Thr"<<"\t"  << "Wfac"<< std::endl;
  std::cout <<"Param: " <<rmin << "\t" << rmax << "\t" << RSel << "\t" << NPointSelMin<<"\t" <<ThresholdPoints << "\t" << Wfac<< std::endl;

  rmin_init=rmin;
  NPointSelMin_init=NPointSelMin;

  //Get Input ROOT File
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");

  bool fullpixel=false;
  //Set Branches and Define Variables
  Int_t nmax=200000;

  Int_t npixel=2304;
  Int_t XMin;
  Int_t XMax;
  Int_t YMin;
  Int_t YMax;


  UInt_t nSc;

  //Pixelsxquer
  Int_t ScNpixels[26];
  Float_t ScIntegral[26];
  Float_t XPix[nmax];
  Float_t YPix[nmax];
  Float_t ZPix[nmax];

  Float_t MCPx[26];
  Float_t MCPy[26];
  Float_t MCPz[26];

  Float_t MCx_beg[26];
  Float_t MCy_beg[26];
  Float_t MCz_beg[26];

  Float_t MCx_end[26];
  Float_t MCy_end[26];
  Float_t MCz_end[26];

  Float_t MC_zhits[700];
  Float_t MC_yhits[700];


  Float_t MC_2D_length[26];

  Float_t ThetaTrue[26];
  Float_t PhiTrue[26];
  Float_t X,Y,Z;


  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_integral",ScIntegral);

  tree->SetBranchAddress("sc_nintpixels",ScNpixels);
  tree->SetBranchAddress("sc_xpixelcoord",YPix);
  tree->SetBranchAddress("sc_ypixelcoord",XPix);
  tree->SetBranchAddress("sc_zpixel",ZPix);

  tree->SetBranchAddress("MC_x_vertex",MCx_beg);

  tree->SetBranchAddress("MC_x_vertex_end",MCx_end);
  tree->SetBranchAddress("MC_y_vertex_end",MCy_end);
  tree->SetBranchAddress("MC_z_vertex_end",MCz_end);

  tree->SetBranchAddress("MC_zhits",MC_zhits);
  tree->SetBranchAddress("MC_yhits",MC_yhits);


  tree->SetBranchAddress("MC_px",MCPx);
  tree->SetBranchAddress("MC_py",MCPy);
  tree->SetBranchAddress("MC_pz",MCPz);

  tree->SetBranchAddress("MC_2D_pathlength",MC_2D_length);

  //Analysis Variables

  std::vector<Int_t> BeginScPix;
  std::vector<Int_t> EndScPix;

  std::vector<Int_t> BeginScSl;
  std::vector<Int_t> EndScSl;

  std::vector<Int_t> OnMainTrack;

  Int_t EvNumber;

  Float_t XBar[26];
  Float_t YBar[26];
  Float_t IPFlag[26];
  Float_t PhiAxisTrack[26];
  Float_t SkewAxisTrack[26];
  Float_t RMSOnLine[26];
  Float_t RMSOnLinePerp[26];
  Float_t RMSRatio[26];

  std::vector<Float_t> XSel;
  std::vector<Float_t> YSel;
  std::vector<Float_t> ZSel;
  std::vector<Float_t> ZSelScaled;

  Float_t PointSkew;
  Float_t PointDistCm;

  Int_t NMinPix=300;

  Int_t NSelPoints[26];
  Float_t XIntPoint[26];
  Float_t YIntPoint[26];

  Float_t XIPPrev;
  Float_t YIPPrev;
  std::vector<Float_t> XIntPointPrev;
  std::vector<Float_t> YIntPointPrev;

  Float_t AngleValFit[26];
  Float_t AngleValGraph[26];
  Float_t AngleValGraphRebin[26];
  Float_t AngleValRMS[26];

  Float_t AngleValFitDist[26];
  Float_t AngleValGraphDist[26];
  Float_t AngleValRMSDist[26];

  Float_t AngleDiff[26];

  Float_t RMSSelOnLine[26];
  Float_t RMSSelOnLinePerp[26];
  Float_t RMSSelRatio[26];

  TH2F* IPDistrib = new TH2F("IPDistrib","IPDistrib",500,0,2304,500,0,2304);
  TH1F* RMSRatioDist = new TH1F("RMSRatioDist","RMSRatioDist",100,0,1);


  TH1F* AngleFitDistrib= new TH1F("AngleFitDistrib","AngleFitDistrib",100,-TMath::Pi(),TMath::Pi());
  TH1F* AngleGraphDistrib= new TH1F("AngleGraphDistrib","AngleGraphDistrib",100,-TMath::Pi()-1,TMath::Pi()+1);
  TH1F* AngleGraphRebinDistrib= new TH1F("AngleGraphRebinDistrib","AngleGraphRebinDistrib",100,-TMath::Pi(),TMath::Pi());
  TH1F* AngleRMSDistrib= new TH1F("AngleRMSDistrib","AngleRMSDistrib",100,-3.16,3.16);

  TH1F* AngleFitDist= new TH1F("AngleFitDist","AngleFitDist",100,-TMath::Pi(),TMath::Pi());
  TH1F* AngleFitGraphDist= new TH1F("AngleFitGraphDist","AngleFitGraphDist",100,-TMath::Pi(),TMath::Pi());
  TH1F* AngleRMSDist= new TH1F("AngleRMSDist","AngleRMSDist",100,-3.5,3.5);
  TH1F* XIPDistrib= new TH1F("XIPDistrib","XIPDistrib",50,-10,10);
  TH1F* YIPDistrib= new TH1F("YIPDistrib","YIPDistrib",50,-10,10);

  //TreeToWriteOn
  TFile* TreeFile = new TFile("Out/Ntuples.root","recreate");
  TTree* treeWrite = new TTree("myTree","");

  treeWrite->Branch("EvNumber", &EvNumber, "EvNumber/I");
  treeWrite->Branch("nSc", &nSc, "nSc/I");
  treeWrite->Branch("ScNpixels",ScNpixels,"ScNpixels[nSc]/I");
  treeWrite->Branch("ScIntegral",ScIntegral,"ScIntegral[nSc]/F");

  treeWrite->Branch("XBar",XBar,"XBar[nSc]/F");
  treeWrite->Branch("YBar",YBar,"YBar[nSc]/F");

  treeWrite->Branch("PhiAxisTrack",PhiAxisTrack,"PhiAxisTrack[nSc]/F");
  treeWrite->Branch("SkewAxisTrack",SkewAxisTrack,"SkewAxisTrack[nSc]/F");
  treeWrite->Branch("RMSOnLine",RMSOnLine,"RMSOnLine[nSc]/F");
  treeWrite->Branch("RMSOnLinePerp",RMSOnLinePerp,"RMSOnLinePerp[nSc]/F");

  treeWrite->Branch("NSelPoints",NSelPoints,"NSelPoints[nSc]/I");
  treeWrite->Branch("XIntPoint",XIntPoint,"XIntPoint[nSc]/F");
  treeWrite->Branch("YIntPoint",YIntPoint,"YIntPoint[nSc]/F");
  treeWrite->Branch("IPFlag",IPFlag,"IPFlag[nSc]/F");

  treeWrite->Branch("MCx_beginnig",MCx_beg,"MCx_beginning[nSc]/F");
  treeWrite->Branch("MCy_beginnig",MCy_beg,"MCy_beginning[nSc]/F");
  treeWrite->Branch("MCz_beginnig",MCz_beg,"MCz_beginning[nSc]/F");

  treeWrite->Branch("MCx_end",MCx_end,"MCx_end[nSc]/F");
  treeWrite->Branch("MCy_end",MCy_end,"MCy_end[nSc]/F");
  treeWrite->Branch("MCz_end",MCz_end,"MCz_end[nSc]/F");

  treeWrite->Branch("MC_2D_length",MC_2D_length,"MC_2D_length[nSc]/F");

  treeWrite->Branch("AngleValFit",AngleValFit,"AngleValFit[nSc]/F");
  treeWrite->Branch("AngleValGraph",AngleValGraph,"AngleValGraph[nSc]/F");
  treeWrite->Branch("AngleValGraphReb",AngleValGraph,"AngleValGraphReb[nSc]/F");
  treeWrite->Branch("AngleValRMS",AngleValRMS,"AngleValRMS[nSc]/F");

  treeWrite->Branch("AngleValFitDist",AngleValFitDist,"AngleValFitDist[nSc]/F");
  treeWrite->Branch("AngleValGraphDist",AngleValGraph,"AngleValGraphDist[nSc]/F");
  treeWrite->Branch("AngleValRMSDist",AngleValRMSDist,"AngleValRMSDist[nSc]/F");
  treeWrite->Branch("ThetaTrue",ThetaTrue,"ThetaTrue[nSc]/F");
  treeWrite->Branch("PhiTrue",PhiTrue,"PhiTrue[nSc]/F");

  treeWrite->Branch("AngleDiffCorr",AngleDiff,"AngleDiffCorr[nSc]/F");

//Analysis
  Int_t Nevents=tree->GetEntries();
  Float_t ZMin;

  Float_t TrackRotationAngle=0;
  Float_t MidImagePix=1152;

  //Nevents=3000;

  for(int k=0;k<Nevents;k++){

    tree->GetEntry(k);


    //EvNumber = k;
    //OnMainTrack=SortSlices(nSc,Sc_nSl,Sc_Profx,Sc_Profy,Sc_Energyprof);
    //std::cout << "Nev: "<< k << std::endl;


    if(fullpixel){
      ZMin=-TMath::MinElement(ScNpixels[0],ZPix);
      std::cout << ZMin << std::endl;
    } else{
      ZMin=0;
    }

    TH2F* Track[nSc];
    //TH2F* Slices[nSc];
    TH2F* SelectedIP[nSc];
    TH2F* SelectionScaled[nSc];
    //TF1* DirectionFit[nSc];
    TF1* DirectionFitRMS[nSc];
    TF1* DirectionFitGraph[nSc];
    //TF1* DirectionFitGraphRebin[nSc];
    TGraphErrors* GraphSelIP[nSc];
    TGraphErrors* GraphSelIPRebin[nSc];

//Distance scaled TH
    TH2F* SelIPDistance[nSc];
    TF1* DirectionFitDist[nSc];
    TGraphErrors* GraphSelIPDistance[nSc];
    TF1* DirectionFitGraphDist[nSc];
    TF1* DirectionFitRMSDist[nSc];

    Float_t rminN[nSc];

    ScIndicesElem(nSc,ScNpixels,&BeginScPix,&EndScPix);
    /*
    for (int i=0;i<nSc;i++){
      std::cout<<"Start"<<std::endl;
      std::cout<<nSc<<"\t CULO \t"<< BeginScPix[i]<<"\t CULO \t"<< EndScPix[i]<<std::endl;
      std::cout<<nSc<<"\t CULO \t"<< XPix[BeginScPix[i]]<<"\t CULO \t"<< XPix[EndScPix[i]]<<std::endl;
      std::cout<<EndScPix[i]-BeginScPix[i]<<std::endl;
      }
    break;
    */
    for(int i= 0;i<nSc;i++){
      XMin = TMath::MinElement(EndScPix[i]-BeginScPix[i],&XPix[BeginScPix[i]]);
      XMax = TMath::MaxElement(EndScPix[i]-BeginScPix[i],&XPix[BeginScPix[i]]);
      YMin = TMath::MinElement(EndScPix[i]-BeginScPix[i],&YPix[BeginScPix[i]]);
      YMax = TMath::MaxElement(EndScPix[i]-BeginScPix[i],&YPix[BeginScPix[i]]);

      if(ScIntegral[i]<10000) continue;

      //std::cout << MC_zhits[0] << "  " << MC_yhits[0] << std::endl;
      //std::cout << MCz_beg[0] << "  " << MCy_beg[0] << std::endl;
      //std::cout << XMin << "\t" << XMax<< "\t" <<YMin<< "\t" <<YMax << std::endl;

      //XMax=(XMax-MidImagePix)*cos(TrackRotationAngle) - (YMax-MidImagePix)*sin(TrackRotationAngle) + MidImagePix;
      //YMax=(YMax-MidImagePix)*cos(TrackRotationAngle) - (XMax-MidImagePix)*sin(TrackRotationAngle) + MidImagePix;

      //XMin=(XMin-MidImagePix)*cos(TrackRotationAngle) - (YMin-MidImagePix)*sin(TrackRotationAngle) + MidImagePix;
      //YMin=(YMin-MidImagePix)*cos(TrackRotationAngle) - (XMin-MidImagePix)*sin(TrackRotationAngle) + MidImagePix;

      Int_t IncrImgSize=30;

      Track[i]= new TH2F(Form("Track%i",k),Form("Track%i",i),XMax-XMin+2*IncrImgSize,XMin-IncrImgSize,XMax+IncrImgSize,YMax-YMin+2*IncrImgSize,YMin-IncrImgSize,YMax+IncrImgSize);
      SelectedIP[i]= new TH2F(Form("SelectedIP%i",k),Form("SelectedIP%i",i),XMax-XMin+2*IncrImgSize,XMin-IncrImgSize,XMax+IncrImgSize,YMax-YMin+2*IncrImgSize,YMin-IncrImgSize,YMax+IncrImgSize);
      SelectionScaled[i]=new TH2F(Form("SelectionScaled%i",k),Form("SelectionScaled%i",i),XMax-XMin+2*IncrImgSize,XMin-IncrImgSize,XMax+IncrImgSize,YMax-YMin+2*IncrImgSize,YMin-IncrImgSize,YMax+IncrImgSize);
      SelIPDistance[i]= new TH2F(Form("SelIPDistance%i",k),Form("SelIPDistance%i",i),XMax-XMin+2*IncrImgSize,XMin-IncrImgSize,XMax+IncrImgSize,YMax-YMin+2*IncrImgSize,YMin-IncrImgSize,YMax+IncrImgSize);


      XSel.clear();
      YSel.clear();
      ZSel.clear();
      ZSelScaled.clear();
      XIntPointPrev.clear();
      YIntPointPrev.clear();

      //Construct Histogram with single Supercluster
      for(int j=BeginScPix[i];j<EndScPix[i];j++){
	Track[i]->SetBinContent(Track[i]->GetXaxis()->FindBin( (XPix[j]-MidImagePix)*cos(TrackRotationAngle) - (YPix[j]-MidImagePix)*sin(TrackRotationAngle) + MidImagePix),Track[i]->GetYaxis()->FindBin( (XPix[j]-MidImagePix)*sin(TrackRotationAngle) + (YPix[j]-MidImagePix)*cos(TrackRotationAngle) + MidImagePix ),ZPix[j]+ZMin);
      //std::cout<<ZPix[j]+ZMin<<std::endl;
      }//chiudo for j (fill histos)
      //Track[i]->Rebin2D(2,2);
      //Track[i]->Smooth();

      //Put to zro evry bin that contain less conent than ThresholdPoints
      ApplyThr(Track[i],ThresholdPoints);
      //check every bin that has content >0. If ALL the nearest (4 points) x+/-1 and y+/-1 has no bin content the initial point is set to zero
      RemoveNoise(Track[i]);

      IPFlag[i]=FlagIP(Track[i]);

      //Calculation of Baricenter, Phi of the line which maximize RMS, Skewness of the track on the Max RMS line
      Barycenter(Track[i],&XBar[i],&YBar[i]);
      PhiAxisTrack[i]= AngleLineMaxRMS(Track[i],XBar[i], YBar[i], &RMSOnLine[i],&RMSOnLinePerp[i]);
      RMSRatio[i]=RMSOnLinePerp[i]/RMSOnLine[i];
      SkewAxisTrack[i]=SkewOnLine(Track[i],XBar[i], YBar[i],PhiAxisTrack[i]);

      //Construction of the histogram with the beginning of the track and IP Calculation
      rmin_init=0.2;
      rmin=rmin_init;
      NPointSelMin=NPointSelMin_init;
      rminN[i]=rmin;
      NSelPoints[i]=800;

      //SelectedIP[i]->Rebin2D(2,2);

      //      std::cout << "here4"	<< std::endl;

      do{
	//std::cout << rminN[i] << "   " << NSelPoints[i] << " " << NPointSelMin <<"  " << rmin <<  std::endl;
	rminN[i]+=0.03;
	XSel.clear();
	YSel.clear();
	ZSel.clear();

	//XIntPointPrev

	Barycenter(SelectedIP[i],&XIPPrev,&YIPPrev);
	XIntPointPrev.push_back(XIPPrev);
	YIntPointPrev.push_back(YIPPrev);

	//	std::cout <<"IPTemp \t" <<XIPPrev << "\t" << YIPPrev << std::endl;

	SelectedIP[i]->Reset();

	for(int j=0;j<Track[i]->GetXaxis()->GetNbins();j++){
	  for(int l=0;l<Track[i]->GetYaxis()->GetNbins();l++){

	    X=Track[i]->GetXaxis()->GetBinCenter(j);
	    Y=Track[i]->GetYaxis()->GetBinCenter(l);
	    Z=Track[i]->GetBinContent(j,l);

	    PointSkew=GetPointSkew(X,Y,XBar[i],YBar[i],PhiAxisTrack[i],SkewAxisTrack[i]);
	    PointDistCm=DistCm(X,Y,XBar[i],YBar[i],RMSOnLine[i]);

	    if(PointSkew>0 && PointDistCm>rminN[i] && Z>0){
	      SelectedIP[i]->SetBinContent(SelectedIP[i]->GetXaxis()->FindBin(X),SelectedIP[i]->GetYaxis()->FindBin(Y),Z);
	      XSel.push_back(X);
	      YSel.push_back(Y);
	      ZSel.push_back(Z);
	    }//chiudo if selection

	  }//chiuso for l
	}//chiudo for j (fill histos)
	NSelPoints[i]=XSel.size();
	//	std::cout << rminN[i] << "\t"<< NSelPoints[i] << std::endl;

      }while(NSelPoints[i]>NPointSelMin);

      //std::cout << "outofwhile" << std::endl;

      //Barycenter of the selected begin of the track
      Barycenter(SelectedIP[i],&XIntPoint[i],&YIntPoint[i]);

      // std::cout << "IPDef \t" << XIntPoint[i] << "\t" << YIntPoint[i] << std::endl; 

      XSel.clear();
      YSel.clear();
      ZSel.clear();

      //SelectionScaled[i]->Rebin2D(2,2);

      //circular selection
      for(int j=0;j<Track[i]->GetXaxis()->GetNbins();j++){
	for(int l=0;l<Track[i]->GetYaxis()->GetNbins();l++){
    X=Track[i]->GetXaxis()->GetBinCenter(j);
    Y=Track[i]->GetYaxis()->GetBinCenter(l);
    Z=Track[i]->GetBinContent(j,l);

    SelectionScaled[i]->SetBinContent(SelectionScaled[i]->GetXaxis()->FindBin(X),SelectionScaled[i]->GetYaxis()->FindBin(Y),Z*SupprFactor(X,Y,XIntPoint[i],YIntPoint[i],Wfac));

	}//chiudo l
      }//chiudo for j

      //recalculation of the barycenter
      //Barycenter(0,XSel.size(),&XSel[0],&YSel[0],&ZSel[0],&XIntPoint[i],&YIntPoint[i]);
      //Barycenter(SelectionScaled[i],&XIntPoint[i],&YIntPoint[i]);
      //Countruction of the rebinned TGraphErrors
      //TGraphFitRebin(SelectedIP[i],GraphSelIP[i],1,8);
      //TGraphFitRebin(SelectedIP[i],GraphSelIPRebin[i],4,8);

      //Fit of the selectedip region
      // DirectionFit[i]=NaiveDirection(SelectedIP[i],XIntPoint[i],YIntPoint[i]);
      // DirectionFit[i]->SetName("DirectionFit");

      //DirectionFitGraph[i]=NaiveDirection(GraphSelIP[i],XIntPoint[i],YIntPoint[i]);
      //DirectionFitGraph[i]->SetName("DirectionFitGraph");

      // DirectionFitGraphRebin[i]=NaiveDirection(GraphSelIPRebin[i],XIntPoint[i],YIntPoint[i]); 
      // DirectionFitGraphRebin[i]->SetName("DirectionFitGraphRebin");

      //Direction of the track with the RMS method
      AngleValRMS[i]=AngleLineMaxRMS(SelectionScaled[i], XIntPoint[i], YIntPoint[i], &RMSSelOnLine[i], &RMSSelOnLinePerp[i]); 
      DirectionFitRMS[i] = new TF1("FitRMS","[0]*(x-[1])+[2]",0,2304);

      /////////////////////////////Analysis with exponential decrease of selected IP////////////////////
      /*for(int j=BeginScPix[i];j<EndScPix[i];j++){
	if(ZPix[j]+ZMin>0){
	  PointSkew=GetPointSkew(XPix[j],YPix[j],XBar[i],YBar[i],PhiAxisTrack[i],SkewAxisTrack[i]);
	  PointDistCm=DistCm(XPix[j],YPix[j],XBar[i],YBar[i],RMSOnLine[i]);
	  if(PointSkew>0 && PointDistCm>rmin && PointDistCm<rmax){
	    ZSelScaled.push_back((ZPix[j]+ZMin)*SupprFactor(XPix[j],YPix[j],XIntPoint[i],YIntPoint[i],Wfac));
	  }//chiudo if selection
	}//chiudo if Z
      }//chiudo for pixels
      */

      //GraphSelIPDistance[i] = new TGraphErrors();
      //GraphSelIPDistance[i]->SetName(Form("GraphSelIPDistance%i",i));
      //GraphSelIPDistance[i]->SetTitle(GraphSelIPDistance[i]->GetName());

      // DirectionFitDist[i]=NaiveDirection(SelIPDistance[i],XIntPoint[i],YIntPoint[i]);
      //TGraphFitRebin(SelIPDistance[i],GraphSelIPDistance[i],1,1);
      //DirectionFitGraphDist[i]=NaiveDirection(GraphSelIPDistance[i],XIntPoint[i],YIntPoint[i]);

      AngleValRMSDist[i]=AngleLineMaxRMS(SelectionScaled[i], XIntPoint[i], YIntPoint[i], &RMSSelOnLine[i],&RMSSelOnLinePerp[i]);
      //std::cout << AngleValRMSDist[i] << std::endl;
      DirectionFitRMSDist[i] = new TF1("FitRMSDist","[0]*(x-[1])+[2]",0,2304);

      //Plot of the first Supercluster (main one) with Analysis components
      //if(i==0 && k%5==0)MakeCanvas(Track[i],XBar[i],YBar[i],PhiAxisTrack[i],k,i,RMSOnLine[i],rminN[i],rmax,XIntPoint[i],YIntPoint[i],AngleValRMS[i],DirectionFitGraph[i]);


      if(ScIntegral[i]>7000){

	//AngleValGraph[i] =TMath::ATan(DirectionFitGraph[i]->GetParameter(0));
	//AngleValGraph[i] = CorrAngle(AngleValGraph[i],XIntPoint[i],XBar[i]);
	//AngleGraphDistrib->Fill(AngleValGraph[i]);

	MCz_beg[i]=1152+MC_yhits[0]/(346/2)*1152;
	MCy_beg[i]=1152+MC_zhits[0]/(346/2)*1152;

	XIPDistrib->Fill((XIntPoint[i]-MCy_beg[i])*0.151);
	YIPDistrib->Fill((YIntPoint[i]-MCz_beg[i])*0.151);


	AngleValRMS[i] = ImprCorrAngle(AngleValRMS[i],XIntPoint[i],YIntPoint[i],XIntPointPrev[XIntPointPrev.size()/2],YIntPointPrev[YIntPointPrev.size()/2]); 
	DirectionFitRMS[i]->FixParameter(1,XIntPoint[i]);
	DirectionFitRMS[i]->FixParameter(2,YIntPoint[i]);
	DirectionFitRMS[i]->FixParameter(0,TMath::Tan(AngleValRMS[i]));
	AngleRMSDistrib->Fill(AngleValRMS[i]);

	//AngleValGraphDist[i]=TMath::ATan(DirectionFitGraphDist[i]->GetParameter(0));
	//AngleValGraphDist[i] = CorrAngle(AngleValGraphDist[i],XIntPoint[i],XBar[i]);
	//AngleFitGraphDist->Fill(AngleValGraphDist[i]);

	AngleValRMSDist[i] = ImprCorrAngle(AngleValRMSDist[i],XIntPoint[i],YIntPoint[i],XIntPointPrev[XIntPointPrev.size()/2-1],YIntPointPrev[YIntPointPrev.size()/2-1]);
	//AngleValRMSDist[i] = CorrAngle(AngleValRMSDist[i],XIntPoint[i],XBar[i]);
	DirectionFitRMSDist[i]->FixParameter(1,XIntPoint[i]);
	DirectionFitRMSDist[i]->FixParameter(2,YIntPoint[i]);
	DirectionFitRMSDist[i]->FixParameter(0,TMath::Tan(AngleValRMSDist[i]));

	ThetaTrue[i]=GetThetaTrue(MCPy[i],MCPz[i]);
	PhiTrue[i]=GetPhiTrue(MCPx[i],MCPy[i],MCPz[i]);
	//std:: cout << "Angle=" << AngleValRMSDist[i]/TMath::Pi()*180 << "\t Real=" << ThetaTrue[i]/TMath::Pi()*180 << "\t Diff=" << AngleValRMSDist[i]/TMath::Pi()*180 - ThetaTrue[i]/TMath::Pi()*180 << std::endl;

	if( abs(AngleValRMSDist[i] - ThetaTrue[i] ) <TMath::Pi() ) {
	  AngleRMSDist->Fill(  (AngleValRMSDist[i] - ThetaTrue[i] ) );
	  AngleDiff[i]=(AngleValRMSDist[i] - ThetaTrue[i] );
	  //std::cout << "Normal Fill" << std::endl;
	}else if((AngleValRMSDist[i] - ThetaTrue[i])>TMath::Pi()) {
	  AngleRMSDist->Fill(  2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i])  );
	  AngleDiff[i]=2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i]);
	  //std::cout << "Greater than Pi diff is 2pi-phi = " << (2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i]))/TMath::Pi()*180 << std::endl;
	}else if((AngleValRMSDist[i] - ThetaTrue[i])<-TMath::Pi()){
	  AngleRMSDist->Fill(  -2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i])  );
	  AngleDiff[i]=-2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i]);
	  //std::cout << "Less than -Pi diff is -2pi-phi = " << (-2*TMath::Pi()- (AngleValRMSDist[i] - ThetaTrue[i]))/TMath::Pi()*180 << std::endl;
	}
	//std::cout << AngleValRMSDist[i] << std::endl;

	IPDistrib->Fill(XIntPoint[i],YIntPoint[i]);
	RMSRatioDist->Fill(RMSRatio[i]);

	////////////////////////////////////Distance scaled plots//////////////////////

      } else {

	//AngleValFit[i]=-9999;
	AngleValGraph[i]=-9999;
	//AngleValGraphRebin[i]=-9999;
	AngleValRMS[i]=-9999;

      }

      //std::cout << XIntPointPrev[XIntPointPrev.size()/2] << std::endl;

      if(i==0 && k%10==0) MakeCanvasDist(Track[i],XBar[i],YBar[i],PhiAxisTrack[i],k,i,RMSOnLine[i],rminN[i],rmax,XIntPoint[i],YIntPoint[i],AngleValRMSDist[i],RSel,XIntPointPrev[(int)(XIntPointPrev.size()/2)],YIntPointPrev[(int)(YIntPointPrev.size()/2)],ThetaTrue[i]/TMath::Pi()*180);
	  }//chiudo for i (loop on the superclusters)

    treeWrite->Fill();

    ////////////////////////Save Plot On File////////////////////////////////////////
    /*if(k%10==0){
      TFile* fout = new TFile(Form("Out/Tracks_Ev%i.root",k),"recreate");
      fout->cd();
      Track[0]->Write();
      SelectedIP[0]->Write();
      SelectionScaled[0]->Write();
      //GraphSelIP[i]->Write();

      //GraphSelIPRebin[i]->Write();
      SelIPDistance[0]->Write();
      //GraphSelIPDistance[i]->Write();
      //DirectionFitRMS[i]->Write();
      DirectionFitRMSDist[0]->Write();

      fout->Save();
      fout->Close();
      }*/

  }//chiudo for k (loop on the events)

//Save Tree
  gStyle->SetOptStat(1111);
  /*
  TreeFile->cd();
  treeWrite->Write();
  TreeFile->Save();
  TreeFile->Close();
  */
//Canvas With Resolution
  TCanvas* pippo = new TCanvas("AngRes","AngRes",1200,1200);
  pippo->Divide(2,3);
  pippo->cd(1);
  AngleGraphDistrib->Draw();
  pippo->cd(2);
  AngleRMSDistrib->Draw();
  pippo->cd(3);
  IPDistrib->Draw("COLZ");
  pippo->cd(4);
  AngleFitGraphDist->Draw();
  pippo->cd(5);
  AngleRMSDist->Draw();
  pippo->SaveAs("Out/AngDist.png");


  TFile* f3 = new TFile("Out/ResoPlot.root","recreate");
  f3->cd();
  AngleFitDistrib->Write();
  AngleGraphDistrib->Write();
  AngleRMSDistrib->Write();
  AngleGraphRebinDistrib->Write();
  AngleFitDist->Write();
  AngleFitGraphDist->Write();
  AngleRMSDist->Write();
  XIPDistrib->Write();
  YIPDistrib->Write();
  IPDistrib->Write();
  RMSRatioDist->Write();
  f3->Save();
  f3->Close();


  TF1* gausfit = new TF1("gausfit","[0]*exp(-0.5*( (x-[1])/[2] )* ( (x-[1])/[2]) )+[3]",-TMath::Pi(),TMath::Pi());
  gausfit->SetParLimits(0,0,700);
  gausfit->SetParLimits(1,-1.5,1.5);
  gausfit->SetParLimits(2,0,1);
  gausfit->SetParLimits(3,0,100);
  gausfit->SetParameters(AngleRMSDist->GetMaximum()-AngleRMSDist->GetBinContent(15),0,0.4,AngleRMSDist->GetBinContent(15));
  AngleRMSDist->Fit(gausfit,"QR");

  TF1* gausfitXIP = new TF1("gausfitXIP","gaus");
  gausfitXIP->SetParLimits(0,0,3000);
  gausfitXIP->SetParLimits(1,-5,5);
  gausfitXIP->SetParLimits(2,0,15);
  XIPDistrib->Fit(gausfitXIP,"Q");

  TF1* gausfitYIP = new TF1("gausfitYIP","gaus");
  gausfitYIP->SetParLimits(0,0,3000);
  gausfitYIP->SetParLimits(1,-5,5);
  gausfitYIP->SetParLimits(2,0,10);
  YIPDistrib->Fit(gausfitYIP,"Q");

  std::cout << "here" << std::endl;

  std::cout<< gausfit->GetParameter(2)<<"\t" << gausfit->GetParError(2) << "\t" << gausfitXIP->GetParameter(2) << "\t"<< gausfitXIP->GetParError(2) << "\t"<< gausfitYIP->GetParameter(2) << "\t"<< gausfitYIP->GetParError(2) << "\t" << AngleRMSDist->Integral(AngleRMSDist->GetXaxis()->FindBin(-1.5707),AngleRMSDist->GetXaxis()->FindBin(1.5707))/AngleRMSDist->Integral()*100<< "\t" << gausfit->GetParameter(3) << "\t" << gausfit->GetParError(3)<< std::endl;

  std::cout << "herelast" << std::endl;

  return 0;

}




void MakeCanvas(TH2F* Plot,Float_t XBar,Float_t YBar,Float_t Phi,Int_t i,Int_t k,Float_t RMS,Float_t rmin, Float_t rmax,Float_t XIP,Float_t YIP,Float_t AngleDir,TF1* DirGraph){

  Double_t RangeX = Plot->ProjectionX()->GetRMS();
  Double_t RangeY = Plot->ProjectionY()->GetRMS();

  TLegend* L = new TLegend();

  TGraph* BarPlot = new TGraph();
  BarPlot->SetName("BarPlot");
  BarPlot->SetTitle("BarPlot");
  
  BarPlot->SetPoint(0,XBar,YBar);
  BarPlot->SetMarkerStyle(8);
  BarPlot->SetMarkerSize(1.);

  TGraph* IPPlot = new TGraph();
  IPPlot->SetName("IPPlot");
  IPPlot->SetTitle("IPPlot");
  
  IPPlot->SetPoint(0,XIP,YIP);
  IPPlot->SetMarkerStyle(8);
  IPPlot->SetMarkerColor(kRed+2);
  IPPlot->SetMarkerSize(1.);

  TF1* MainAxis = new TF1("MainAxis","[0]*(x-[1])+[2]",0,2304);
  MainAxis->SetParameter(0,TMath::Tan(Phi));
  MainAxis->SetParameter(1,XBar);
  MainAxis->SetParameter(2,YBar);

  TArc* arcMin = new TArc(XBar,YBar,rmin*sqrt(RMS));
  arcMin->SetFillColorAlpha(kBlue, 0);

  TArc* arcMax = new TArc(XBar,YBar,rmax*sqrt(RMS));
  arcMax->SetFillColorAlpha(kBlue, 0);

  
  //DirGraph->SetLineColor(kGray);
      
  TF1* DirRMS = new TF1("DirRMS","[0]*(x-[1])+[2]",0,2304);
  DirRMS->SetParameter(0,TMath::Tan(AngleDir));
  DirRMS->SetParameter(1,XIP);
  DirRMS->SetParameter(2,YIP);
  DirRMS->SetLineColor(kOrange);
  
  TCanvas* Canv = new TCanvas("Canv","Canv",1200,1200);
  Canv->cd();

  //Plot->Rebin2D(2,2);
  
  //Plot->GetXaxis()->SetRangeUser(XBar-6*(RangeX+RangeY)/2,XBar+6*(RangeX+RangeY)/2);
  //Plot->GetYaxis()->SetRangeUser(YBar-6*(RangeX+RangeY)/2,YBar+6*(RangeX+RangeY)/2);

  L->AddEntry(DirRMS,"Max RMS method");
  //L->AddEntry(DirGraph,"Graph with weights");
  L->AddEntry(MainAxis,"Track main axis");

  gStyle->SetOptStat(0000);
  
  Plot->Draw("COLZ");
  MainAxis->Draw("SAME");
  BarPlot->Draw("SAMEP");
  IPPlot->Draw("SAMEP");
  DirRMS->Draw("SAME");
  //DirGraph->Draw("SAME");
  arcMin->Draw("SAME");
  arcMax->Draw("SAME");
  L->Draw("SAME");
  Canv->SaveAs(Form("Out/TrackWithBarycenter_Ev%i_NTr%i.png",i,k));

  
  /*
  TFile* f = new TFile(Form("Out/File_Ev%i_Tr%i.root",i,k),"recreate");
  f->cd();
  Plot->Write();
  MainAxis->Write();
  BarPlot->Write();
  IPPlot->Write();
  DirRMS->Write();
  arcMin->Write();
  arcMax->Write();
  f->Save();
  f->Close();
  */

  
  delete Canv;
  delete BarPlot;
  delete MainAxis;
  delete arcMin;
  delete arcMax;
}

void MakeCanvasDist(TH2F* Plot,Float_t XBar,Float_t YBar,Float_t Phi,Int_t i,Int_t k,Float_t RMS,Float_t rmin, Float_t rmax,Float_t XIP,Float_t YIP,Float_t AngleDir,Float_t Rsel,Float_t XIPPrev,Float_t YIPPrev,Float_t ThTrue){

  Double_t RangeX = Plot->ProjectionX()->GetRMS();
  Double_t RangeY = Plot->ProjectionY()->GetRMS();
  
  TLegend* L = new TLegend(0.2,0.1,0.2,0.1);

  TGraph* BarPlot = new TGraph();
  BarPlot->SetName("BarPlot");
  BarPlot->SetTitle("BarPlot");
  
  BarPlot->SetPoint(0,XBar,YBar);
  BarPlot->SetMarkerStyle(8);
  BarPlot->SetMarkerSize(1.);

  TGraph* IPPlot = new TGraph();
  IPPlot->SetName("IPPlot");
  IPPlot->SetTitle("IPPlot");
  
  IPPlot->SetPoint(0,XIP,YIP);
  IPPlot->SetMarkerStyle(8);
  IPPlot->SetMarkerColor(kRed+2);
  IPPlot->SetMarkerSize(1.);

  TGraph* IPPreviusPlot = new TGraph();
  IPPreviusPlot->SetName("IPPreviusPlot");
  IPPreviusPlot->SetTitle("IPPreviusPlot");
  
  IPPreviusPlot->SetPoint(0,XIPPrev,YIPPrev);
  IPPreviusPlot->SetMarkerStyle(8);
  IPPreviusPlot->SetMarkerColor(kBlue+2);
  IPPreviusPlot->SetMarkerSize(1.);

  TF1* MainAxis = new TF1("MainAxis","[0]*(x-[1])+[2]",0,2304);
  MainAxis->SetParameter(0,TMath::Tan(Phi));
  MainAxis->SetParameter(1,XBar);
  MainAxis->SetParameter(2,YBar);

  TArc* arcMin = new TArc(XBar,YBar,rmin*sqrt(RMS));
  arcMin->SetFillColorAlpha(kBlue, 0);

  TArc* arcMax = new TArc(XBar,YBar,rmax*sqrt(RMS));
  arcMax->SetFillColorAlpha(kBlue, 0);

  TArc* arcSel = new TArc(XIP,YIP,Rsel);
  arcSel->SetFillColorAlpha(kBlue, 0);

  //DirGraph->SetLineColor(kGray);
      
  TF1* DirRMS = new TF1("DirRMS","[0]*(x-[1])+[2]",0,2304);
  DirRMS->SetParameter(0,TMath::Tan(AngleDir));
  DirRMS->SetParameter(1,XIP);
  DirRMS->SetParameter(2,YIP);
  DirRMS->SetLineColor(kOrange+1);
  
  TCanvas* Canv = new TCanvas("Canv","Canv",1200,1200);
  Canv->cd();

  //Plot->Rebin2D(2,2);

  //Plot->GetXaxis()->SetRangeUser(XBar-6*(RangeX+RangeY)/2,XBar+6*(RangeX+RangeY)/2);
  //Plot->GetYaxis()->SetRangeUser(YBar-6*(RangeX+RangeY)/2,YBar+6*(RangeX+RangeY)/2);
  
  //L->AddEntry(DirRMS,"Max RMS method");
  //L->AddEntry(DirGraph,"Graph with weights");
  //L->AddEntry(MainAxis,"Track main axis");

  L->AddEntry((TObject*)0, Form("%f",AngleDir/TMath::Pi()*180));
  L->AddEntry((TObject*)0, Form("%f",ThTrue));
  
  gStyle->SetOptStat(0000);
  
  Plot->Draw("COLZ");
  MainAxis->Draw("SAME");
  BarPlot->Draw("SAMEP");
  IPPlot->Draw("SAMEP");
  DirRMS->Draw("SAME");
  //DirGraph->Draw("SAME");
  arcMin->Draw("SAME");
  arcMax->Draw("SAME");
  arcSel->Draw("SAME");
  L->Draw("SAME");
  Canv->SaveAs(Form("Out/Track_Ev%i_NTr%i.png",i,k));

  TFile* f = new TFile(Form("Out/FileScaled_Ev%i_Tr%i.root",i,k),"recreate");
  f->cd();
  Plot->Write();
  MainAxis->Write();
  BarPlot->Write();
  IPPlot->Write();
  IPPreviusPlot->Write();
  DirRMS->Write();
  arcMin->Write();
  arcMax->Write();
  f->Save();
  f->Close();

  
  delete Canv;
  delete BarPlot;
  delete MainAxis;
  delete arcMin;
  delete arcMax;

}



///Slices Functions


Int_t GetIndexMinDist(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount){
  int index=0;
  Float_t dist=999;
  Float_t d;
  
  for(int i=0;i<nSl;i++){
    d=DistPit(x1,y1,X[parcount+i],Y[parcount+i]);
    if(d<dist){
      index=i;
      dist=d;
    }
  }

  return index;
}



Int_t GetIndexMinDistTemp(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl){
  int index=0;
  Float_t dist=999;
  Float_t d;
  
  for(int i=0;i<nSl;i++){
    d=DistPit(x1,y1,X[i],Y[i]);
    if(d<dist){
      index=i;
      dist=d;
    }
  }

  return index;
}


Int_t HasNeighbour(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount,Int_t radius){
  Int_t N=0;

  for(int i=0;i<nSl;i++){
    if(DistPit(x1,y1,X[parcount+i],Y[parcount+i])<radius) N++;
  }

  return N;
}

Int_t HasNeighbourTemp(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t radius){
  Int_t N=0;

  for(int i=0;i<nSl;i++){
    if(DistPit(x1,y1,X[i],Y[i])<radius) N++;
  }

  return N;
}

Float_t DistPit(Float_t x1,Float_t y1,Float_t x2,Float_t y2) {
  Float_t d;
  d=sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) );
  if(d>0){
    return d;
  } else {
    return 9999; ///////TIENILA POSITIVA, NON CAMBIARE MAI
  }
}

Int_t GetIndexMaxDist(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,Int_t parcount){
  int index=0;
  Float_t dist=0;
  Float_t d;
  
  for(int i=0;i<nSl;i++){
    d=DistPit(x1,y1,X[parcount+i],Y[parcount+i]);
    if(d>dist && X[parcount+i]>0){
      index=i;
      dist=d;
    }
    //    std::cout << i <<" \t d=" <<  d  << std::endl;
  }

  return index;
}


Int_t GetIndexMinChiOnLine(Float_t x1,Float_t y1,Float_t* X, Float_t* Y,Float_t nSl,TGraph* plot,TF1* line,Int_t parcount, Int_t radius){

  std::vector<int> index;

  for(int i=0;i<nSl;i++){
    if(DistPit(x1,y1,X[i],Y[i])< 50){
      index.push_back(i);
      //std::cout << "near index found:" << i+1 << std::endl;
    }
  }
  
  Int_t minDistIndex;
  Float_t Chi = 999999;
  Float_t NewChi;

  //line->FixParameter(0,line->GetParameter(0));
  //line->FixParameter(1,line->GetParameter(1));
  
  for(int i=0;i<index.size();i++){
    plot->SetPoint(3,X[index[i]],Y[index[i]]);
    plot->Fit(line,"QW");
    
    NewChi=line->GetChisquare();
    //std::cout << "indexChiCal:"<< index[i]+1 <<"\t chi2=" << NewChi << std::endl;
    if(NewChi<Chi){
      Chi=NewChi;
      minDistIndex=index[i];
    }
    plot->RemovePoint(3);
  }

  
  return minDistIndex;  
}



std::vector<int> SortSlices(Int_t nSc,Float_t* nSl, Float_t* X, Float_t* Y, Float_t* E){
  
  int parcount=0;
  int index,indexFirstNeigh;
  int nNeigh, nNeighFirstNeigh;
  Float_t xmax,ymax;

  std::vector<Int_t> OnMainTrackVec;
  
  Int_t r1=45;
  Int_t r2=55;
  
  bool Control=false;
  bool ControlOnTrack=false;
  
  for(int i=0;i<nSc;i++){
      
    std::vector<Float_t> Xtemp;
    std::vector<Float_t> Ytemp;
    std::vector<Float_t> Etemp;

    Float_t MeanX=0;
    Float_t MeanY=0;
    
    xmax=0;
    index=0;

    for(int j=0; j<nSl[i]; j++){
      Xtemp.push_back(X[parcount+j]);
      Ytemp.push_back(Y[parcount+j]);
      Etemp.push_back(E[parcount+j]);
            
      MeanX+=X[parcount+j];
      MeanY+=Y[parcount+j];
    }

    MeanX/=nSl[i];
    MeanY/=nSl[i];

  if(Control)  std:: cout << MeanX << "  " << MeanY <<std::endl;

    bool check=1;
    Int_t controlLoop=0;
    
    if(nSl[i]>2){
      do{
	index=GetIndexMaxDist(MeanX,MeanY,X,Y,nSl[i],parcount);    
	nNeigh=HasNeighbour(X[parcount+index],Y[parcount+index],X,Y,nSl[i],parcount,r1); //search first neighbour in radius r1
	if(Control) std::cout << "IndexMaxDistance=" << index << "\t nNeigh=" << nNeigh << std::endl; 
	if(nNeigh>0){
	  indexFirstNeigh=GetIndexMinDist(X[parcount+index],Y[parcount+index],X,Y,nSl[i],parcount);
	  
	  nNeighFirstNeigh=HasNeighbour(X[parcount+indexFirstNeigh],Y[parcount+indexFirstNeigh],X,Y,nSl[i],parcount,r2); //search neighbour first neighbour in radius r2
	  if(Control) std::cout << "indexFirstNeigh=" << indexFirstNeigh <<"\t nNeighFirstNeigh=" <<nNeighFirstNeigh <<std::endl;
	  if(nNeighFirstNeigh>1){
	    check=0;
	  } //chiudo if
	}//chiudo if
	
	X[parcount+index]=-999;
	Y[parcount+index]=-999;
	controlLoop++;
	if(controlLoop==10){
	  if(Control)std::cout << "ALLERT LOOP" << std::endl;
	  break;
	}
	
      }while(check);

    
      
      for(int j=0;j<nSl[i];j++){
	X[parcount+j]=Xtemp[j];
	Y[parcount+j]=Ytemp[j];
      }

      bool OnMainTrack=true;
      
      if(Control) std::cout << "Xseed=" << Xtemp[index] << "\t Yseed=" << Ytemp[index] <<"\t indexSeed="<< index+1 << std::endl;
      
      X[parcount+0]=Xtemp[index];
      Y[parcount+0]=Ytemp[index];
      E[parcount+0]=Etemp[index];
      
      Xtemp[index]=-9999;	    
      Ytemp[index]=-9999;	    
      Etemp[index]=-9999;

      OnMainTrackVec.push_back(0);

      if(ControlOnTrack)std::cout << "index=" << index << "\t X=" << X[parcount+0] << "\t Y=" << Y[parcount+0] << "\t control:" << OnMainTrack <<"\t IsOnTrack=" <<OnMainTrackVec[parcount+0]<< "\t dist=" << std::endl; 
      if(Control) std::cout << "nextindex" << index+1 << "\t will go in-> " << parcount+0 << std::endl;
      
      index=GetIndexMinDistTemp(X[parcount+0],Y[parcount+0],&Xtemp[0],&Ytemp[0],nSl[i]);
      
      X[parcount+1]=Xtemp[index];
      Y[parcount+1]=Ytemp[index];
      E[parcount+1]=Etemp[index];
      
      Xtemp[index]=-9999;	    
      Ytemp[index]=-9999;	    
      Etemp[index]=-9999;

      OnMainTrackVec.push_back(0);

      if(ControlOnTrack)std::cout << "index=" << index << "\t X=" << X[parcount+1] << "\t Y=" << Y[parcount+1] << "\t control:" << OnMainTrack <<"\t IsOnTrack=" <<OnMainTrackVec[parcount+1]<< "\t dist=" << std::endl;
      if(Control) std::cout << "nextindex" << index +1 << "\t will go in-> " << parcount+1 << std::endl;
      
      for(int j=2;j<nSl[i];j++){
		
	nNeigh=HasNeighbourTemp(X[parcount+j-1],Y[parcount+j-1],&Xtemp[0],&Ytemp[0],nSl[i],50); //BEST RAD 45 Mango 55 Lemon //search neighbour first neighbour to apply the line fit condition for next point
	if(Control)std::cout << "nNeigh"<< nNeigh << std::endl; 
	if(nNeigh<=1){

	  index=GetIndexMinDistTemp(X[parcount+j-1],Y[parcount+j-1],&Xtemp[0],&Ytemp[0],nSl[i]);

	} else if(nNeigh>1){

	  TGraph* Graph2pt = new TGraph();
	  Graph2pt->SetPoint(0,X[parcount+j-1],Y[parcount+j-1]);
	  Graph2pt->SetPoint(1,X[parcount+j-2],Y[parcount+j-2]);
	  Graph2pt->SetPoint(2,X[parcount+j-3],Y[parcount+j-3]);
	  //Graph2pt->SetPoint(3,X[parcount+j-4],Y[parcount+j-4]);
	  	  
	  TF1* Fit2pt = new TF1("Fit2pt","pol1");
	  Graph2pt->Fit(Fit2pt,"Q");
	  index=GetIndexMinChiOnLine(X[parcount+j-1],Y[parcount+j-1],&Xtemp[0],&Ytemp[0],nSl[i],Graph2pt,Fit2pt,parcount,40); //radius of search for the neighbour to fit with
	  
	}//chiudo else

	if(Control) std::cout << "nextindex=" << index+1 << "\t will go in-> " << parcount+j << std::endl;
	
	X[parcount+j]=Xtemp[index];
	Y[parcount+j]=Ytemp[index];
	E[parcount+j]=Etemp[index]; 
	
	Xtemp[index]=-9999;
	Ytemp[index]=-9999;
	Etemp[index]=-9999;
	
	
	if(DistPit(X[parcount+j],Y[parcount+j],X[parcount+j-1],Y[parcount+j-1])>54){ // radius to exclude outcomes
	  OnMainTrack=false;
	}//chiudo if
	
	if(OnMainTrack==true){
	  OnMainTrackVec.push_back(0);
	  if(ControlOnTrack) std::cout << "OnTrack:" << OnMainTrack <<"\t messo" <<OnMainTrackVec.size()<< std::endl; 
	} else {
	  OnMainTrackVec.push_back(1);
	  if(ControlOnTrack) std::cout << "NotOnTrack:" << OnMainTrack<<"\t messo" <<OnMainTrackVec.size() << std::endl;
	}//chiudo else
	
	if(ControlOnTrack) std::cout <<"nSc="<<i <<"\t index=" << parcount+j << "\t X=" << X[parcount+j] << "\t Y=" << Y[parcount+j] << "\t control:" << OnMainTrack <<"\t IsOnTrack=" <<OnMainTrackVec[parcount+j]<< "\t dist="<<DistPit(X[parcount+j],Y[parcount+j],X[parcount+j-1],Y[parcount+j-1]) << std::endl;
      }//chiudi for
      
    }//chiudo if nSl
    else { for(int k=0;k<nSl[i];k++){OnMainTrackVec.push_back(0);}}

    //for(int k=0;k<OnMainTrackVec.size();k++) std::cout << OnMainTrackVec[k] << std::endl;
    
    if(Control) std::cout << "\n" << std::endl;
    parcount+=nSl[i];
    
  }//chiudo for nSC
    
  return OnMainTrackVec;
}//chiudo funzione
