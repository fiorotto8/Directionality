#include <cstdio>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <Riostream.h>
#include <TFile.h>
#include "Analyzer.h"
#include "TAxis.h"
#include <TH1.h>
#include <TF1.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TSpectrum.h>
#include <TMatrixD.h>
#include "TDecompSVD.h"


double RMSOnLine(double XBar, double YBar, double Phi);

//default constructor
Analyzer::Analyzer():
fminx(0.),
fminy(0.),
fmaxx(0.),
fmaxy(0.),
fintegral(0.),
fradius(0.),
fheight(0.),
fxcentr(0.),
fycentr(0.),
fTrack(nullptr),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fnpixelx(0),
fnpixely(0),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fPhiMainAxis(0.),
fLineMaxRMS(nullptr),
fRMSOnMainAxis(0.),
fSkewOnLine(0.),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr), 
fPhiDir(0.),
fLineDirection(nullptr)
{
}

//Standard copyconstructor
Analyzer::Analyzer(const Analyzer& source):
fminx(source.fminx),
fminy(source.fminy),
fmaxx(source.fmaxx),
fmaxy(source.fmaxy),
fintegral(source.fintegral),
fradius(source.fradius),
fheight(source.fheight),
fxcentr(source.fxcentr),
fycentr(source.fycentr),
fnpixelx(source.fnpixelx),
fnpixely(source.fnpixely),
fNPIP(source.fNPIP),
fwScal(source.fwScal),
fXbar(source.fXbar),
fYbar(source.fYbar),
fPhiMainAxis(source.fPhiMainAxis),
fRMSOnMainAxis(source.fRMSOnMainAxis),
fSkewOnLine(source.fSkewOnLine),
fXIPPrev(source.fXIPPrev),
fYIPPrev(source.fYIPPrev),
fXIP(source.fXIP),
fYIP(source.fYIP),
fPhiDir(source.fPhiDir)
{
  if(source.fTrack==nullptr)	   fTrack=nullptr;
  else fTrack=(TH2F*)source.fTrack->Clone("Trackcopy");
  if(source.fTrackTail==nullptr)	   fTrackTail=nullptr;
  else fTrackTail=(TH2F*)source.fTrackTail->Clone("TrackTailcopy");
  if(source.fScaledTrack==nullptr)	   fScaledTrack=nullptr;
  else fScaledTrack=(TH2F*)source.fScaledTrack->Clone("TrackScaledcopy");
  if(source.fBarPlot==nullptr)	   fBarPlot=nullptr;
  else fBarPlot= new TGraph(*source.fBarPlot);
  if(source.fIPPlot==nullptr)	   fIPPlot=nullptr;
  else fIPPlot=new TGraph(*source.fIPPlot);
  if(source.fLineMaxRMS==nullptr)	   fLineMaxRMS=nullptr;
  else fLineMaxRMS=(TF1*)source.fLineMaxRMS->Clone("LineMaxRMScopy");
  if(source.fLineDirection==nullptr)	   fLineDirection=nullptr;
  else fLineDirection=(TF1*)source.fLineDirection->Clone("LineDirectioncopy");

}
////constructor: given a 2D histogram generate the object with the same track but resized
//// npixel= number of pixel of the new track, Tracklarge= original 2D histo with the track, npixelorig= number of pixels of the original 2D histogram
////integral not evaluated in this constructor
Analyzer::Analyzer(const char* nometh2, int npixel, TH2F *Tracklarge, int npixelorig):
fintegral(0.),
fradius(0.),
fheight(0.),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fPhiMainAxis(0.),
fLineMaxRMS(nullptr),
fRMSOnMainAxis(0.),
fSkewOnLine(0.),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr), 
fPhiDir(0.),
fLineDirection(nullptr)
{
	fycentr=Tracklarge->GetMaximumBin()/(npixelorig+2);
	fxcentr=Tracklarge->GetMaximumBin()%(npixelorig+2);

	fminx=npixelorig+100;
	fminy=npixelorig+100;
	fmaxx=-1;
	fmaxy=-1;
	for(int i=fxcentr-npixel/2;i<fxcentr+npixel/2;i++)
	{
		for(int j=fycentr-npixel/2;j<fycentr+npixel/2;j++)
		{
			double x=Tracklarge->GetBinContent(i,j);
			if(i<fminx && x>0)  fminx=i;
			if(i>fmaxx && x>0)  fmaxx=i;
			if(j<fminy && x>0)  fminy=j;
			if(j>fmaxy && x>0)  fmaxy=j;
		}
	}

	fnpixelx=fmaxx-fminx+1+10;
	fnpixely=fmaxy-fminy+1+10;
	fTrack=new TH2F(nometh2,nometh2,fnpixelx,0,fnpixelx,fnpixely,0,fnpixely);
	for(int j=1;j<=fnpixelx;j++)
	{
		for(int k=1;k<=fnpixely;k++)
		{
			fTrack->SetBinContent(j,k,Tracklarge->GetBinContent(fminx+j-1-5,fminy+k-1-5));
		}
	}
	

	TH1D *tax=fTrack->ProjectionX();
	tax->Fit("gaus","Q");
	fxcentr=tax->GetFunction("gaus")->GetParameter(1);
	TH1D *tay=fTrack->ProjectionY();
	tay->Fit("gaus","Q");
	fycentr=tay->GetFunction("gaus")->GetParameter(1);
	delete tax;
	delete tay;
	
	for(int j=1;j<=fnpixelx;j++)
	{
		for(int k=1;k<=fnpixely;k++)
		{
			if(fTrack->GetBinContent(j,k)>0)
			{
				double x=sqrt((fycentr-k)*(fycentr-k) + (fxcentr-j)*(fxcentr-j) );
				if(x>fradius)  fradius=x;
			}
		}
	}
	
}

////////////////////This constructor may not work
Analyzer::Analyzer(const char* nometh2, TH2F *Tracklarge):
fradius(0.),
fheight(0.),
fxcentr(0),
fycentr(0),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr), 
fPhiDir(0.),
fLineDirection(nullptr)
{
	fintegral=0.;
	
	fminx=3000;
	fminy=3000;
	fmaxx=0;
	fmaxy=0;

	for(int i=0;i<Tracklarge->GetXaxis()->GetNbins();i++){
	  for(int j=0;j<Tracklarge->GetYaxis()->GetNbins();j++){

	    double x=Tracklarge->GetBinContent(i,j);
	    
	    if(i<fminx && x>0)  fminx=i;
	    if(i>fmaxx && x>0)  fmaxx=i;
	    if(j<fminy && x>0)  fminy=j;
	    if(j>fmaxy && x>0)  fmaxy=j;
	    
	  }
	}

	fmaxx=fmaxx+30;
	fminx=fminx-30;
	fmaxy=fmaxy+30;
	fminy=fminy-30;
	
	fnpixelx=fmaxx-fminx;
	fnpixely=fmaxy-fminy;

	//std::cout << fnpixelx<<"\t"<<fminx <<"\t"<< fmaxx << "\t" << fnpixely <<"\t"<< fminy <<"\t"<< fmaxy << std::endl;
	
	fTrack=new TH2F(Form("A%s",nometh2),Form("A%s",nometh2),fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
	int np=0;
	  
	for(int j=0;j<=fnpixelx;j++)
	{
		for(int k=0;k<=fnpixely;k++)
		{
		  if(Tracklarge->GetBinContent(fminx+j,fminy+k)>0){
		    fTrack->SetBinContent(j,k,Tracklarge->GetBinContent(fminx+j,fminy+k));
		    fintegral+=Tracklarge->GetBinContent(fminx+j,fminy+k);
		  }
		}
	}
	
		
	fPhiMainAxis=AngleLineMaxRMS();
	BuildLineMaxRMS();		//defines fLineMaxRMS
	fRMSOnMainAxis=GetRMSOnMainAxis();
	fSkewOnLine=SkewOnMainAxis();
	
	
					
	TCanvas* OriTrack = new TCanvas();
	Tracklarge->GetXaxis()->SetRangeUser(fminx,fmaxx);
	Tracklarge->GetYaxis()->SetRangeUser(fminy,fmaxy);
	Tracklarge->Draw("COLZ");
	OriTrack->SaveAs(Form("Tracks/%s.png",Tracklarge->GetName()));
	delete OriTrack;

	fTrack->Rebin2D(2,2);
}

////Constructor 
//////Constructor taking the arrays of x,y,z coord of the track. B and E are delimeters for the track hits
////fTrack is rebbined
Analyzer::Analyzer(const char* nometh2,float* X,float* Y,float* Z,int B,int E):
fradius(0.),
fheight(0.),
fxcentr(0),
fycentr(0),
fTrackTail(nullptr),
fScaledTrack(nullptr),
fNPIP(0),
fwScal(0.),
fXbar(0.),
fYbar(0.),
fBarPlot(nullptr),
fXIPPrev(0.),
fYIPPrev(0.),
fXIP(0.),
fYIP(0.),
fIPPlot(nullptr), 
fPhiDir(0.),
fLineDirection(nullptr)
{
	fintegral=0.;
	
	fminx = TMath::MinElement(E-B,&X[B]);
	fmaxx = TMath::MaxElement(E-B,&X[B]);
	fminy = TMath::MinElement(E-B,&Y[B]);
	fmaxy = TMath::MaxElement(E-B,&Y[B]);
	
	
	fmaxx=fmaxx+30;
	fminx=fminx-30;
	fmaxy=fmaxy+30;
	fminy=fminy-30;
	
	fnpixelx=fmaxx-fminx;
	fnpixely=fmaxy-fminy;

	//std::cout << fnpixelx<<"\t"<<fminx <<"\t"<< fmaxx << "\t" << fnpixely <<"\t"<< fminy <<"\t"<< fmaxy << std::endl;
	
	fTrack=new TH2F(Form("A%s",nometh2),Form("A%s",nometh2),fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
		  
	for(int i=B;i<E;i++){
	  fTrack->SetBinContent(fTrack->GetXaxis()->FindBin(X[i]),fTrack->GetYaxis()->FindBin(Y[i]),Z[i]);
	  fintegral+=Z[i];
	}

	//fTrack->Rebin2D(2,2);
	
	fPhiMainAxis=AngleLineMaxRMS();
	BuildLineMaxRMS();
	fRMSOnMainAxis=GetRMSOnMainAxis();
	fSkewOnLine=SkewOnMainAxis();	
}

//Destructor
Analyzer::~Analyzer()
{
	if(fTrack!=nullptr)	 	 	 delete fTrack;
	if(fTrackTail!=nullptr)	 	 delete fTrackTail;
	if(fScaledTrack!=nullptr)	 delete fScaledTrack;
	if(fBarPlot!=nullptr)	 	 delete fBarPlot;
	if(fLineMaxRMS!=nullptr)	 delete fLineMaxRMS;
	if(fIPPlot!=nullptr)	 	 delete fIPPlot;
	if(fLineDirection!=nullptr)	 delete fLineDirection;
}


///////////////Other functions
//Initializer for fLineMaxRMS
void Analyzer::BuildLineMaxRMS(){

  fLineMaxRMS= new TF1("LineMaxRMS","[0]*(x-[1])+[2]",0,2304);

  fLineMaxRMS->SetParameter(1,fXbar);
  fLineMaxRMS->SetParameter(2,fYbar);
  fLineMaxRMS->SetParameter(0,TMath::Tan(fPhiMainAxis));
  
}

//Initializer for fLineDirection
void Analyzer::BuildLineDirection(){

  fLineDirection= new TF1("LineDirection","[0]*(x-[1])+[2]",0,2304);

  fLineDirection->SetParameter(1,fXIP);
  fLineDirection->SetParameter(2,fYIP);
  fLineDirection->SetParameter(0,TMath::Tan(fPhiDir));
  
}

//Resets the variables
void Analyzer::Reset()
{
	fminx=0;
	fminy=0;
	fmaxx=0;
	fmaxy=0;
	fintegral=0;
	fheight=0;
	fradius=0;
	fxcentr=0;
	fycentr=0;
	fnpixelx=0;
	fnpixely=0;
	fNPIP=0;
	fwScal=0;
	fXbar=0;
    fYbar=0;
    fPhiMainAxis=0;
    fRMSOnMainAxis=0;
    fSkewOnLine=0;
    fXIPPrev=0;
    fYIPPrev=0;
    fXIP=0;
    fYIP=0;
    fPhiDir=0;
    
	
	delete fTrack;
	fTrack=nullptr;
    delete fTrackTail;
    fTrackTail=nullptr;
    delete fScaledTrack;
	fScaledTrack=nullptr;
	delete fLineDirection;
	fLineDirection=nullptr;
	delete fIPPlot;
	fIPPlot=nullptr;
	delete fLineMaxRMS;
	fLineMaxRMS=nullptr;
	delete fBarPlot;
	fBarPlot=nullptr;
}

//Calculation of the integral
void Analyzer::Integral()
{
	fintegral=0.;
	for(int i=fminx;i<fmaxx;i++)
	{
		for(int j=fminy;j<fmaxy;j++)
		{
			fintegral+=fTrack->GetBinContent(i,j);
		}
	}
}

//Saves histogram to root file
void Analyzer::SavetoFile(const char* nometh2)
{
	fTrack->SetName(nometh2);
	fTrack->Write();
	return;
}

//SavePic
void Analyzer::SavePic(const char* nomepic)
{
  TCanvas* canv = new TCanvas("canv","canv",1500,1500);
  
  TGraph* g = new TGraph();
  g->SetPoint(0,fXbar,fYbar);
  g->SetMarkerStyle(8);

  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0,Form("NPx=%i",fnpixelx*fnpixely)); 	//assuming the track is a rectangle... doubtful
  l->AddEntry((TObject*)0,Form("Int=%f",fintegral));
  l->AddEntry((TObject*)0,Form("Dens=%f",fintegral/(fnpixelx*fnpixely)));
  l->AddEntry((TObject*)0,Form("Skew=%f",fSkewOnLine));
  l->AddEntry((TObject*)0,Form("SkewNorm=%f",fSkewOnLine/fintegral));
  
  l->AddEntry((TObject*)0,Form("RMS=%f",fRMSOnMainAxis));
  l->AddEntry((TObject*)0,Form("RMSNorm=%f",fRMSOnMainAxis/fintegral));
  
  //fTrack->SetName(nomepic);
  fTrack->Draw("COLZ");
  g->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");
  l->Draw("SAME");
  
  canv->SaveAs(Form("Tracks/%s",nomepic));

  delete canv;
  
  return;
}

void Analyzer::SavePicDir(const char* nomepic){

  TCanvas* canv = new TCanvas("canv","canv",3500,1500);
  canv->Divide(3,1);
  
  
  
  TLegend* l = new TLegend();
  l->AddEntry((TObject*)0, Form("%f",fPhiDir/TMath::Pi()*180));
  
  canv->cd(1);
  fTrack->Draw("COLZ");
  fBarPlot->Draw("SAMEP");
  fLineMaxRMS->Draw("SAME");

  canv->cd(2);
  fTrackTail->Draw("COLZ");
  fIPPlot->Draw("SAMEP");

  canv->cd(3);
  fScaledTrack->Draw("COLZ");
  fIPPlot->Draw("SAMEP");
  fLineDirection->Draw("SAME");
  l->Draw("SAME");

  canv->SaveAs(Form("Tracks/%s",nomepic));

  delete canv;
}

//Saves the output root file
void Analyzer::SaveRootFile(const char* nomefile){
  
  TFile* f = new TFile(Form("Tracks/%s",nomefile),"recreate");
  f->cd();

  fTrack->Write();
  fTrackTail->Write();
  fScaledTrack->Write();
  fBarPlot->Write();
  fIPPlot->Write();
  fLineDirection->Write();

  f->Save();
  f->Close();
}


void Analyzer::Barycenter()
{
  double Xb=0;
  double Yb=0;
  double Z=0;
  double ChargeTot=0;
  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z>0)
	  {
	    Xb+=(Z*fTrack->GetXaxis()->GetBinCenter(i));
	    Yb+=(Z*fTrack->GetYaxis()->GetBinCenter(j));
	    ChargeTot+=Z;
	  }
	}
  }
  
  Xb/=ChargeTot;
  Yb/=ChargeTot;

  fXbar=Xb;
  fYbar=Yb;

  return;
}

//Calculates the barycentre of the track
void Analyzer::Barycenter(TH2F* Tr,double *X,double *Y)
{
  double Xb=0;
  double Yb=0;
  double Z=0;
  double ChargeTot=0;
  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=Tr->GetBinContent(i,j);
	  if(Z>0)
	  {
	    Xb+=(Z*Tr->GetXaxis()->GetBinCenter(i));
	    Yb+=(Z*Tr->GetYaxis()->GetBinCenter(j));
	    ChargeTot+=Z;
	  }
	}
  }
  
  Xb/=ChargeTot;
  Yb/=ChargeTot;

  *X=Xb;
  *Y=Yb;

  return;
}

//Find main line (maximize RMS on this line)
double Analyzer::AngleLineMaxRMS()			
{
  double Sum1=0;
  double Sum2=0;
  double Z=0.;
  double Phi;
  double RmsAng;
  double RmsAngPerp;
  
  Barycenter(fTrack,&fXbar,&fYbar);

  fBarPlot=new TGraph();
  fBarPlot->SetName("BarycenterPlot");
  fBarPlot->SetPoint(0,fXbar,fYbar);
  fBarPlot->SetMarkerStyle(8);

  
  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z>0)
	  {
	    Sum1+= Z*(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*(fTrack->GetYaxis()->GetBinCenter(j)-fYbar);
	    Sum2+= Z*( (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*(fTrack->GetYaxis()->GetBinCenter(j)-fYbar) - (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)  );
	  }
	}
  }

  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);
  
  RmsAng=RMSOnLine(Phi);
  RmsAngPerp=RMSOnLine(Phi+TMath::Pi()/2);
  
  if( RmsAng > RmsAngPerp )
  {  
      fRMSOnMainAxis=RmsAng;
      return Phi;
  } 
  else 
  {
      fRMSOnMainAxis=RmsAngPerp;
      if(Phi+TMath::Pi()/2>TMath::Pi()/2)     return Phi+TMath::Pi()/2-TMath::Pi();
      else     return Phi+TMath::Pi()/2;
  }
  
}


//Called by AngleLineMaxRMS
double Analyzer::RMSOnLine(double Phi)
{
  double RMS=0;
  double ChargeTot=0;
  double Z=0.;


  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z!=0)
	  {
	    RMS+= Z*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(Phi) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(Phi) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(Phi) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(Phi) );
	    ChargeTot+=Z;
	  }
	}
  }
  
  return RMS;
}





//
double Analyzer::SkewOnMainAxis(){

  Float_t Skew=0;
  Float_t ChargeTot=0;
  Float_t Z;
  

  for(int i=1;i<fnpixelx;i++)
    {
      for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z>0)
	    {
	      Skew+= Z*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) )*( (fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*cos(fPhiMainAxis) + (fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*sin(fPhiMainAxis) );
	      ChargeTot+=Z;
	    }
	}
    }
  
  fSkewOnLine=Skew;
  
  return Skew;
  
}

//
void Analyzer::RemoveNoise(){
  
  std::vector<Int_t> ToRemoveXBin;
  std::vector<Int_t> ToRemoveYBin;

  Float_t ETemp;

  for(int i=1; i<fnpixelx-1;i++){
    for(int j=1;j<fnpixely-1;j++){
      if(fTrack->GetBinContent(i,j)>0){
	ETemp=0;
	if(fTrack->GetBinContent(i+1,j)>=0) ETemp+=fTrack->GetBinContent(i+1,j);
	if(fTrack->GetBinContent(i-1,j)>=0) ETemp+=fTrack->GetBinContent(i-1,j);
	if(fTrack->GetBinContent(i,j+1)>=0) ETemp+=fTrack->GetBinContent(i,j+1);
	if(fTrack->GetBinContent(i,j-1)>=0) ETemp+=fTrack->GetBinContent(i,j-1);
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
    fTrack->SetBinContent(ToRemoveXBin[i],ToRemoveYBin[i],0);
  }

}

//
void Analyzer::ApplyThr(){
  Int_t XBinMin=fTrack->GetXaxis()->GetFirst();
  Int_t XBinMax=XBinMin+fTrack->GetXaxis()->GetNbins();

  Int_t YBinMin=fTrack->GetYaxis()->GetFirst();
  Int_t YBinMax=YBinMin+fTrack->GetYaxis()->GetNbins();
  
  Float_t z;
  
  for(int i=XBinMin; i<XBinMax;i++){
    for(int j=YBinMin;j<YBinMax;j++){
      z=fTrack->GetBinContent(i,j);
      if(z>0 && z<0.5){
	fTrack->SetBinContent(i,j,0);
      }
    }
  }
 
}
//
void Analyzer::ImpactPoint(const char* nometh2){
  
  fTrackTail = new TH2F(nometh2,nometh2,fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);

  fTrackTail->Rebin2D(2,2);
  
  Float_t rminN=0.7;
  Double_t X,Y,Z;
  Int_t NSelPoints;
  Float_t PointSkew,PointDistCm;

  std::vector<float> XIntPointPrev;
  std::vector<float> YIntPointPrev;
  
  
  do{
    rminN+=0.03;
    NSelPoints=0;
    
    fTrackTail->Reset();
    
    for(int j=0;j<fTrack->GetXaxis()->GetNbins();j++){
      for(int l=0;l<fTrack->GetYaxis()->GetNbins();l++){
	
	X=fTrack->GetXaxis()->GetBinCenter(j);
	Y=fTrack->GetYaxis()->GetBinCenter(l);
	Z=fTrack->GetBinContent(j,l);
	
	PointSkew=GetPointSkew(X,Y);
	PointDistCm=PDistCm(X,Y);
	
	if(PointSkew>0 && PointDistCm> rminN && Z>0){
	  fTrackTail->SetBinContent(fTrackTail->GetXaxis()->FindBin(X),fTrackTail->GetYaxis()->FindBin(Y),Z);
	  NSelPoints++;
	}//chiudo if selection
	
      }//chiuso for l
    }//chiudo for j (fill histos)

    Barycenter(fTrackTail,&fXIPPrev,&fYIPPrev);

    XIntPointPrev.push_back(fXIPPrev);
    YIntPointPrev.push_back(fYIPPrev);
    
  }while(NSelPoints>fNPIP);
  
  Barycenter(fTrackTail,&fXIP,&fYIP);
  fXIPPrev=XIntPointPrev[(int)(XIntPointPrev.size()/2)];
  fYIPPrev=YIntPointPrev[(int)(YIntPointPrev.size()/2)];

  fIPPlot = new TGraph();
  fIPPlot->SetName("IPPLot");
  fIPPlot->SetPoint(0,fXIP,fYIP);
  fIPPlot->SetMarkerStyle(8);
  
}

//
void Analyzer::ScaledTrack(const char* nometh2){

  fScaledTrack = new TH2F(nometh2,nometh2,fnpixelx,fminx,fmaxx,fnpixely,fminy,fmaxy);
  fScaledTrack->Rebin2D(2,2);
  
  double X,Y,Z;

  for(int j=0;j<fnpixelx;j++){
    for(int l=0;l<fnpixely;l++){
      if(fTrack->GetBinContent(j,l)>0){
	X=fTrack->GetXaxis()->GetBinCenter(j);
	Y=fTrack->GetYaxis()->GetBinCenter(l);
	Z=fTrack->GetBinContent(j,l);

	fScaledTrack->SetBinContent(fScaledTrack->GetXaxis()->FindBin(X),fScaledTrack->GetYaxis()->FindBin(Y),Z*exp(-( sqrt( (X-fXIP)*(X-fXIP)+(Y-fYIP)*(Y-fYIP) ) )/fwScal ) );
      }
    }//chiudo l
  }//chiudo for j
  
  
}

//
void Analyzer::Direction()			
{
  double Sum1=0;
  double Sum2=0;
  double Z=0.;
  double Phi;
  double RmsAng;
  double RmsAngPerp;
  
  
  for(int i=1;i<fnpixelx;i++)
    {
      for(int j=1;j<fnpixely;j++)
	{  
	  Z=fScaledTrack->GetBinContent(i,j);
	  if(Z>0)
	    {
	      Sum1+= Z*(fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)*(fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP);
	      Sum2+= Z*( (fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP)*(fScaledTrack->GetYaxis()->GetBinCenter(j)-fYIP) - (fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)*(fScaledTrack->GetXaxis()->GetBinCenter(i)-fXIP)  );
	    }
	}
    }
  
  Phi=-0.5*TMath::ATan(2*Sum1/Sum2);
  
  RmsAng=RMSOnLine(Phi);
  RmsAngPerp=RMSOnLine(Phi+TMath::Pi()/2);
  
  if( RmsAng > RmsAngPerp )
    {  
      fPhiDir=Phi;
    } 
  else 
    {
      fRMSOnMainAxis=RmsAngPerp;
      if(Phi+TMath::Pi()/2>TMath::Pi()/2)     fPhiDir = Phi+TMath::Pi()/2-TMath::Pi();
      else     fPhiDir= Phi+TMath::Pi()/2;
    }
    
}

//
void Analyzer::ImprCorrectAngle(){

  Float_t AnglePerp = fPhiDir+TMath::Pi()/2;
  
  Float_t qIP= fYIP-TMath::Tan(AnglePerp)*fXIP;   
  Float_t qPIP= fYIPPrev-TMath::Tan(AnglePerp)*fXIPPrev;

  if(fPhiDir>0){
    if(qIP<qPIP) {
      return;
    } else {
      fPhiDir= fPhiDir-TMath::Pi();
    }//chiudo else qIP
  } else {
    if(qIP<qPIP){
      fPhiDir= fPhiDir+TMath::Pi();
    } else {
      return;
    }//chiudo else qIP
  }//chiudo esle angolo
  
}

//Leftmost and rightmost points in track (used in longitudinal and transverse profile) 
void Analyzer::Edges(double &Xl, double &Yl, double &Xr, double &Yr, double slope) 
{
  double Xp, Yp, Zp;
  double dist, tempdist_r=0, tempdist_l=0;
  double ii=0, jj=0;

  Barycenter(fTrack, &fXbar, &fYbar);

  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  
	  Zp=fTrack->GetBinContent(i,j);
	  if(Zp!=0)
	  {
	    ii=fTrack->GetXaxis()->GetBinCenter(i);
	    jj=fTrack->GetYaxis()->GetBinCenter(j);
	    Xp = (1./(1+pow(slope,2)))*(ii+fXbar*pow(slope,2)+slope*(jj-fYbar));
	    Yp = fYbar+(slope/(1+pow(slope,2)))*(ii-fXbar+slope*(jj-fYbar));
	    dist = sqrt(pow((Xp-fXbar),2)+pow((Yp-fYbar),2));
	    if(dist>tempdist_r && Xp>fXbar)
	    {
			tempdist_r = dist;
			Xr = Xp; Yr = Yp;
	    }
	    else if(dist>tempdist_l && Xp<fXbar)
	    {
			tempdist_l = dist;
			Xl = Xp; Yl = Yp;
	    }
	  }

	}
  }

  return;
}


//Profiling over main axis of the track
//if longitudinal is 1 it profiles on longitudinal axis (along principal axis), else it profiles on transverse axis (perpendicular) 
TH1D* Analyzer::FillProfile(bool longitudinal, float x1, float x2)
{

  double xl,yl,xr,yr;
  double slope;
  double ii=0, jj=0;
  
  if(longitudinal) slope =  tan(AngleLineMaxRMS());
  else slope = tan(TMath::Pi()/2.+AngleLineMaxRMS());
  
  Edges(xl,yl,xr,yr,slope);
  int binmax = (int)sqrt(pow((xl-xr),2)+pow((yl-yr),2));
  TH1D* TrackProfile=new TH1D("TrackProf","TrackProf",binmax+2,0,binmax+2);
  
  double Xp, Yp, Zp;

  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  
	  Zp=fTrack->GetBinContent(i,j);
	  if(Zp!=0 )
	  {
	    ii=fTrack->GetXaxis()->GetBinCenter(i);
	    jj=fTrack->GetYaxis()->GetBinCenter(j);
            if(ii>x1 && ii<x2){
	      Xp = (1/(1+pow(slope,2)))*(ii+fXbar*pow(slope,2)+slope*(jj-fYbar));
	      Yp = fYbar+(slope/(1+pow(slope,2)))*(ii-fXbar+slope*(jj-fYbar));
  	      TrackProfile->Fill(sqrt(pow((Xp-xl),2)+pow((Yp-yl),2)),Zp);
            }
      }
    }
  }

  return TrackProfile;
}

//To be called with longitudinal or transverse profile to cut it around the main energy deposit
TH1D* Analyzer::CutProfile(TH1D* profile, double height)
{ 
//height: cut bins with content lower than height*maximum intensity in the profile  (default is 0.25%)

  int binmin = profile->FindFirstBinAbove(height*profile->GetMaximum(),1);
  int binmax = profile->FindLastBinAbove(height*profile->GetMaximum(),1); 
  int bins_cut = 0;

  std::string title = profile->GetTitle();
  title += "_cut";
  TH1D* profile_cut = new TH1D(title.c_str(),title.c_str(),binmax-binmin,0,binmax-binmin);

  for(int bins=binmin; bins<binmax; bins++){
    profile_cut->SetBinContent(bins_cut,profile->GetBinContent(bins));
    bins_cut++;
  }

  return profile_cut;

}

//Profile along x direction
TH1D* Analyzer::FillProfileX() 
{
  double Zp;
  
  TH1D* TrackProfileX=new TH1D("TrackProfX","TrackProfX",fnpixelx,fminx,fmaxx);

  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  
	  Zp=fTrack->GetBinContent(i,j);
	  if(Zp!=0)
	  {
  	    TrackProfileX->Fill(fTrack->GetXaxis()->GetBinCenter(i),Zp);
      }
    }
  }

  return TrackProfileX;

}

//Profile along y direction
TH1D* Analyzer::FillProfileY() 
{
  double Zp;

  TH1D* TrackProfileY=new TH1D("TrackProfY","TrackProfY",fnpixely,fminy,fmaxy);

  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{ 
	  Zp=fTrack->GetBinContent(i,j);
	  if(Zp!=0)
	  {
  	    TrackProfileY->Fill(fTrack->GetYaxis()->GetBinCenter(j),Zp);
      }
     }
  }

  return TrackProfileY;

}


void Analyzer::AnglePCA(double &ang)
{
TMatrixD M(2,2);
TMatrixD V(2,2);
double x,y,x2,y2,xy;
double Z;

  for(int i=0; i<fnpixelx; i++)
  {
    for(int j=0; j<fnpixely; j++)
    {
      Z=fTrack->GetBinContent(i,j);
      if(Z!=0){
        x=(fTrack->GetXaxis()->GetBinCenter(i)-fXbar)*Z;
        y=(fTrack->GetYaxis()->GetBinCenter(j)-fYbar)*Z;
        x2 += x*x;
        y2 += y*y;
        xy += x*y;
      }
    }
  }

M(0,0) = x2;
M(0,1) = M(1,0) = xy;
M(1,1) = y2;

TDecompSVD SVDmatrix(M);
if(SVDmatrix.Decompose()) V = SVDmatrix.GetV();

ang = TMath::ATan(V(1,0)/V(0,0));

}

void Analyzer::FindNPeaks(TH1D* h, std::vector<std::pair<double,double>> &foundpeaks)
{

  TSpectrum* s = new TSpectrum();

  int npeaks;
  double* peaksPos(nullptr);

  std::vector< std::pair<double,double> > peaks;
  std::vector<double> peaks_tocompare;

  for(int i=2; i<16; i=i+1){ //scan for different sigma

    npeaks = s->Search(h,i,"nobackground",0.1); //number of peaks with current sigma (npeaks2=number of peaks with previous sigma)
    peaksPos = s->GetPositionX(); //positions of peaks with current sigma (peaks2=positions of peaks with previous sigma)
    if (npeaks == 0){continue;} //if no peaks are found, go to next sigma
    else{ //if peaks are found
      //if any of the peaks found in the previous iteration is equal (within one sigma) to the current one, ignore it; otherwise, save the position of the new peak and increase the total number of peaks
      for(int j=0; j<npeaks; j++){ //loop over number of new peaks

        if(!peaks_tocompare.empty()){ 

          auto it = std::find_if(peaks_tocompare.begin(), peaks_tocompare.end(), [&](double p){ return (p>(peaksPos[j]-(double)i) && p<(peaksPos[j]+(double)i)); });
          if(it != peaks_tocompare.end()){ //if it's already stored
            continue;
          }
	  else{ //it's a new peak, save it
            peaks.push_back( std::make_pair(peaksPos[j],(double)i) );
	    continue;
          }
        }
        peaks.push_back(std::make_pair(peaksPos[j],(double)i));
      } //end loop sui picchi trovati
    } //end if peaks were found
  
  peaks_tocompare.clear();
  for(int k=0; k<npeaks; k++){peaks_tocompare.push_back(peaksPos[k]);}
  
  } //end scan on different sigma

foundpeaks = peaks; 

 delete s;
 
 return;

}

//coordinates of maximum intensity pixel
void Analyzer::FindPeak(double &xpeak, double &ypeak, double &xpeak_rebin, double &ypeak_rebin) 
{

  int maxbin = fTrack->GetMaximumBin();
  int x,y,z;
  fTrack->GetBinXYZ(maxbin, x, y, z);

  xpeak=fTrack->GetXaxis()->GetBinCenter(x);
  ypeak=fTrack->GetYaxis()->GetBinCenter(y);

  TH2F* TrackRebin = (TH2F*)fTrack->Clone();
  TrackRebin->Rebin2D(2,2);
  maxbin = TrackRebin->GetMaximumBin();
  TrackRebin->GetBinXYZ(maxbin, x, y, z);
  xpeak_rebin=TrackRebin->GetXaxis()->GetBinCenter(x);
  ypeak_rebin=TrackRebin->GetYaxis()->GetBinCenter(y);

}

//
/*TF1* Analyzer::LeastSquareLine() 
{
  double sum1=0, sum2=0, sum3=0, sum4=0, sum5=0, a=0 , b=0;
  double Z=0;

  for(int i=1;i<fnpixelx;i++)
  {
	for(int j=1;j<fnpixely;j++)
	{  
	  Z=fTrack->GetBinContent(i,j);
	  if(Z!=0)
	  {
		  sum1+= Z*(fTrack->GetXaxis()->GetBinCenter(i))*(fTrack->GetXaxis()->GetBinCenter(i));
		  sum2+= Z*(fTrack->GetYaxis()->GetBinCenter(j));
		  sum3+= Z*(fTrack->GetXaxis()->GetBinCenter(i));
		  sum4+= Z*(fTrack->GetXaxis()->GetBinCenter(i))*(fTrack->GetYaxis()->GetBinCenter(j));
		  sum5+= Z;
	  }
	}
  }

a = (sum1*sum2-sum3*sum4)/(sum5*sum1-(sum3*sum3));
b = (sum5*sum4-sum3*sum2)/(sum5*sum1-(sum3*sum3));

TF1* line = new TF1("leastsquare","[0]*x+[1]",1000,2000);
line->SetParameters(b,a);

return line; 

}*/


//Minifuction launching Atul's script
int Analyzer::Execute_Atul_script(std::string pyvers, std::string inputfile, std::string outfolder, int entries, bool plot, bool text ) const
{
	int i=system("python"+pyvers+" discriminating_vars_BaData.py -I "+inputfile+ " -E "+entries+ " -P "+plot+" -T "+text+"  -O "+ outfolder);
	return i;		//0 if the command was successful, 1 if not
}
