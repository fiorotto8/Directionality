//This in particular compile using  g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum
//Then use as ./nameprog path_to_rootfile

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

using namespace std;



void ScIndicesElem(int nSc,int* ScNelements, vector<int>& B, vector<int>& E){
  B.clear();
  E.clear();

  int parcount=0;

  for(int i=0;i<nSc;i++){
    B.push_back(parcount);
    E.push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}



int main(int argc, char** argv){

  ////////////////////////////////////Get File //////////////////////////////////////
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");

  ///////////////////////////////////Set Branches and Define Variables////////////////////////////////////
  int nmax=250000;
  int nscmax=50;
  int npixel=2304;
  int npixelsmall=250;
  float slimnesslimit=0.6;
  unsigned int nSc;
  int run;
  int event;
  //Pixels

  vector<int> scID;
  scID.reserve(nmax);
  vector<int> scIDall;
  scIDall.reserve(nmax);
  vector<int> ScNpixels;
  ScNpixels.reserve(nscmax);
  vector<float> XPix;
  XPix.reserve(nmax);
  vector<float> YPix;
  YPix.reserve(nmax);
  vector<float> ZPix;
  ZPix.reserve(nmax);
  vector<int> ScNpixelsall;
  ScNpixelsall.reserve(nscmax);
  vector<float> XPixall;
  XPixall.reserve(nmax);
  vector<float> YPixall;
  YPixall.reserve(nmax);
  vector<float> ZPixall;
  ZPixall.reserve(nmax);
  vector<float> width;
  width.reserve(nscmax);
  vector<float> length;
  length.reserve(nscmax);

  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event);
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_ID",scID.data());
  tree->SetBranchAddress("sc_nintpixels",ScNpixels.data());
  tree->SetBranchAddress("sc_ypixelcoord",XPix.data());
  tree->SetBranchAddress("sc_xpixelcoord",YPix.data());
  tree->SetBranchAddress("sc_zpixel",ZPix.data());
  tree->SetBranchAddress("sc_IDall",scIDall.data());
  tree->SetBranchAddress("sc_nallintpixels",ScNpixelsall.data());
  tree->SetBranchAddress("sc_yallpixelcoord",XPixall.data());
  tree->SetBranchAddress("sc_xallpixelcoord",YPixall.data());
  tree->SetBranchAddress("sc_zallpixel",ZPixall.data());
  tree->SetBranchAddress("sc_width",width.data());
  tree->SetBranchAddress("sc_length",length.data());


  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  vector<int> BeginScPix;
  vector<int> EndScPix;
  vector<int> BeginScallPix;
  vector<int> EndScallPix;

  /////////////////////////////////Analysis //////////////////////////////////////////////////////////
  int counter=0;
  int counterall=0;

  for(int k=0;k<tree->GetEntries();k++)
  {
    tree->GetEntry(k);
    cout << "Nev: "<< k << "\n";
    cout << "nSc:  " << nSc << endl;

    ScIndicesElem(nSc,ScNpixels.data(),BeginScPix,EndScPix);
    ScIndicesElem(nSc,ScNpixelsall.data(),BeginScallPix,EndScallPix);

	//Start the cycle on the supercluster of the event
    for(int i=0;i<nSc;i++)
    {
	   if(scID[BeginScPix[i]]>=0)  ////////////preselection that can be changed to one's desire
    {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)	////////////preselection that can be changed to one's desire
      {
			Analyzer Traccia("Name_of_the_histogram",XPix.data(),YPix.data(),ZPix.data(),BeginScPix[i],EndScPix[i]);

			//////////////Do your things//////////////
			counter++;
      }
		}		//end if on overthreshold sc

	   /////////NOW THE TRACK WITH ALL THE PIXELS

	   if(scIDall[BeginScallPix[i]]>=0)			////////////preselection that can be changed to one's desire
    {
		   if(width[i]/length[i]>0.6 && width[i]/length[i]<1)	////////////preselection that can be changed to one's desire
      {

			Analyzer Traccia("Name_of_the_histogram_allpixels",XPixall.data(),YPixall.data(),ZPixall.data(),BeginScallPix[i],EndScallPix[i]);

			////////////////Do your things with the track with also under theshold pixels//////////////

			counterall++;
      }
		}

	}	//end if on supercluster number
  }//end for k (loop on the events)

  return 0;
}
