//This in particular compile using  g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum
//Then use as ./nameprog path_to_rootfile

#include <iostream>
#include <string>
#include <vector>
#include "Analyzer.h"
#include <TTree.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLegendEntry.h>

using namespace std;



/*void ScIndicesElem(int nSc,int* ScNelements, vector<int>& B, vector<int>& E){
  B.clear();
  E.clear();

  int parcount=0;

  for(int i=0;i<nSc;i++){
    B.push_back(parcount);
    E.push_back(parcount+ScNelements[i]);

    parcount+=ScNelements[i];
  }
}*/

//for reduced pixels
void ScIndicesElem(int nSc, UInt_t npix, float* sc_redpixID, int &nSc_red, vector<int>& B, vector<int>& E){
  B.clear();
  E.clear();

  vector<float> sc_redpix_start = {0};

  int parcount=0;

  for(int i=0; i<npix; i++){
    cout<<i<<endl;
    if(sc_redpixID[i]>0) sc_redpix_start.push_back(sc_redpixID[i]);
  }

  nSc_red = sc_redpix_start.size();

  sc_redpix_start.push_back(npix);

  for(int i=0;i<sc_redpix_start.size()-1;i++){
    B.push_back(sc_redpix_start[i]);
    E.push_back(sc_redpix_start[i+1]);
    //std::cout<<B[i]<<" "<<E[i]<<endl;
  }

  sc_redpix_start.clear();

}

double distance(double x1, double y1, double x2, double y2){

  return sqrt(pow((x2-x1),2)+pow((y2-y1),2));

}

int closest_divisor(int A, float B){ //Divisor of A closest to B

  vector<float>  div_dist;
  vector<int> div_candidate;

  div_dist.clear();
  div_candidate.clear();

  for(int i=1; i<=A; i++){
    if(A%i==0){ //if it's a divisor
      div_dist.push_back(abs((float)i-B));
      div_candidate.push_back(i);
    }
  }

  auto result = min_element(div_dist.begin(), div_dist.end());  

  return div_candidate[result-div_dist.begin()];

}

int main(int argc, char** argv){

  ////////////////////////////////////Get File //////////////////////////////////////
  TFile* f = TFile::Open(Form("%s",argv[1]));
  TTree* tree = (TTree*)f->Get("Events");
  if (argv[3]=="MC") cout<<"saving MC variables"<<endl;
  ///////////////////////////////////Set Branches and Define Variables////////////////////////////////////
  int nmax=500000;
  int nscmax=100;
  int npixel=2304;
  int npixelsmall=250;
  //float slimnesslimit=0.6;

  unsigned int nSc=0;
  int nSc_red=0;
  UInt_t Nredpix=0;
  Int_t sc_npix=0;
  int run;
  int event;

  //Pixels

  vector<float> sc_redpixID;
  sc_redpixID.reserve(nmax);
  vector<UInt_t> ScNpixels;
  ScNpixels.reserve(nscmax);
  vector<int> XPix;
  XPix.reserve(nmax);
  vector<int> YPix;
  YPix.reserve(nmax);
  vector<float> ZPix;
  ZPix.reserve(nmax);
  vector<float> xmean;
  xmean.reserve(nscmax);
  vector<float> ymean;
  ymean.reserve(nscmax);
  vector<float> ymin;
  ymin.reserve(nscmax);
  vector<float> ymax;
  ymax.reserve(nscmax);
  vector<float> xmin;
  xmin.reserve(nscmax);
  vector<float> xmax;
  xmax.reserve(nscmax);
  vector<float> scsize; 
  scsize.reserve(nscmax);
  vector<float> scnhits;
  scnhits.reserve(nscmax);
  vector<float> v_sc_rms;
  v_sc_rms.reserve(nscmax);
  vector<float> v_sc_tgausssigma;
  v_sc_tgausssigma.reserve(nscmax);
  vector<float> v_sc_theta;
  v_sc_theta.reserve(nscmax);

  vector<float> width;
  width.reserve(nscmax);
  vector<float> length;
  length.reserve(nscmax);
  vector<float> integral;
  integral.reserve(nscmax);

  //MC variables
  Float_t MC_Xvertex=0;
  Float_t MC_Yvertex=0;
  Float_t MC_Zvertex=0;
  Int_t MC_energy=0;
  Float_t MC_phi=0;
  Float_t MC_theta=0;
  Int_t MC_particleID=0;
  Float_t MC_Xvertex_out=0;
  Float_t MC_Yvertex_out=0;
  Float_t MC_Zvertex_out=0;
  Float_t MC_energy_out=0;
  Float_t MC_phi_out=0;
  Float_t MC_theta_out=0;
  Int_t MC_particleID_out=0;
    
  tree->SetBranchAddress("run",&run);
  tree->SetBranchAddress("event",&event); 
  /////////////Reco variables//////////////     
  tree->SetBranchAddress("nSc",&nSc);
  tree->SetBranchAddress("sc_redpixIdx",sc_redpixID.data());
  //tree->SetBranchAddress("sc_size",&sc_npix);
  tree->SetBranchAddress("nRedpix",&Nredpix);
  tree->SetBranchAddress("redpix_ix",YPix.data());
  tree->SetBranchAddress("redpix_iy",XPix.data());
  tree->SetBranchAddress("redpix_iz",ZPix.data());
  tree->SetBranchAddress("sc_width",width.data());
  tree->SetBranchAddress("sc_length",length.data());
  tree->SetBranchAddress("sc_integral",integral.data());
  tree->SetBranchAddress("sc_xmean",xmean.data());
  tree->SetBranchAddress("sc_ymean",ymean.data());
  tree->SetBranchAddress("sc_xmin",xmin.data());
  tree->SetBranchAddress("sc_ymin",ymin.data());
  tree->SetBranchAddress("sc_xmax",xmax.data());
  tree->SetBranchAddress("sc_ymax",ymax.data());
  tree->SetBranchAddress("sc_size",scsize.data());
  tree->SetBranchAddress("sc_nhits",scnhits.data());
  tree->SetBranchAddress("sc_theta",v_sc_theta.data());
  tree->SetBranchAddress("sc_rms",v_sc_rms.data());
  tree->SetBranchAddress("sc_tgausssigma",v_sc_tgausssigma.data());
  //if (argv[3]=="MC"){
        tree->SetBranchAddress("MC_particle_type",&MC_particleID);
        tree->SetBranchAddress("MC_cutexposure_energy",&MC_energy);
        tree->SetBranchAddress("MC_phi",&MC_phi);
        tree->SetBranchAddress("MC_theta",&MC_theta);
        tree->SetBranchAddress("MC_x_vertex",&MC_Xvertex);
        tree->SetBranchAddress("MC_y_vertex",&MC_Yvertex);
        tree->SetBranchAddress("MC_z_vertex",&MC_Zvertex);

   // }

  /////////////////////////////////Analysis Variables ////////////////////////////////////////////////
  vector<int> BeginScPix;
  vector<int> EndScPix;
  //vector<int> BeginScallPix;
  //vector<int> EndScallPix;

  double phi_PCA=0, phi_RMS=0, phi=0;
  double phi_HT_skew=0, phi_HT_maxpx=0, phi_HT_maxpxRebin=0, phi_HT_peakprof=0, phi_HT_maxprof=0, phi_HT_integral=0;
  double scint = 0;
  double skewness_L=0, skewness_T=0, kurtosis_L=0, kurtosis_T=0;
  double xbar=0, ybar=0;
  double xl=0, yl=0, xr=0, yr=0;
  double reco_theta=0, x_mean=0, y_mean=0, x_min=0, x_max=0, y_min=0, y_max=0;
    
  //impact point and directionality
  Int_t NPIP;
  Float_t wFac;
  double xIP=0, yIP=0;

  vector<std::pair<double,double>> peakslong;
  vector<std::pair<double,double>> peakstrans;
  int npeaks=0;
  double peak_density=0, tracklength=0, track_width=0, recolength=0, recowidth=0, reco_sc_rms=0, reco_sc_tgausssigma=0;
  vector<double> peak_distance, peak_width;

//temporanee; per check vari
int htcheck_skew=0, htcheck_peak=0, htcheck_macropeak=0, htcheck_integral=0, htcheck_profpeak=0, htcheck_profmax=0;
double xpr, ypr, xp, yp;

  /////////////////////////////////Histograms/////////////////////////////////////////////////////////
  TH1F* h_Integral = new TH1F("h_Integral","Integral",2000,0,200000);
  TH1F* h_AngDistribution = new TH1F("h_AngDistribution","AngDistribution",180,0,360);
  TH2F* polar_hist = new TH2F("polar_hist","polar_hist", 90, 0, 360, 45, 1000, 10000);

  TH1F* h_npeaks = new TH1F("h_npeaks","npeaks",1000,0,1000);
  TH1F* h_sigma_peaks = new TH1F("h_sigma_peaks","sigma_peaks",80,0,40);
  TH1F* h_peak_density = new TH1F("h_peak_density","peak_density",500,0,0.5);
  TH1F* h_peak_distance = new TH1F("h_peak_distance","peak_distance",150,0,150);

  
  /////////////////////////////////Analysis //////////////////////////////////////////////////////////
  int counter=0;
  int ired=0;
  //int counterall=0;

  ///////output file////////
  string filename(argv[1]);
  filename = filename.substr(filename.find_last_of("/\\")+1);
cout<<filename<<endl;
  TFile* fout = new TFile(Form("%s/Analysis_%s",argv[2],filename.c_str()),"recreate");
cout<<Form("%s/Analysis_%s",argv[2],filename.c_str())<<endl;
  fout->cd();
  fout->mkdir("Tracks");

  TTree* out_tree = new TTree("track_info","track_info");
  out_tree->Branch("run",&run);
  out_tree->Branch("event",&event);
  out_tree->Branch("nSc",&nSc);
  out_tree->Branch("nSc_red",&nSc_red);
  out_tree->Branch("Integral",&scint);
  out_tree->Branch("ScSize",&sc_npix);
  out_tree->Branch("AnglePCA",&phi_PCA);
  out_tree->Branch("AngleRMS",&phi_RMS);
  out_tree->Branch("RecoTheta",&reco_theta);
  out_tree->Branch("RecoScRMS",&reco_sc_rms);
  out_tree->Branch("RecoScTGaussSigma",&reco_sc_tgausssigma);
  //out_tree->Branch("X_ImpactPoint",&xIP);
  //out_tree->Branch("Y_ImpactPoint",&yIP);
  out_tree->Branch("Xmean",&x_mean);
  out_tree->Branch("Ymean",&y_mean);
  out_tree->Branch("Xmin",&x_min);
  out_tree->Branch("Ymin",&y_min);
  out_tree->Branch("Xmax",&x_max);
  out_tree->Branch("Ymax",&y_max);
  out_tree->Branch("XBar",&xbar);
  out_tree->Branch("YBar",&ybar);
  out_tree->Branch("Npeaks",&npeaks);
  out_tree->Branch("PeakDensity",&peak_density);
  out_tree->Branch("TrackLength",&tracklength);
  out_tree->Branch("TrackWidth",&track_width);
  out_tree->Branch("RecoLength",&recolength);
  out_tree->Branch("RecoWidth",&recowidth);
  out_tree->Branch("PeakDistance",&peak_distance);
  out_tree->Branch("Skewness",&skewness_L);
  out_tree->Branch("Kurtosis",&kurtosis_T);
  //if (argv[3]=="MC"){
        out_tree->Branch("MC_particle_type",&MC_particleID_out);
        out_tree->Branch("MC_energy",&MC_energy_out);
        out_tree->Branch("MC_phi",&MC_phi_out);
        out_tree->Branch("MC_theta",&MC_theta_out);
        out_tree->Branch("MC_Xvertex",&MC_Xvertex_out);
        out_tree->Branch("MC_Yvertex",&MC_Yvertex_out);
        out_tree->Branch("MC_Zvertex",&MC_Zvertex_out);
       
   // }
  out_tree->Branch("phi_HT_skew",&phi_HT_skew);
  out_tree->Branch("phi_HT_maxpx",&phi_HT_maxpx);
  out_tree->Branch("phi_HT_maxpxRebin",&phi_HT_maxpxRebin);    
  out_tree->Branch("phi_HT_peakprof",&phi_HT_peakprof);
  out_tree->Branch("phi_HT_maxprof",&phi_HT_maxprof);
  out_tree->Branch("phi_HT_integral",&phi_HT_integral);
    
  cout<<"this run has "<<tree->GetEntries()<<" entries"<<endl;
  for(int k=0;k<tree->GetEntries();k++)
  {
    cout<<"getting entries..."<<endl;
    tree->GetEntry(k);
    cout << "Nev: "<< k << "\nnSc:  " << nSc << " event "<< event <<endl;
    //for reduced pixels:
    ScIndicesElem(nSc,Nredpix,sc_redpixID.data(),nSc_red,BeginScPix,EndScPix);
      

    cout<<"reduced sc "<<nSc_red<<endl;
      //cout<<"MC energy "<<MC_energy<<endl;
	//Start the cycle on the supercluster of the event
      ired=0;
      for(int i=0;i<nSc;i++)
      {   

	//cuts:
	//cosmic selection: cut1=length[i]>500 || width[i]/length[i]<0.3; cut2=integral[i]/ScNpixels[i]<5
	//59keV photons from Am: cut1 = integral[i]/ScNpixels[i] > (12.-length[i]/50.) && integral[i]/ScNpixels[i] < (16.-length[i]/50.); cut2 = length[i] > 90. && length[i] < 250.
	//AmBe selection: if not cosmic && not 59keV photons && density>=10 (50% efficiency on signal, 1% on background)

//	bool cut1 = integral[i]/ScNpixels[i] > (12.-length[i]/50.) && integral[i]/ScNpixels[i] < (16.-length[i]/50.) && length[i] > 90. && length[i] < 250.;
//	bool cut2 = integral[i]/ScNpixels[i] >= 10;
bool cut1=0;
//bool cut2=xmean[i]>800 && xmean[i]<1800 && ymean[i]>500 && ymean[i]<1800;
bool cut2=1;
          
		  scint = integral[i];
		  recolength=length[i];
		  recowidth=width[i];
		  sc_npix = scnhits[i];
          reco_theta=v_sc_theta[i];
          x_mean=xmean[i];
          y_mean=ymean[i];
          x_min=xmin[i];
          x_max=xmax[i];
          y_min=ymin[i];
          y_max=ymax[i];
          reco_sc_rms=v_sc_rms[i];
          reco_sc_tgausssigma=v_sc_tgausssigma[i];
          //cout<<i<<" Sc out of "<<nSc<<endl;

       cut2=reco_sc_rms>6 && scint>1000 && reco_sc_tgausssigma*0.152>0.5 && x_min>400 && x_max<1900 && y_min>400 && y_max<1900;
          
	   if(sc_redpixID[i]!=-1)  ////////////preselection that can be changed to one's desire
	   {
		   if(!cut1) {	////////////preselection that can be changed to one's desire

		   if(cut2) {	////////////preselection that can be changed to one's desire
	
			////Construct the track
            cout<<"constructing the track"<<endl;

			//cout<<BeginScPix[ired]<<" - "<<EndScPix[ired]<<endl;
			if (EndScPix[ired]-BeginScPix[ired]<=0) continue;

            if (TMath::MinElement(EndScPix[ired]-BeginScPix[ired],&XPix[BeginScPix[ired]])>400 && TMath::MinElement(EndScPix[ired]-BeginScPix[ired],&YPix[BeginScPix[ired]])>400 && TMath::MaxElement(EndScPix[ired]-BeginScPix[ired],&XPix[BeginScPix[ired]])<1900 && TMath::MaxElement(EndScPix[ired]-BeginScPix[ired],&YPix[BeginScPix[ired]])<1900) {
                
			Analyzer Traccia(Form("Track%i_event%i_%i",counter,event,run),XPix.data(),YPix.data(),ZPix.data(),BeginScPix[ired],EndScPix[ired]);
            //Analyzer Traccia(Form("Track%i_event%i_%i",counter,event,run),YPix.data(),XPix.data(),ZPix.data(),BeginScPix[ired],EndScPix[ired]);

            cout<<Traccia.GetXmin()<<" "<<Traccia.GetXmax()<<" "<<Traccia.GetYmin()<<" "<<Traccia.GetYmax()<<endl;
            //if (Traccia.GetXmin()>400 && Traccia.GetYmin()>400 && Traccia.GetXmax()<1900 && Traccia.GetYmax()<1900){
			Traccia.RemoveNoise();
            

			///////Get some properties//////
			/*scint = integral[i];
			recolength=length[i];
			recowidth=width[i];
			sc_npix = scnhits[i];*/
			phi_RMS = Traccia.AngleLineMaxRMS();
            phi_PCA = Traccia.AnglePCA();
			phi = Traccia.AnglePCA(); //to apply the head tail
			xbar = Traccia.GetXbar(); ybar = Traccia.GetYbar();
			Traccia.Edges(xl,yl,xr,yr,tan(phi));

               
            MC_Xvertex_out=MC_Xvertex;
            MC_Yvertex_out=MC_Yvertex;
            MC_Zvertex_out=MC_Zvertex;
            MC_energy_out=MC_energy;
            MC_phi_out=MC_phi;
            MC_theta_out=MC_theta;
            MC_particleID_out=MC_particleID;
            
               
			///////Profiles///////
			TH1D* LongProfile = Traccia.FillProfile(1);
			LongProfile->SetTitle(Form("TrackLongProf%i_event%i_%i",counter,event,run)); LongProfile->SetName(LongProfile->GetTitle());
			TH1D* TransProfile = Traccia.FillProfile(0);
			TransProfile->SetTitle(Form("TrackTransProf%i_event%i_%i",counter,event,run)); TransProfile->SetName(TransProfile->GetTitle());
			
			TH1D* XProfile = Traccia.FillProfileX();
			XProfile->SetTitle(Form("TrackXProf%i_event%i_%i",counter,event,run)); XProfile->SetName(XProfile->GetTitle());
			TH1D* YProfile = Traccia.FillProfileY();
			YProfile->SetTitle(Form("TrackYProf%i_event%i_%i",counter,event,run)); YProfile->SetName(YProfile->GetTitle());


			/////Asymmetry//////
			//actual edges of track
			int binmin_L = LongProfile->FindFirstBinAbove(0.0025*LongProfile->GetMaximum(),1); //first bin above 0.25% of max intensity
			int binmax_L = LongProfile->FindLastBinAbove(0.0025*LongProfile->GetMaximum(),1); //last bin above 0.25% of max intensity
			int binmin_T = TransProfile->FindFirstBinAbove(0.0025*TransProfile->GetMaximum(),1); //first bin above 0.25% of max intensity
			int binmax_T = TransProfile->FindLastBinAbove(0.0025*TransProfile->GetMaximum(),1); //last bin above 0.25% of max intensity

			TH1D* LongProfile_cut = Traccia.CutProfile(LongProfile,1);
			TH1D* TransProfile_cut = Traccia.CutProfile(TransProfile,0);
			LongProfile_cut->SetTitle(Form("TrackLongProf_cut%i_event%i_%i",counter,event,run)); 	LongProfile_cut->SetName(LongProfile_cut->GetTitle());
			TransProfile_cut->SetTitle(Form("TrackTransProf_cut%i_event%i_%i",counter,event,run)); 	TransProfile_cut->SetName(TransProfile_cut->GetTitle());

			double xmin_online = xl+LongProfile->GetXaxis()->GetBinCenter(binmin_L)/sqrt(1+tan(phi)*tan(phi));
			double ymin_online = yl+tan(phi)*LongProfile->GetXaxis()->GetBinCenter(binmin_L)/sqrt(1+tan(phi)*tan(phi));
			double xmax_online = xl+LongProfile->GetXaxis()->GetBinCenter(binmax_L)/sqrt(1+tan(phi)*tan(phi));
			double ymax_online = yl+tan(phi)*LongProfile->GetXaxis()->GetBinCenter(binmax_L)/sqrt(1+tan(phi)*tan(phi));

			tracklength = distance(xmin_online,ymin_online,xmax_online,ymax_online);

			skewness_L=LongProfile_cut->GetSkewness(1);
			skewness_T=TransProfile_cut->GetSkewness(1);
			kurtosis_L=LongProfile_cut->GetKurtosis(1);
			kurtosis_T=TransProfile_cut->GetKurtosis(1);


			Traccia.FindPeak(xp,yp,xpr,ypr); //max intensity and max intensity from rebin
			double xpr_online = (1/(1+pow(tan(phi),2)))*(xpr+xbar*pow(tan(phi),2)+tan(phi)*(ypr-ybar));
			double ypr_online = ybar+(tan(phi)/(1+pow(tan(phi),2)))*(xpr-xbar+tan(phi)*(ypr-ybar));
			double xp_online = (1/(1+pow(tan(phi),2)))*(xpr+xbar*pow(tan(phi),2)+tan(phi)*(ypr-ybar));
			double yp_online = ybar+(tan(phi)/(1+pow(tan(phi),2)))*(xpr-xbar+tan(phi)*(ypr-ybar));

			peakslong.clear();
			peak_distance.clear();
			peak_width.clear();

			Traccia.FindProfilePeaks(LongProfile_cut, peakslong);  //peak from profile
			Traccia.FindProfilePeaks(TransProfile_cut, peakstrans);
			
			TF1* fittrans = new TF1("fittrans","gaus",0,TransProfile_cut->GetNbinsX()-1);
			TransProfile_cut->Fit("fittrans","Q","",peakstrans[0].first-5*peakstrans[0].second,peakstrans[0].first+5*peakstrans[0].second);
			track_width = 4*fittrans->GetParameter(2);

			npeaks = peakslong.size();
               
            if (npeaks==0) {
			delete LongProfile;
			delete TransProfile;
			delete LongProfile_cut;
			delete TransProfile_cut;
			delete XProfile;
			delete YProfile;
            continue;}
               
			h_npeaks->Fill(npeaks);
			peak_density = npeaks/tracklength;
			h_peak_density->Fill(peak_density);
			sort(peakslong.begin(), peakslong.end());
			TF1* fitfunc = new TF1("fitfunc","gaus",0,LongProfile_cut->GetNbinsX()-1);

			for(int p=0; p<npeaks; p++){
				fitfunc->SetParameters(peakslong[p].second*LongProfile_cut->GetBinContent(LongProfile_cut->FindBin(peakslong[p].first)),peakslong[p].first,peakslong[p].second);
				LongProfile_cut->Fit("fitfunc","Q","",peakslong[p].first-5*peakslong[p].second,peakslong[p].first+5*peakslong[p].second);
				peak_width.push_back(fitfunc->GetParameter(2));
				h_sigma_peaks->Fill(fitfunc->GetParameter(2));
				if(p+1<peakslong.size()) {
				    peak_distance.push_back(peakslong[p+1].first-peakslong[p].first);
				    h_peak_distance->Fill(peakslong[p+1].first-peakslong[p].first);
				}
				//if(LongProfile_cut->GetBinContent(picco)<LongProfile_cut->GetBinContent(LongProfile_cut->FindBin(fitfunc->GetParameter(1)))) picco=fitfunc->GetParameter(1);
			}
			double xpos_peak_online = xmin_online+peakslong[0].first/sqrt(1+tan(phi)*tan(phi));
			double ypos_peak_online = ymin_online+tan(phi)*peakslong[0].first/sqrt(1+tan(phi)*tan(phi));
			double prof_max = LongProfile_cut->GetXaxis()->GetBinCenter(LongProfile_cut->GetMaximumBin()); //max from profile
			double xprof_max_online = xmin_online+prof_max/sqrt(1+tan(phi)*tan(phi));
			double yprof_max_online = ymin_online+tan(phi)*prof_max/sqrt(1+tan(phi)*tan(phi));

			///////Head-tail//////// assuming a low energy NR (the majority of energy loss is at the beginning of track)
			//from skewness : 1 
			if (skewness_L<0) { //track from R to L, I expect negative skewness
			  phi_HT_skew = phi + TMath::Pi();
			} 
			else if (skewness_L>=0) { //track from L to R, I expect positive skewness
			  if(phi < 0) phi_HT_skew = phi + 2*TMath::Pi();
              else phi_HT_skew = phi;
			}

			//from max intensity pixel : 2
			if(distance(xp_online,yp_online,xmin_online,ymin_online) < 0.5*tracklength) { //from L to R
			  if(phi < 0) phi_HT_maxpx = phi + 2*TMath::Pi();
              else phi_HT_maxpx = phi;
			}
			else if(distance(xp_online,yp_online,xmin_online,ymin_online) >= 0.5*tracklength) {//from R to L
			  phi_HT_maxpx = phi + TMath::Pi();
			} 
			

            //from max intensity macro-pixel : 3
			if(distance(xpr_online,ypr_online,xmin_online,ymin_online) < 0.5*tracklength) { //from L to R
			  if(phi < 0) phi_HT_maxpxRebin = phi + 2*TMath::Pi();
              else phi_HT_maxpxRebin = phi;
			}
			else if(distance(xpr_online,ypr_online,xmin_online,ymin_online) >= 0.5*tracklength) { //from R to L
			  phi_HT_maxpxRebin = phi + TMath::Pi();
			}

			//from peak found in profile : 4
			if(distance(xpos_peak_online,ypos_peak_online,xmin_online,ymin_online) < 0.5*tracklength) { //from L to R
			  if(phi < 0) phi_HT_peakprof = phi + 2*TMath::Pi();
              else phi_HT_peakprof = phi;
			}
			else if(distance(xpos_peak_online,ypos_peak_online,xmin_online,ymin_online) >= 0.5*tracklength) { //from R to L
			  phi_HT_peakprof = phi + TMath::Pi();
			}
			//from max found in profile : 5
			if(distance(xprof_max_online,yprof_max_online,xmin_online,ymin_online) < 0.5*tracklength) { //from L to R
			  if(phi < 0) phi_HT_maxprof = phi + 2*TMath::Pi();
              else phi_HT_maxprof = phi;
			}
			else if(distance(xprof_max_online,yprof_max_online,xmin_online,ymin_online) >= 0.5*tracklength) { //from R to L
			  phi_HT_maxprof = phi + TMath::Pi();
			}
			//from integral : 6
			int center_bin = LongProfile_cut->FindBin(0.5*tracklength);
			if((LongProfile_cut->Integral(1,center_bin) > LongProfile_cut->Integral()-LongProfile_cut->Integral(1,center_bin))) { //from L to R
			  if(phi < 0) phi_HT_integral = phi + 2*TMath::Pi();
              else phi_HT_integral = phi;
			}
			else if ((LongProfile_cut->Integral(1,center_bin) <= LongProfile_cut->Integral()-LongProfile_cut->Integral(1,center_bin))) { //from R to L
			  phi_HT_integral = phi + TMath::Pi();
			}


			///////Fill histograms////////
			
			h_Integral->Fill(scint); 

			h_AngDistribution->Fill(180*phi/TMath::Pi());
			polar_hist->Fill(180*phi/TMath::Pi(), integral[i]);

//			TF1* TransProfGauss = new TF1("TransProfGauss","gaus",0,TransProfile_cut->GetNbinsX()-1);
//			TransProfile_cut->Fit("TransProfGauss","Q R");
//			h_TransProf_width->Fill(TransProfGauss->GetParameter("Sigma"),MC_drift[i]);

			///////Draw and save the tracks//////////
			//fout->cd();
			fout->cd("Tracks");

			TCanvas *c1=new TCanvas(Form("TrackDirection%i_event%i_%i",counter,event,run),Form("TrackDirection%i_event%i_%i",counter,event,run),1000,1000);

			TH2F *TracciaHisto=Traccia.GetHistoTrack();
			TracciaHisto->Draw("colz");
			TF1* lineRMS = (TF1*)Traccia.GetLineMaxRMS();
			lineRMS->Draw("same");

            TF1* linePCA = new TF1("LinePCA","[0]*(x-[1])+[2]",0,2304);

            linePCA->SetParameter(1,xbar);
            linePCA->SetParameter(2,ybar);
            linePCA->SetParameter(0,TMath::Tan(phi));
            linePCA->SetLineColor(kBlue);
            linePCA->Draw("same");
  

			//Left extreme point along main line
			//TMarker ml(xl,yl,20);
			TMarker ml(xmin_online,ymin_online,20);
			ml.SetMarkerColor(kBlack);
			ml.SetMarkerSize(1);
			//Right extreme point along main line
			//TMarker mr(xr,yr,20);
			TMarker mr(xmax_online,ymax_online,20);
			mr.SetMarkerColor(kBlack);
			mr.SetMarkerSize(1);
			//barycenter
			TMarker m(xbar,ybar,20);
			m.SetMarkerColor(kBlack);
			m.SetMarkerSize(1);
			//max intensity pixel
			TMarker mp(xp_online,yp_online,20);
			mp.SetMarkerColor(kGreen);
			mp.SetMarkerSize(1);
			//max intensity macropixel
			TMarker mpr(xpr_online,ypr_online,20);
			mpr.SetMarkerColor(kRed);
			mpr.SetMarkerSize(1);
			//peak from profile
			TMarker mprof_peak(xpos_peak_online,ypos_peak_online,20);
			mprof_peak.SetMarkerColor(kRed+2);
			mprof_peak.SetMarkerSize(1);
			//max from profile
			TMarker mprof_max(xprof_max_online,yprof_max_online,20);
			mprof_max.SetMarkerColor(kGreen+2);
			mprof_max.SetMarkerSize(1);

			ml.Draw("same");
			mr.Draw("same");
			m.Draw("same");
			mp.Draw("same");
			mpr.Draw("same");
			mprof_max.Draw("same");
			mprof_peak.Draw("same");

			//BuildLegend
			string title;
			auto legend = new TLegend(0.1,0.7,0.48,0.9);
			TLegendEntry *e2 = legend->AddEntry(&m,"Barycenter","p");
			TLegendEntry *e3 = legend->AddEntry(&mp,"Intensity peak","p");
			TLegendEntry *e1 = legend->AddEntry(&mpr,"Intensity peak - 2x2","p");
			TLegendEntry *e0 = legend->AddEntry(&mprof_peak,"Intensity peak - from profile","p");
			TLegendEntry *e6 = legend->AddEntry(&mprof_max,"Intensity max - from profile","p");
			title = "RMS phi = "+to_string((int)(phi_RMS*180/TMath::Pi()))+"\xB0";
			TLegendEntry *e4 = legend->AddEntry(lineRMS,&title[0],"l");
			title = "PCA phi = "+to_string((int)(phi*180/TMath::Pi()))+"\xB0";
			TLegendEntry *e5 = legend->AddEntry(linePCA,&title[0],"l");

			legend->Draw();

            if (counter%500==0){
                c1->Write();
                LongProfile->Write();
                LongProfile_cut->Write();
                TransProfile->Write();
                TransProfile_cut->Write();
                XProfile->Write();
                YProfile->Write();
            }



			delete c1;


			delete LongProfile;
			delete TransProfile;
			delete LongProfile_cut;
			delete TransProfile_cut;
			delete XProfile;
			delete YProfile;
               
            //cout<<"finally computing IP"<<endl;
    /*        Traccia.SetWScal(wFac);
            Traccia.SetNPIP(NPIP);
            Traccia.ApplyThr();
            Traccia.ImpactPoint(Form("Track%i_IP_event%i_%i",counter,event,run));
            xIP=Traccia.GetXIP();
            yIP=Traccia.GetYIP();
               cout<<"got IP"<<endl;*/
            cout<<"fine traccia"<<endl;

			out_tree->Fill();
			counter++;
			ired++;
			
            }
			//cout<<"tree filled"<<endl;
		   } //cut 1	
		   } //cut 2
		}		//end if on overthreshold sc
        //out_tree->Fill();
	
	}	//end if on supercluster number
	
//if(counter>=300){break;}

  }//end for k (loop on the events)

  fout->cd();
  out_tree->Write();

  /////////Fit histograms/////////////

  /////////Save histograms and graphs to output///////
  h_Integral->GetXaxis()->SetTitle("Integral");
  h_AngDistribution->GetXaxis()->SetTitle("Angle phi [degrees]");

  h_npeaks->GetXaxis()->SetTitle("number of peaks");
  h_sigma_peaks->SetTitle("sigma of peaks");
  h_peak_density->SetTitle("peak density (N of peaks/length)");
  h_peak_distance->SetTitle("peak distance");

  //h_Integral->Write();
  //h_AngDistribution->Write();
  //h_Calibration->Write();
  //h_Skewness_Long->Write();
  //h_Skewness_Trans->Write();
  //h_Kurtosis_Long->Write();
  //h_Kurtosis_Trans->Write();
  //h_TransProf_width->Write();
  //h_HTefficiency->Write();
  //h_npeaks->Write();
  //h_sigma_peaks->Write();
  //h_peak_density->Write();
  //h_peak_distance->Write();

  /*delete h_Integral;
  delete h_AngDistribution;
  delete h_npeaks;
  delete h_sigma_peaks;
  delete h_peak_density;
  delete h_peak_distance;
*/

  //polar_hist->Write();

  fout->Close();


  return 0;
}
