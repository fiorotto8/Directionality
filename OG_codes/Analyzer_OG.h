#ifndef ANALYZER_H
#define ANALYZER_H

#include <vector>
#include <cmath>
#include <string>
#include "TMatrixD.h"
#include "TDecompSVD.h"

class TH2F;
class TF1;
class TGraph;
class TH1D;

class Analyzer {

 public:
  
  //constructors and destructor
  Analyzer();
  Analyzer(const Analyzer& source);
  Analyzer(const char* nometh2, int npixel ,TH2F* Tracklarge, int npixelorig );		//constructor npixel= number of pixel of the new track, Tracklarge= original 2D histo with the track, npixelorig= number of pixels of the original 2D histogram
  Analyzer(const char* nometh2, TH2F *Tracklarge); 		//constructor for generic track
  Analyzer(const char* nometh2,float* X,float* Y,float* Z, int B,int E);	//Constructor taking the arrays of x,y,z coord of the track. B and E are delimeters for the track hits 
  virtual ~Analyzer();

  //initializer
  void BuildLineMaxRMS();
  void BuildLineDirection();

  //Simple returners
  inline double GetIntegral() const {return fintegral;}
  inline double GetHeight() const {return fheight;}
  inline double GetRadius() const {return fradius;}
  
  inline TH2F* GetHistoTrack() const {return fTrack;}
  inline TH2F* GetHistoTrackTail() const {return fTrackTail;}
  inline TH2F* GetHistoScaledTrack() const {return fScaledTrack;}
  
  inline void GetCentre(double &x, double &y) const {x=fxcentr; y=fycentr;};
  inline double GetXbar() const {return fXbar;}
  inline double GetYbar() const  {return fYbar;}
  inline double GetRMSOnMainAxis() const {return fRMSOnMainAxis;}
  inline double GetSkewOnMainAxis() const {return fSkewOnLine;}
  inline double GetPointSkew(double X, double Y) const {return ( (X-fXbar)*cos(fPhiMainAxis) + (Y-fYbar)*sin(fPhiMainAxis) )/fSkewOnLine;}
  inline TF1* GetLineMaxRMS() const  {return fLineMaxRMS;}
  inline double GetXIP() const {return fXIP;}
  inline double GetYIP() const {return fYIP;}
  inline double GetDir() const {return fPhiDir;}
  inline double PDistCm(double X, double Y) const {return sqrt( ( (X-fXbar)*(X-fXbar) + (Y-fYbar)*(Y-fYbar) ));}

  
  
  //Setters for directionality
  inline void SetWScal(float a) {fwScal=a;}
  inline void SetNPIP(int a) {fNPIP=a;}
  

  //More complex operations
  void Reset();
  void Integral();
  void SavetoFile(const char* nometh2);
  void SavePic(const char* nometh2);
  void SavePicDir(const char* nomepic);
  void SaveRootFile(const char* nomefile);
  
  void Barycenter();
  void Barycenter(TH2F* Tr,double *X, double *Y);
  double AngleLineMaxRMS();
  double RMSOnLine(double Phi);
  double SkewOnMainAxis();
  void RemoveNoise();
  void ApplyThr();
  void ImpactPoint(const char* nometh2);
  void ScaledTrack(const char* nometh2);
  
  void Direction();
  void ImprCorrectAngle();

  void Edges(double &Xl, double &Yl, double &Xr, double &Yr, double slope);
  TH1D* FillProfile(bool longitudinal, float x1=0, float x2=2304);
  TH1D* CutProfile(TH1D* profile, double height=0.0025);
  TH1D* FillProfileX();
  TH1D* FillProfileY();
  void AnglePCA(double &ang);
  void FindNPeaks(TH1D* h, std::vector<std::pair<double,double>> &foundpeaks);
  void FindPeak(double &xpeak, double &ypeak, double &xpeak_rebin, double &ypeak_rebin);
  //TF1* LeastSquareLine();

  int Execute_Atul_script(std::string pyvers, std::string inputfile, std::string outfolder, int entries, bool plot=false, bool text=false ) const;
  
 private:
  
  //histogram parameters
  int fminx,fminy,fmaxx,fmaxy;
  
  
  ///////////////
  double fintegral;
  double fradius;
  double fheight;
  double fxcentr;
  double fycentr;
	
  TH2F *fTrack;
  TH2F *fTrackTail;
  TH2F *fScaledTrack;
  //TH2 Parameters
  int fnpixelx;
  int fnpixely;
  
  //directionaliy parameters
  
  int fNPIP;
  float fwScal;
  
  double fXbar;
  double fYbar;
  TGraph* fBarPlot;
  double fPhiMainAxis;
  TF1* fLineMaxRMS;
  double fRMSOnMainAxis;
  double fSkewOnLine;
  
  double fXIPPrev;
  double fYIPPrev;
  double fXIP;
  double fYIP;
  TGraph* fIPPlot;
  
  
  double fPhiDir;
  TF1* fLineDirection;
  
  
};

#endif

