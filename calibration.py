import uproot
import numpy as np
import ROOT
import awkward as ak
import argparse

parser = argparse.ArgumentParser(description="Get calibration from Fe and Cd reco datasets")
parser.add_argument("Fe_dataset", help="Input 55Fe reco dataset.")
parser.add_argument("Cd_dataset", help="Input 109Cd reco dataset.")
parser.add_argument("outfile",help="Name of the output root file")
args = parser.parse_args()

cuts_offset="(sc_rms>5) & (sc_tgausssigma>2.632) & (sc_width/sc_length>0.8) & (sc_integral<100000)& (sc_integral>1000)"
def get_var(var,file,cuts=None):
    try:
        events=uproot.open(file+":Events")
    except:
        print("Failed to open (maybe empty)",file)
        return 0
    if cuts is None: cutVAR=events.arrays([var],cuts_offset)
    else: cutVAR=events.arrays([var],f"({var}>{cuts[0]}) & ({var}<{cuts[1]}) & {cuts_offset}")
    VARTemp = [next(iter(d.values())) for d in ak.to_list(cutVAR)]
    VAR = [item for sublist in VARTemp for item in sublist]
    return np.array(VAR,dtype="d")
def hist(data, x_name, channels=100, linecolor=4, linewidth=4,write=True):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, 0.99*np.min(array), 1.01*np.max(array))
    # Use numpy vectorization to fill histogram
    for x in array:
        hist.Fill(x)
    # Set visual attributes and axis titles
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    # Set maximum digits on axes to manage display
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetMaxDigits(3)
    if write:
        hist.Write()
    return hist
def grapherr(x,y,ex,ey,x_string, y_string,name=None, color=4, markerstyle=22, markersize=2,write=True):
    plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") , np.array(   ex   ,dtype="d"),np.array( ey   ,dtype="d"))
    if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
    else: plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)#blue
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write==True: plot.Write()
    return plot

# Open the file using uproot, much faster for reading
#fileFe = "LongRun/ArCF4_Mar24/55FeAr8020.root"
#fileCd = "LongRun/ArCF4_Mar24/109CdAr8020.root"
fileFe = args.Fe_dataset
fileCd = args.Cd_dataset
#main = ROOT.TFile("Ar80-20calib.root", "RECREATE")
main = ROOT.TFile(args.outfile, "RECREATE")


energy,adc,eadc=[5.9,8],[],[]
gaussian=ROOT.TF1("gaussian", "gaus(0)",20000,40000)

sc_int_Fe=get_var("sc_integral",fileFe)
hFe = hist(sc_int_Fe,"sc_integral 55Fe",write=False,channels=500)
#gaussian.SetParameters(70,3500,300,40,9000,400)
hFe.Fit("gaussian","RQ")
hFe.Write()
adc.append(gaussian.GetParameter(1))
eadc.append(gaussian.GetParError(1))

gaussian=ROOT.TF1("gaussian", "gaus(0)",30000,45000)
gausPol=ROOT.TF1("gausPol", "gaus(0)+pol1(3)",30000,45000)

sc_int_Cd=get_var("sc_integral",fileCd)
hCd = hist(sc_int_Cd,"sc_integral 109Cd",write=False,channels=200)
hCd.Fit("gaussian","RQ")
gausPol.SetParameters(gaussian.GetParameter(0),gaussian.GetParameter(1),gaussian.GetParameter(2),0,0)
hCd.Fit("gausPol","RQ")
hCd.Write()
adc.append(gaussian.GetParameter(1))
eadc.append(gaussian.GetParError(1))


# if there are also the double peaks
"""
energy,adc,eadc=[5.9,5.9*2,8,8*2],[],[]
#energy,adc,eadc=[5.9,5.9*2,8*2],[],[]

gaussian=ROOT.TF1("gaussian", "gaus(0)+gaus(3)",3000,15000)

sc_int_Fe=get_var("sc_integral",fileFe)
hFe = hist(sc_int_Fe,"sc_integral 55Fe",write=False,channels=500)
gaussian.SetParameters(70,3500,300,40,9000,400)
hFe.Fit("gaussian","RQ")
hFe.Write()
adc.append(gaussian.GetParameter(1))
adc.append(gaussian.GetParameter(4))
eadc.append(gaussian.GetParError(1))
eadc.append(gaussian.GetParError(4))


sc_int_Cd=get_var("sc_integral",fileCd)
hCdTOT = hist(sc_int_Cd,"sc_integral 109Cd TOT",channels=1000)
hCd = hist(sc_int_Cd,"sc_integral 109Cd",write=False,channels=500)
gaussian.SetParameters(50,3500,200,30,9000,300)
hCd.Fit("gaussian","RQ")
hCd.Write()
adc.append(gaussian.GetParameter(1))
adc.append(gaussian.GetParameter(4))
eadc.append(gaussian.GetParError(1))
eadc.append(gaussian.GetParError(4))
"""





#CALIBRATION CURVE
line=ROOT.TF1("line", "[0]+x*[1]",0,50000)
line.FixParameter(0,0)
line.SetParameters(0,1.2E-3)

calib=grapherr(adc,energy,eadc,[1E-10,1E-10,1E-10,1E-10],"sc_integral","Energy(keV)",write=False)
calib.Fit("line","RQ")
calib.Write()

print("Calibration line: ")
print("Offset: ",line.GetParameter(0),"+/-",line.GetParError(0))
print("Slope: ",line.GetParameter(1),"+/-",line.GetParError(1))


# Create a TTree
calibration_tree = ROOT.TTree("calibration_data", "Calibration Data")

# Create variables to hold the parameter values and their errors
offset = np.array([0.0], dtype="float64")
offset_error = np.array([0.0], dtype="float64")
slope = np.array([0.0], dtype="float64")
slope_error = np.array([0.0], dtype="float64")

# Create branches in the tree for each parameter and its error
calibration_tree.Branch("offset", offset, "offset/D")
calibration_tree.Branch("offset_error", offset_error, "offset_error/D")
calibration_tree.Branch("slope", slope, "slope/D")
calibration_tree.Branch("slope_error", slope_error, "slope_error/D")

# Assign the fit results to the variables
offset[0] = line.GetParameter(0)
offset_error[0] = line.GetParError(0)
slope[0] = line.GetParameter(1)
slope_error[0] = line.GetParError(1)

# Fill the tree with the variables' current contents
calibration_tree.Fill()

# Write the tree to the open ROOT file
calibration_tree.Write()

# After writing, you can proceed to close the file as usual
main.Close()