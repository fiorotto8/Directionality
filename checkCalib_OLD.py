import uproot
import pandas
import numpy as np
import os
import ROOT
from tqdm import tqdm
import awkward as ak
import re
import sys
import pandas as pd
import os
import multiprocessing as mp

ROOT.gROOT.SetBatch(True)

#cuts_offset="(sc_rms>5) & (sc_tgausssigma>2.632) & (sc_width/sc_length>0.8)"
#cuts_offset="(sc_rms>5) & (sc_integral >15000) & (sc_integral <45000) & (sc_length<80) & (sc_tgausssigma>3)& (sc_tgausssigma<10)"


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
def graph(x,y,x_string, y_string,name=None, color=4, markerstyle=22, markersize=2,write=True):
    plot = ROOT.TGraphErrors(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") )
    if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
    else: plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)#blue
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write==True: plot.Write()
    return plot
def extract_numbers(filename):
    # Pattern to match the numbers after 'G' and before 'D', and after 'D'
    pattern = r"G(\d+)D(\d+)"
    # Search for the pattern in the filename
    match = re.search(pattern, filename)
    # If a match is found, extract the numbers
    if match:
        num_after_G = match.group(1)
        num_after_D = match.group(2)
        return int(num_after_G), int(num_after_D)
    else:
        return None, None
def create_fill_TH2(name,x_name,y_name,z_name, x_vals, y_vals, weights=None, x_bins=20, y_bins=20,write=True):
    # Convert x_bins and y_bins to integers explicitly
    x_bins = int(x_bins)
    y_bins = int(y_bins)
    hist = ROOT.TH2F(name, name, x_bins, 0.99*np.min(x_vals), 1.01*np.max(x_vals), y_bins, 0.99*np.min(y_vals), 1.01*np.max(y_vals))
    # Fill the histogram with the data
    if weights is None:
        for x, y in zip(x_vals, y_vals):
            hist.Fill(x, y)
    else:
        for x, y, weight in zip(x_vals, y_vals, weights):
            hist.Fill(x, y, weight)
    # Set axis titles
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle(y_name)
    hist.GetZaxis().SetTitle(z_name)  # Set the z-axis title
    # Draw the histogram
    hist.Draw("COLZ")  # Use the "COLZ" option to draw with a color palette
    if write==True: hist.Write()
    return hist
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
def plot_tgraph2d(x, y, z, title="3D Plot", x_title="X axis", y_title="Y axis", z_title="Z axis",write=True, color=4, markerstyle=22, markersize=2,plot=False,log=False):
    """
    Create and draw a TGraph2D plot.

    Parameters:
    - x, y, z: Arrays of x, y, and z coordinates of the points.
    - title: Title of the plot.
    - x_title, y_title, z_title: Titles for the X, Y, and Z axes.
    - draw_option: Drawing option as a string (e.g., "P" for points, "TRI" for triangles, etc.).

    Returns:
    - The TGraph2D object.
    """
    # Convert input data to numpy arrays if they aren't already
    x_array = np.array(x, dtype="float64")
    y_array = np.array(y, dtype="float64")
    z_array = np.array(z, dtype="float64")
    # Create the TGraph2D object
    graph = ROOT.TGraph2D(len(x_array), x_array, y_array, z_array)
    # Set titles
    graph.SetNameTitle(title,title)
    graph.GetXaxis().SetTitle(x_title)
    graph.GetYaxis().SetTitle(y_title)
    graph.GetZaxis().SetTitle(z_title)
    graph.SetMarkerColor(color)#blue
    graph.SetMarkerStyle(markerstyle)
    graph.SetMarkerSize(markersize)
    # Draw the graph
    graph.Draw("COLZ")
    if write==True: graph.Write()
    if plot==True:
        can1=ROOT.TCanvas("Chi2 values","Chi2 values", 1000, 1000)
        can1.SetFillColor(0);
        can1.SetBorderMode(0);
        can1.SetBorderSize(2);
        can1.SetLeftMargin(0.15);
        can1.SetRightMargin(0.2);
        can1.SetTopMargin(0.1);
        can1.SetBottomMargin(0.1);
        can1.SetFrameBorderMode(0);
        can1.SetFrameBorderMode(0);
        can1.SetFixedAspectRatio();
        if log==True:
            can1.SetLogz()
        graph.Draw("colz")
        can1.Update()
        can1.SaveAs(f"./{title}.png")

    return graph

# Open the file using uproot, much faster for reading
file = "LongRun/long60-40.root"
main = ROOT.TFile("check_60-40.root", "RECREATE")

##################################for 8keV##################################
main.mkdir("8keV")
main.cd("8keV")
#NON MALE cuts_offset="(sc_rms>5) &(sc_integral>15000) & (sc_integral/sc_length>300) & (sc_integral/sc_length<650) & (sc_length<65) & (sc_length>45) & (sc_width<100) & (sc_tgausssigma>3)  & (sc_width/sc_length>0.65)"
#radius in de/dx vs dx elliptic cut
x0,ex0=55,500
xlen,exlen=350,10000
cuts_offset = f"((((sc_length-{x0})*(sc_length-{x0}))/{xlen})+(((sc_integral/sc_length)-{ex0})*((sc_integral/sc_length)-{ex0}))/{exlen})<1"
variables = {"sc_integral": [], "sc_length": [], "sc_width":[],"sc_tgausssigma":[],"sc_rms":[] }
for var in variables:
    values=get_var(var,file)
    variables[var]=values
    #print(len(values))
    hTemp = hist(values, f"60-40_{var}")

create_fill_TH2( f"60-40_dEdx","length","sc_int/length","Entries",variables["sc_length"],variables["sc_integral"]/variables["sc_length"], x_bins=15, y_bins=15)
graph(variables["sc_length"],variables["sc_integral"]/variables["sc_length"],"length","sc_int/length")

hTemp = hist(variables["sc_width"]/variables["sc_length"], "sc_width/sc_length")

##################################for 22keV##################################
main.mkdir("All")
main.cd("All")
#cuts_offset = f"((((sc_length-{x0})*(sc_length-{x0}))/{xlen})+(((sc_integral/sc_length)-{ex0})*((sc_integral/sc_length)-{ex0}))/{exlen})>1"
cuts_offset ="(sc_rms>5)"
variables = {"sc_integral": [], "sc_length": [], "sc_width":[],"sc_tgausssigma":[],"sc_rms":[] }
for var in variables:
    values=get_var(var,file)
    variables[var]=values
    #print(len(values))
    hTemp = hist(values, f"60-40_{var}",channels=1000)

create_fill_TH2( f"60-40_dEdx","length","sc_int/length","Entries",variables["sc_length"],variables["sc_integral"]/variables["sc_length"], x_bins=100, y_bins=100)
graph(variables["sc_length"],variables["sc_integral"]/variables["sc_length"],"length","sc_int/length")

hTemp = hist(variables["sc_width"]/variables["sc_length"], "sc_width/sc_length")





















"""
OPTIMIZATION FOR 8KeV
main.mkdir("optimization")
main.cd("optimization")

x0, ex0 = 55, 500
steps = 5000
gaus = ROOT.TF1("gaus", "gaus(0)", 1000, 70000)

LY, entries,chi2,xlen_values,exlen_values = [],[],[],[],[]
events = uproot.open(file)["Events"]

for i in tqdm(range(steps)):
    xlen=np.random.uniform(200,500)
    exlen=np.random.uniform(8000,12000)
    cuts = f"((((sc_length-{x0})*(sc_length-{x0}))/{xlen})+(((sc_integral/sc_length)-{ex0})*((sc_integral/sc_length)-{ex0}))/{exlen})<1"
    cutVAR = events.arrays("sc_integral", cuts)
    # Utilize Awkward Array's functionality to directly obtain numpy array
    values = ak.to_numpy(ak.flatten(cutVAR["sc_integral"]))
    hTemp = hist(values, f"Sc_int_xlen{xlen}_exlen{exlen}", write=False)
    hTemp.Fit("gaus", "RQ")
    if gaus.GetParameter(1)>27500 and len(values)>800 and (gaus.GetChisquare()/gaus.GetNDF())<1.2 and (gaus.GetChisquare()/gaus.GetNDF())>0.8:
        LY.append(gaus.GetParameter(1))
        chi2.append(gaus.GetChisquare()/gaus.GetNDF())
        entries.append(len(values))
        xlen_values.append(xlen)
        exlen_values.append(exlen)
    if i<5:
        hTemp.Write()

hist(LY,"LY hist")
hist(chi2,"chi2 hist")
hist(entries,"entries hist")

bins=10
plot_tgraph2d( xlen_values, exlen_values, chi2,"chi2 vs xlen and exlen","xlen","exlen","chi2",plot=True)
plot_tgraph2d( xlen_values, exlen_values, LY,"LY vs xlen and exlen","xlen","exlen","LY",plot=True)
plot_tgraph2d( xlen_values, exlen_values, entries,"Entries vs xlen and exlen","xlen","exlen","LY",plot=True)


"""