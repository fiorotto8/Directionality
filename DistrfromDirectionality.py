import uproot
import ROOT
import numpy as np
import pandas as pd
from array import array
import math
import argparse

parser = argparse.ArgumentParser(description="Split a ROOT TTree into multiple subsets.")
parser.add_argument("input_file", help="root file with the directionality analyzed data")
parser.add_argument("calib_file", help="root file with the calibration parameters")
parser.add_argument("output_file", help="Name of the output file")
args = parser.parse_args()

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

#import calibration
#filecalib="LongRun/ArCF4_Mar24/Ar80-20calib.root"
filecalib=args.calib_file
# Open the ROOT file
with uproot.open(filecalib) as file:
    # Access the TTree
    tree = file["calibration_data"]
    # Read the branches as arrays
    offset = tree["offset"].array(library="np")
    offset_error = tree["offset_error"].array(library="np")
    slope = tree["slope"].array(library="np")
    slope_error = tree["slope_error"].array(library="np")

# Now you can use the data just like any other numpy arrays
print("Offset: ", offset[0], "+/-", offset_error[0])
print("Slope: ", slope[0], "+/-", slope_error[0])

#analyzed file
#fileanal="LongRun/ArCF4_Mar24/Analysis_Ar8020.root"
fileanal=args.input_file
# Open the ROOT file
with uproot.open(fileanal) as file:
    # Access the TTree
    tree = file["track_info"]
    # Convert the TTree into a Pandas DataFrame
    df = tree.arrays(library="pd")

#cut on impact point
condition = (df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1900) & (df['Y_ImpactPoint'] > 950) & (df['Y_ImpactPoint'] < 1300) & (df['Ymin'] >550)

df_cut=df[condition]

adc,adc_error=df_cut["Integral"],1E-20*np.ones(len(df_cut["Integral"]))

#I have 
# Calculate energy
energy = offset + adc * slope
# Propagate the error
energy_error = np.sqrt((offset_error)**2 + (adc * slope_error)**2 + (slope * adc_error)**2)

# Adding new columns
df_cut['energy'] = energy
df_cut['err_energy'] = energy_error

#DF CUT CONTAINS everything with cuts
#outfile="Finalmaybe_Ar8020.root"
outfile=args.output_file
# Convert DataFrame to a dictionary for uproot
data_to_write = df_cut.to_dict(orient='list')

# Write the DataFrame to a ROOT file as a TTree
with uproot.recreate(outfile) as f:
    f["mytree"] = data_to_write

main = ROOT.TFile(outfile, "UPDATE")
hist(adc,"sc_integral")
hist(energy,"Energy (keV)")
create_fill_TH2("IP positions","X","Y","Entries",df_cut['X_ImpactPoint'],df_cut['Y_ImpactPoint'],x_bins=20,y_bins=20)
#create_fill_TH2("MIN and MAX positions","X","Y","Entries",pd.concat([df_cut['Xmin'],df_cut['Xmax']]).tolist(),pd.concat([df_cut['Ymin'],df_cut['Ymax']]).tolist(),x_bins=20,y_bins=20)


####DRAW LINES
# Assuming df_cut is your pandas DataFrame
length = 2304

# Create a TCanvas
canvas = ROOT.TCanvas("canvas", "Lines Drawn at Angles", 800, 600)

# Create a TH2F object to define the axes
# Parameters: name, title, n_bins_x, x_min, x_max, n_bins_y, y_min, y_max
frame = ROOT.TH2F("frame", "Lines with Axes;X;Y", 1, 0, 3000, 1, 0, 2304)
frame.Draw()  # Draw the frame to set up the axes

# Dictionaries to store lines and markers
lines = {}
markers = {}

for i, (x, y, angle) in enumerate(zip(df_cut['X_ImpactPoint'], df_cut['Y_ImpactPoint'], df_cut['AnglePCA']), start=1):
    angle_radians = angle  # Assuming angle is already in radians
    slope = math.tan(angle_radians)

    # Calculate Y at X=0 and X=3000
    y_start = slope * (0 - x) + y
    y_end = slope * (3000 - x) + y

    # Create and store the line and marker
    lines[f'line{i}'] = ROOT.TLine(0, y_start, 3000, y_end)
    lines[f'line{i}'].Draw("same")

    markers[f'marker{i}'] = ROOT.TMarker(x, y, 20)
    markers[f'marker{i}'].SetMarkerColor(ROOT.kRed)
    markers[f'marker{i}'].Draw("same")

# Update the canvas to display all lines and axes
canvas.Update()
# Optional: Save the canvas as an image
canvas.SaveAs("lines_at_angles_with_axes.png")
canvas.Write()