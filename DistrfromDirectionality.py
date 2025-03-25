import uproot
import ROOT
import numpy as np
import pandas as pd
from array import array
import math
import argparse
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

parser = argparse.ArgumentParser(description="Split a ROOT TTree into multiple subsets.")
parser.add_argument("input_file", help="root file with the directionality analyzed data")
parser.add_argument("calib_file", help="root file with the calibration parameters")
parser.add_argument("output_file", help="Name of the output file")
parser.add_argument('--intr_file', type=str, help='Path to the intrinsic distribution file, if passed it will align the energy distribution',default=None)
parser.add_argument("--sim",help="Flag for simulated datasets",action="store_true")
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
def nparr(list):
    return np.array(list,dtype="d")

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
if args.sim:
    #! In simulation no cut is needed since it is done in preprocessing in the ConvertForDigi
    condition = (df['X_ImpactPoint'] > 0) & (df["Integral"]>1000)
    #condition = (df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1800) & (df['Y_ImpactPoint'] > 850) & (df['Y_ImpactPoint'] < 1400)
    #condition = (df['X_ImpactPoint'] > 0) 
    df_cut_temp=df[condition]

    #condition = (df['X_ImpactPoint'] > 0)
    #condition = (df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1800) & (df['Y_ImpactPoint'] > 850) & (df['Y_ImpactPoint'] < 1400) & (df['Ymin'] >500) & (df['Ymax'] <2304-500) & (df["Xmin"]>500)
    #condition = (df['X_ImpactPoint'] > 1900) & (df['X_ImpactPoint'] < 2000) & (df['Y_ImpactPoint'] > 600) & (df['Y_ImpactPoint'] < 1600) & (df['Ymin'] >500) & (df['Ymax'] <2304-500) & (df["Xmin"]>500)
    df_cut=df[condition]
else:
    condition = (df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1800) & (df['Y_ImpactPoint'] > 850) & (df['Y_ImpactPoint'] < 1400)
    df_cut_temp=df[condition]

    condition = (df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1800) & (df['Y_ImpactPoint'] > 850) & (df['Y_ImpactPoint'] < 1400) & (df['Ymin'] >500) & (df['Ymax'] <2304-500) & (df["Xmin"]>500)
    df_cut=df[condition]

print("Events before IP selection:", len(df['X_ImpactPoint']),"Events after IP selection:", len(df_cut_temp['X_ImpactPoint']),"Events afer containement",len(df_cut['X_ImpactPoint']))
print("containment fraction after IP selection:", len(df_cut['X_ImpactPoint'])/len(df_cut_temp['X_ImpactPoint']))

adc,adc_error=df_cut["Integral"],1E-20*np.ones(len(df_cut["Integral"]))

#I have
# Calculate energy
energy = offset + adc * slope
# Propagate the error
energy_error = np.sqrt((offset_error)**2 + (adc * slope_error)**2 + (slope * adc_error)**2)

#DF CUT CONTAINS everything with cuts
outfile=args.output_file

#Realign the energy distribution with the intrinsic one onl if the intrinsic distribution file is passed
if args.intr_file is None:
    # Adding new columns
    df_cut['energy'] = energy
    df_cut['err_energy'] = energy_error
else:
    fileintr=args.intr_file
    # Open the ROOT file
    with uproot.open(fileintr) as file:
        # Access the TTree
        intrTree = file["angles"]
        angDeg_intr=intrTree["AngleDegree"].array(library="np")
        energy_intr=1000*(intrTree["EDep"].array(library="np"))
    #measured
    energy_meas=energy

    # Calculate histograms
    hist_intr, bins_intr = np.histogram(energy_intr, bins=100, density=True)
    hist_meas, bins_meas = np.histogram(energy_meas, bins=100, density=True)

    # Find the peak positions
    peak_intr=np.median(energy_intr)
    peak_meas=np.median(energy_meas)

    shift = peak_intr - peak_meas
    #print(peak_intr,peak_meas,shift)

    # Shift peak_meas
    energy_meas_shift = energy_meas + shift

    # Calculate new histogram for the shifted array2
    hist_meas_shifted, bins_meas_shift = np.histogram(energy_meas_shift, bins=100, density=True)

    # Plotting results
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.bar(bins_intr[:-1], hist_intr, width=np.diff(bins_intr), alpha=0.5, label='Intrisic Distribution')
    plt.bar(bins_meas[:-1], hist_meas, width=np.diff(bins_meas), alpha=0.5, label='Measured Distribution')
    plt.legend()
    plt.title('Before Shift')

    plt.subplot(1, 2, 2)
    plt.bar(bins_intr[:-1], hist_intr, width=np.diff(bins_intr), alpha=0.5, label='Intrisic Distribution')
    plt.bar(bins_meas_shift[:-1], hist_meas_shifted, width=np.diff(bins_meas_shift), alpha=0.5, label='Measured Distribution Shifted')
    plt.legend()
    plt.title('After Shift')

    plt.tight_layout()
    plt.show()

    # Adding new columns
    df_cut['energy'] = energy_meas_shift
    df_cut['err_energy'] = energy_error



#Angluar ditribution in degreees
# Example array of angles in radians
angles_radians = nparr(df_cut["AnglePCA"])
# Convert radians to degrees
angles_degrees = angles_radians * (180 / np.pi)
angle_shifted=np.empty(len(angles_degrees))
for i,angle in enumerate(angles_degrees):
    if angle>0:
        angle_shifted[i]=angle-180
    else:
        angle_shifted[i]=angle+180
df_cut["AngleDegree"]=angle_shifted


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

hist(nparr(df_cut["Integral"])/nparr(df_cut["ScSize"]),"delta (sc_int/sc_size)")

angDistr=hist(angle_shifted,"Angular distribution (\circ)",write=False)
# Define the combined fit function
gaussian=ROOT.TF1("gaussian", "gaus(0)",-180,180)
gausPol=ROOT.TF1("gausPol", "gaus(0)+pol0(3)",-180,180)
angDistr.Fit("gaussian","RQ")
gausPol.SetParameters(gaussian.GetParameter(0),gaussian.GetParameter(1),gaussian.GetParameter(2),0,0)
fitResult = angDistr.Fit(gausPol, "S")
angDistr.Write()



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

    # Initialize intersection points
    intersection_points = []

    # Check intersections with y = 0 and y = 2304, if they are within x-bounds
    if slope != 0:
        x_at_y0 = x - y / slope
        if 0 <= x_at_y0 <= 3000:
            intersection_points.append((x_at_y0, 0))

        x_at_y2304 = x + (2304 - y) / slope
        if 0 <= x_at_y2304 <= 3000:
            intersection_points.append((x_at_y2304, 2304))

    # Check intersections with x = 0 and x = 3000, if they are within y-bounds
    y_at_x0 = slope * (0 - x) + y
    if 0 <= y_at_x0 <= 2304:
        intersection_points.append((0, y_at_x0))

    y_at_x3000 = slope * (3000 - x) + y
    if 0 <= y_at_x3000 <= 2304:
        intersection_points.append((3000, y_at_x3000))

    # Determine the best valid intersection points
    if len(intersection_points) >= 2:
        # Sort points by proximity to the original point (x, y)
        intersection_points.sort(key=lambda p: math.hypot(p[0] - x, p[1] - y))
        start_point = intersection_points[0]
        end_point = intersection_points[-1]
    else:
        # Fallback to the original point if no valid intersections
        start_point = (x, y)
        end_point = (x, y)

    # Create and store the line and marker
    lines[f'line{i}'] = ROOT.TLine(start_point[0], start_point[1], end_point[0], end_point[1])
    lines[f'line{i}'].Draw("same")

    markers[f'marker{i}'] = ROOT.TMarker(x, y, 20)
    markers[f'marker{i}'].SetMarkerColor(ROOT.kRed)
    markers[f'marker{i}'].Draw("same")

# Update the canvas to display all lines and axes
canvas.Update()
# Optional: Save the canvas as an image
canvas.SaveAs("lines_at_angles_with_axes.png")
canvas.Write()