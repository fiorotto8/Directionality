import uproot
import numpy as np
import ROOT
import argparse

parser = argparse.ArgumentParser(description="Split a ROOT TTree into multiple subsets.")
parser.add_argument("input_file", help="root file with the directionality analyzed data")
parser.add_argument("output_file", help="Name of the output file")
args = parser.parse_args()

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

# Open the original ROOT file
with uproot.open(args.input_file) as file:
    # Access the TTree
    tree = file["track_info"]
    # Read the data from the branch you're interested in
    X_IP = tree["X_ImpactPoint"].array(library="np")
    Y_IP = tree["Y_ImpactPoint"].array(library="np")

main = ROOT.TFile(args.output_file, "RECREATE")

# Check for NaN values
nan_indices_X = np.isnan(X_IP)
nan_indices_Y = np.isnan(Y_IP)
# Filter out NaN values
X_IP_filtered = X_IP[~nan_indices_X]
Y_IP_filtered = Y_IP[~nan_indices_Y]

IP=create_fill_TH2("IP positions","X","Y","Entries",X_IP_filtered,Y_IP_filtered,x_bins=20,y_bins=20)
