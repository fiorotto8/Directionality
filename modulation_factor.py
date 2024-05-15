import numpy as np
import argparse
import uproot
import ROOT
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad
import pandas as pd
import math as m

parser = argparse.ArgumentParser(description="Obtain modulation factor from the intrinsic angular resolution.")
parser.add_argument("deconvoluted", help="root file from the deconvultion")
parser.add_argument("efficiency", help="root file from the efficiency")
parser.add_argument("outfile",help="Name of the output root file")
args = parser.parse_args()

# Define the cosine squared function
def cos_squared(x):
    return np.cos(np.radians(x))**2
# Define a function to create a Gaussian
def gaussian(x, mean, sigma):
    return norm.pdf(x, mean, sigma)
def gaussian_cos2_convolution(x_array, sigma):
    """
    Calculate the convolution of a Gaussian with a cos^2 function for an array of x values.

    Parameters:
    x_array (numpy array): Array of x values.
    sigma (float): Standard deviation of the Gaussian function.

    Returns:
    numpy array: The convolution result for each x value.
    """
    # Constants
    coeff = (sigma * np.sqrt(2 * np.pi)) / 2
    exp_factor = np.exp(-2 * sigma**2)
    
    # Compute the convolution for each x
    result = coeff * (1 + exp_factor * np.cos(2 * x_array))
    
    return result
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
def style_and_draw_graph(graph, canvas_name="Modulation Factor",color=ROOT.kBlue ,log=False):
    # Apply styles to the graph
    graph.SetNameTitle(canvas_name, "")
    graph.GetXaxis().SetTitle("Energy (keV)")
    graph.GetYaxis().SetTitle(canvas_name)
    graph.SetMarkerColor(color)  # blue
    graph.SetLineColor(color)
    graph.SetLineWidth(3)
    graph.SetMarkerStyle(8)
    graph.SetMarkerSize(2)
    # Create a canvas
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 1000, 1000)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    # Set log scale on the y-axis
    if log==True: canvas.SetLogy()
    # Draw the graph on the canvas
    graph.Draw("APZ")  # "A" for axis, "P" for points
    # Set the y-axis range
    #graph.GetYaxis().SetRangeUser(1e-5, 1e-1)
    # Update the canvas to display the graph
    canvas.Update()
    # Save the canvas if needed
    canvas.SaveAs(f"{canvas_name}.png")
    canvas.Write()
    return canvas
def analytical_gaussian_cos2_convolution(x_array_degrees, sigma_degrees, sigma_sigma_degrees):
    """
    Analytical convolution of a Gaussian with a cos^2 function for an array of x values in degrees,
    including propagated error from the statistical error of sigma.

    Parameters:
    x_array_degrees (numpy array): Array of x values (angles in degrees).
    sigma_degrees (float): Standard deviation of the Gaussian function in degrees.
    sigma_sigma_degrees (float): Statistical error of the standard deviation of the Gaussian in degrees.

    Returns:
    tuple: (convolved values, propagated errors)
    """
    # Convert degrees to radians for the cosine calculation
    x_array_radians = np.radians(x_array_degrees)
    sigma_radians = np.radians(sigma_degrees)
    sigma_sigma_radians = np.radians(sigma_sigma_degrees)
    exp_factor = np.exp(-2 * sigma_radians**2)
    convolved_values = 0.5 * (1 + exp_factor * np.cos(2 * x_array_radians))
    # Partial derivative of the convolved function with respect to sigma
    partial_sigma = -2 * sigma_radians * exp_factor * np.cos(2 * x_array_radians)
    # Propagated errors
    propagated_errors = np.abs(partial_sigma) * sigma_sigma_radians

    return convolved_values, propagated_errors


# Open the ROOT file and access the TTree
file = uproot.open(args.deconvoluted)
tree = file['Resolutions']

# Read the branches into NumPy arrays
energy = tree['energy'].array()
e_energy = tree['e_energy'].array()
angRes = tree['angRes'].array()
errAngres = (tree['errAngresPLUS'].array()+tree['errAngresMINUS'].array())/2 #takes the average of the two since we suppose it is gaussian
angRes_gaus = tree['angRes_gaus'].array()
errAngres_gaus = tree['errAngres_gaus'].array()

modulation_factors=np.empty(len(energy))
modulation_errors=np.empty(len(energy))
graphs=[]
# Example usage
x_values = np.linspace(-360, 360, 1000)
num_samples=1000
# Example usage
for i in tqdm(range(len(energy))):
    mean = 0
    sigma,e_sigma = angRes[i],errAngres[i]

    convolved_mean,convolved_err=analytical_gaussian_cos2_convolution(x_values,sigma,e_sigma)

    # Calculate modulation factor and its error using standard error propagation
    A = np.max(convolved_mean)
    B = np.min(convolved_mean)
    # Errors on A and B (assuming they are convolved_sem)
    sigma_A = convolved_err[np.argmax(convolved_mean)]
    sigma_B = convolved_err[np.argmin(convolved_mean)]
    # Modulation factor
    modulation_factor = (A - B) / (A + B)
    # Error propagation
    partial_A = 2 * B / (A + B)**2
    partial_B = -2 * A / (A + B)**2
    sigma_mod = np.sqrt((partial_A * sigma_A)**2 + (partial_B * sigma_B)**2)

    modulation_factors[i]=(modulation_factor)
    modulation_errors[i]=(sigma_mod)

    graphs.append(grapherr(x_values,convolved_mean,1E-20*np.ones(len(x_values)),convolved_err,"Angle (deg)",'Convolution Result',write=False ))

    # Plotting the convolved function with error bars
    plt.plot(x_values, convolved_mean, label=f'{energy[i]} keV')
    plt.fill_between(x_values, convolved_mean - convolved_err, convolved_mean + convolved_err, alpha=0.3)

plt.xlabel('Angle (deg)')
plt.ylabel('Convolution Result')
plt.title('Convolution of Gaussian with cos^2 Function')
plt.legend()
# Set the range of the x-axis
plt.xlim(-180, 180)
plt.savefig('Cos2Conv.png', bbox_inches='tight', dpi=1000)
#plt.show()

# Now the quality factor
# Open the ROOT file and access the TTree
file_eff = uproot.open(args.efficiency)
tree_eff = file_eff["efficiency"]
# Convert the TTree to a DataFrame
df_eff = tree_eff.arrays(library="pd")

# List to store results
averaged_efficiencies = []
err_averaged_efficiencies = []

for en, e_en in zip(energy, e_energy):
    # Define the range
    lower_bound = en - e_en
    upper_bound = en + e_en

    # Filter the DataFrame
    filtered_df = df_eff[(df_eff['energy'] >= lower_bound) & (df_eff['energy'] <= upper_bound)]

    # Calculate the average efficiency
    if not filtered_df.empty:
        average_efficiency = filtered_df['efficiency'].mean()
        err_average_efficiency = filtered_df['efficiency'].std()/np.sqrt(len(filtered_df['efficiency']))
    else:
        average_efficiency = np.nan  # or any other value indicating no data in range

    # Store the result
    averaged_efficiencies.append(average_efficiency)
    err_averaged_efficiencies.append(err_average_efficiency)

# Create a DataFrame
data = {
    'energy': energy,
    'modulation_factors': modulation_factors,
    'e_energy': e_energy,
    'modulation_errors': modulation_errors
}
df = pd.DataFrame(data)
# Convert DataFrame to a dictionary of numpy arrays for uproot
arrays = {col: df[col].values for col in df.columns}

# Write the DataFrame to a ROOT file
with uproot.recreate(args.outfile) as file:
    file["tree"] = arrays


main = ROOT.TFile(args.outfile, "UPDATE")
for g in graphs: g.Write()
mod=grapherr(energy,modulation_factors,e_energy,modulation_errors,"Photon Energy (keV)", "Modulation Factor")
style_and_draw_graph(mod,canvas_name="Modulation Factor",color=ROOT.kBlue)

eff=grapherr(energy,averaged_efficiencies,e_energy,err_averaged_efficiencies,"Photon Energy (keV)", "Efficiency")
style_and_draw_graph(eff,canvas_name="Efficiency",log=True,color=ROOT.kViolet-1)

fom=modulation_factors*np.sqrt(averaged_efficiencies)
#fom=modulation_factors*(averaged_efficiencies)
e_fom=np.sqrt(
    (np.sqrt(averaged_efficiencies) * modulation_errors)**2 +
    (modulation_factors / (2 * np.sqrt(averaged_efficiencies)) * err_averaged_efficiencies)**2
)

mod_qual=grapherr(energy,fom,e_energy,e_fom,"Photon Energy (keV)", "figure-of-merit")
style_and_draw_graph(mod_qual,canvas_name="figure-of-merit",color=ROOT.kRed)


