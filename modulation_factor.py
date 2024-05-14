import numpy as np
import argparse
import uproot
import ROOT
from tqdm import tqdm
import os
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.integrate import quad

parser = argparse.ArgumentParser(description="Obtain modulation factor from the intrinsic angular resolution.")
parser.add_argument("deconvoluted", help="root file from the deconvultion")
parser.add_argument("outfile",help="Name of the output root file")
args = parser.parse_args()

# Define the cosine squared function
def cos_squared(x):
    return np.cos(np.radians(x))**2
# Define a function to create a Gaussian
def gaussian(x, mean, sigma):
    return norm.pdf(x, mean, sigma)
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

main = ROOT.TFile(args.outfile, "RECREATE")

modulation_factors=np.empty(len(energy))
modulation_errors=np.empty(len(energy))
graphs=[]
# Example usage
x_values = np.linspace(-180, 180, 1000)
num_samples=100
# Example usage
for i in tqdm(range(len(energy))):
    mean = 0
    sigma = angRes[i]

    # Generate multiple samples of angRes with noise
    convolved_samples = []
    for _ in range(num_samples):
        sampled_sigma = np.random.normal(sigma, errAngres[i])
        gauss = gaussian(x_values, mean, sampled_sigma)
        cos2 = cos_squared(x_values)
        convolved_values = np.convolve(gauss, cos2, mode='same')
        convolved_samples.append(convolved_values)

    # Compute the mean and standard deviation of the convolved results
    convolved_samples = np.array(convolved_samples)
    convolved_mean = np.mean(convolved_samples, axis=0)
    convolved_std = np.std(convolved_samples, axis=0)
    convolved_sem = convolved_std / np.sqrt(num_samples)

    # Calculate modulation factor and its error using standard error propagation
    A = np.max(convolved_mean)
    B = np.min(convolved_mean)
    # Errors on A and B (assuming they are convolved_sem)
    sigma_A = convolved_sem[np.argmax(convolved_mean)]
    sigma_B = convolved_sem[np.argmin(convolved_mean)]
    # Modulation factor
    modulation_factor = (A - B) / (A + B)
    # Error propagation
    partial_A = 2 * B / (A + B)**2
    partial_B = -2 * A / (A + B)**2
    sigma_mod = np.sqrt((partial_A * sigma_A)**2 + (partial_B * sigma_B)**2)

    modulation_factors[i]=(modulation_factor)
    modulation_errors[i]=(sigma_mod)

    graphs.append(grapherr(x_values,convolved_mean,1E-20*np.ones(len(x_values)),convolved_sem,"Angle (deg)",'Convolution Result' ))

    # Plotting the convolved function with error bars
    plt.plot(x_values, convolved_mean, label=f'{energy[i]} keV')
    plt.fill_between(x_values, convolved_mean - convolved_sem, convolved_mean + convolved_sem, alpha=0.3)

plt.xlabel('Angle (deg)')
plt.ylabel('Convolution Result')
plt.title('Convolution of Gaussian with cos^2 Function')
plt.legend()
plt.savefig('Cos2Conv.png', bbox_inches='tight', dpi=1000)
#plt.show()

grapherr(energy,modulation_factors,e_energy,modulation_errors,"Photon Energy (keV)", "Modulation Factor")