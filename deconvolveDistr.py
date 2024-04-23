import numpy as np
from scipy.fft import fft, ifft
import matplotlib.pyplot as plt
import argparse
import uproot
import ROOT
from scipy.signal import wiener
from tqdm import tqdm
import os

parser = argparse.ArgumentParser(description="Deconvolve f=g*h to get g.")
parser.add_argument("intrinsic", help="root file from the Geant4 intrinsics spread")
parser.add_argument("measured", help="root file with measured directions")
parser.add_argument("output_file", help="Name of the output file")
parser.add_argument("-d","--deconvolution",type=int, help="Deconvolution method 0 is Richardson-Lucy 1 is Fourier-Wiener, 2 is Tikhonov",default=0)
parser.add_argument("-t","--test", help="Using gaussian instead of real distr to test algorithms",default=False, action='store_true')
parser.add_argument("-v","--verbose", help="Set Verbose",default=False, action='store_true')
args = parser.parse_args()

def compute_mean_and_std(histogram, bin_edges):
    # Calculate midpoints of bins
    bin_midpoints = (bin_edges[:-1] + bin_edges[1:]) / 2
    
    # Ensure the histogram sums to 1 (if not already)
    histogram_normalized = histogram / np.sum(histogram)
    
    # Calculate mean
    mean = np.sum(bin_midpoints * histogram_normalized)
    
    # Calculate variance
    variance = np.sum((bin_midpoints - mean) ** 2 * histogram_normalized)
    
    # Standard deviation is the square root of variance
    std_dev = np.sqrt(variance)
    
    return mean, std_dev
def hist(data, x_name, channels=100, linecolor=4, linewidth=4,write=True, normalize=False):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, 0.99*np.min(array), 1.01*np.max(array))
    # Use numpy vectorization to fill histogram
    for x in array:
        hist.Fill(x)
    # Normalize if required
    if normalize:
        if hist.Integral() > 0:  # Check to avoid division by zero
            hist.Scale(1.0 / hist.Integral())
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
def hist_with_canvas(data, x_name, channels=200, linecolor=4, linewidth=4, write=True, normalize=False):
    # Convert list directly to numpy array to avoid redundant loop
    array = np.array(data, dtype="d")
    # Create histogram
    hist = ROOT.TH1D(x_name, x_name, channels, -180, 180)
    # Use numpy vectorization to fill histogram
    for x in array:
        hist.Fill(x)
    # Normalize if required
    if normalize:
        if hist.Integral() > 0:  # Check to avoid division by zero
            hist.Scale(1.0 / hist.Integral())
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
    
    # Create canvas
    canvas = ROOT.TCanvas(x_name, x_name, 1000, 1000)
    # Draw histogram on canvas
    hist.Draw()
    # Update canvas
    canvas.Update()
    canvas.SaveAs(f"{x_name}.png")
    # Return canvas
    return canvas
def richardson_lucy(histogram, psf, iterations):
    # histogram is the measured distribution f
    # psf is the point spread function (intrinsic distribution h)
    # iterations is the number of iterations
    rl_estimate = np.copy(histogram)
    for i in range(iterations):
        relative_blur = histogram / np.convolve(rl_estimate, psf, mode='same')
        rl_estimate *= np.convolve(relative_blur, psf[::-1], mode='same')
    return rl_estimate
def tikhonov_deconvolution(histogram, psf, regularization_param):
    F = fft(histogram)
    H = fft(psf)
    # Tikhonov regularization
    G = F / (H + regularization_param)
    return ifft(G).real
def bootstrap_deconvolution(angDeg_meas, angDeg_intr, bin_edges, num_bootstrap, deconv_method=args.deconvolution):
    means = []
    stds = []
    for _ in range(num_bootstrap):
        # Resample the original data
        resampled_meas = np.random.choice(angDeg_meas, size=len(angDeg_meas), replace=True)
        resampled_intr = np.random.choice(angDeg_intr, size=len(angDeg_intr), replace=True)

        # Histogram the resampled data
        hist_meas, _ = np.histogram(resampled_meas, bins=bin_edges, density=True)
        hist_intr, _ = np.histogram(resampled_intr, bins=bin_edges, density=True)

        # Deconvolve based on the chosen method
        if deconv_method == 0:  # Richardson-Lucy
            deconvolved = richardson_lucy(hist_meas, hist_intr, 50)
        elif deconv_method == 1:  # Fourier-Wiener
            F = fft(hist_meas)
            H = fft(hist_intr)
            epsilon = 1e-20
            H[H == 0] = epsilon
            G = F / H
            deconvolved = ifft(G).real
            deconvolved = wiener(deconvolved, noise=0.1)
        elif deconv_method == 2:  # Tikhonov
            regularization_param = 0.01
            deconvolved = tikhonov_deconvolution(hist_meas, hist_intr, regularization_param)

        # Calculate mean and std of the deconvolved result
        mean, std = compute_mean_and_std(deconvolved, bin_edges)
        means.append(mean)
        stds.append(std)

    return means, stds
def deconvolution(angDeg_meas, angDeg_intr, bin_edges, deconv_method=args.deconvolution):
    # Histogram the data
    hist_meas, _ = np.histogram(angDeg_meas, bins=bin_edges, density=True)
    hist_intr, _ = np.histogram(angDeg_intr, bins=bin_edges, density=True)

    # Deconvolve based on the chosen method
    if deconv_method == 0:  # Richardson-Lucy
        deconvolved = richardson_lucy(hist_meas, hist_intr, 100)
    elif deconv_method == 1:  # Fourier-Wiener
        F = fft(hist_meas)
        H = fft(hist_intr)
        epsilon = 1e-20
        H[H == 0] = epsilon
        G = F / H
        deconvolved = ifft(G).real
        deconvolved = wiener(deconvolved, noise=0.1)
    elif deconv_method == 2:  # Tikhonov
        regularization_param = 0.01
        deconvolved = tikhonov_deconvolution(hist_meas, hist_intr, regularization_param)
    """ 
    # Plot the results
    plt.figure(figsize=(10, 10))
    plt.subplot(3, 1, 1)
    plt.plot(bin_edges[:-1], hist_meas, label='Measured Distribution f')
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(bin_edges[:-1], hist_intr, label='Intrinsic Distribution h')
    plt.legend()

    plt.subplot(3, 1, 3)
    if deconv_method == 0:  # Richardson-Lucy
        plt.plot(bin_edges[:-1], deconvolved, label='Deconvolved Distribution g Richardson-Lucy Deconvolution')
    elif deconv_method == 1:  # Fourier-Wiener
        plt.plot(bin_edges[:-1], deconvolved, label='Deconvolved Distribution g Fourier-Wiener Deconvolution')
    elif deconv_method == 2:  # Tikhonov
        plt.plot(bin_edges[:-1], deconvolved, label='Deconvolved Distribution g Tikhonov Deconvolution')
    plt.legend()
    plt.tight_layout()
    plt.show()
    """
    mean_dec, std_dev_dec = compute_mean_and_std(deconvolved, bin_edges)

    return mean_dec, std_dev_dec, deconvolved
def write_numpy_to_root(np_histogram, bin_edges, hist_name, linecolor=4, linewidth=4,write=True):
    """
    Writes a numpy histogram array to a ROOT TH1D histogram and saves it to a ROOT file.
    """
    # Calculate the number of bins (subtract one because bin_edges has an extra edge)
    num_bins = len(bin_edges) - 1

    # Create the ROOT histogram
    hist = ROOT.TH1D(hist_name, hist_name, num_bins, bin_edges)

    # Fill the ROOT histogram
    for i in range(num_bins):
        # Set the bin content
        # NumPy bin content is assumed to be the integral count in the bin,
        # so for ROOT, which expects the integral count, we use the value directly.
        hist.SetBinContent(i + 1, np_histogram[i])

    # Set visual attributes and axis titles
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(hist_name)
    hist.GetYaxis().SetTitle("Entries")
    # Set maximum digits on axes to manage display
    hist.GetYaxis().SetMaxDigits(3)
    hist.GetXaxis().SetMaxDigits(3)
    if write:
        hist.Write()
    return hist
def grapherrAsym(x, y, ex, ey_low, ey_high, x_string, y_string, name=None, color=4, markerstyle=22, markersize=2, write=True):
    # Convert input data to numpy arrays of type double
    x_array = np.array(x, dtype="d")
    y_array = np.array(y, dtype="d")
    ex_array = np.array(ex, dtype="d")
    ey_low_array = np.array(ey_low, dtype="d")
    ey_high_array = np.array(ey_high, dtype="d")

    # Create a TGraphAsymmErrors object
    plot = ROOT.TGraphAsymmErrors(len(x), x_array, y_array, ex_array, ex_array, ey_low_array, ey_high_array)

    # Set the title and axis labels
    if name is None:
        plot.SetNameTitle(y_string + " vs " + x_string, y_string + " vs " + x_string)
    else:
        plot.SetNameTitle(name, name)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)

    # Set marker properties
    plot.SetMarkerColor(color)  # Set the color (default is blue)
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)

    # Write the plot to the current directory or file if specified
    if write:
        plot.Write()

    return plot
def draw_on_canvas(graph, canvas_size=(1000, 1000), save_path=None):
    """
    Draw a graph on a ROOT canvas and optionally save it to a file.

    Parameters:
        graph (ROOT.TGraphAsymmErrors): The graph to be drawn.
        title (str): Title of the canvas.
        canvas_size (tuple): Size of the canvas as (width, height).
        save_path (str): If provided, the path to save the canvas as an image file.

    Returns:
        ROOT.TCanvas: The canvas with the graph drawn.
    """
    # Create a canvas
    canvas = ROOT.TCanvas("canvas", graph.GetTitle(), canvas_size[0], canvas_size[1])
    #canvas.SetGrid()  # Enable grid on the canvas

    # Draw the graph
    graph.Draw("AP")  # Draw the graph with axes and points
    canvas.Update()  # Update the canvas to reflect changes

    # Set axis titles style
    graph.GetXaxis().SetTitleSize(0.04)
    graph.GetXaxis().SetLabelSize(0.035)
    graph.GetYaxis().SetTitleSize(0.04)
    graph.GetYaxis().SetLabelSize(0.035)

    # Optionally save the canvas to a file
    if save_path:
        canvas.SaveAs(save_path)

    return canvas
def GetStdErr(arr):
    mean=np.mean(arr)
    std=np.std(arr)
    N=len(arr)
    D4=np.sum(np.square(np.square(arr-mean)))/N
    #print(D4)
    return np.sqrt( (D4-std**4)/(N-1) )/(2*std)
def propagate_error_std(s_meas, s_intr, sigma_meas, sigma_intr):
    f = np.sqrt(s_meas**2 - s_intr**2)
    if s_meas**2 <= s_intr**2:
        raise ValueError("s_meas^2 should be larger than s_intr^2 to have a real result")

    df_ds_meas = s_meas / f
    df_ds_intr = -s_intr / f

    sigma_f = np.sqrt((df_ds_meas**2 * sigma_meas**2) + (df_ds_intr**2 * sigma_intr**2))
    return f, sigma_f
def ensure_directory_exists(directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path):
        # Create the directory, also creating any necessary intermediate directories
        os.makedirs(directory_path, exist_ok=True)
        #print(f"Directory '{directory_path}' was created.")
def draw_multigraph_with_legend(name, graphs, graph_labels,y_title,x_title,out_dir="./"):
    # Ensure each canvas has a unique name to avoid conflicts
    canvas_name = name
    canvas = ROOT.TCanvas(canvas_name, canvas_name, 1000, 1000)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)

    multigraph = ROOT.TMultiGraph()
    legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.95)  # Adjust the position as needed

    # Define colors and markers (add more if you have many graphs)
    colors = [ROOT.kBlue+1, ROOT.kRed+1, ROOT.kOrange+7, ROOT.kSpring-1, ROOT.kPink+1, ROOT.kTeal-1]
    markers = [20, 21, 22, 23, 29,33]

    for i, (graph, label) in enumerate(zip(graphs, graph_labels)):
        color = colors[i % len(colors)]
        marker = markers[i % len(markers)]

        graph.SetMarkerColor(color)
        graph.SetLineColor(color)
        graph.SetMarkerStyle(marker)
        multigraph.Add(graph)

        # Add entry to the legend
        legend.AddEntry(graph, label, "lp")

    multigraph.Draw("AP")
    legend.Draw()  # Draw the legend

    multigraph.GetXaxis().SetTitle(x_title)
    multigraph.GetYaxis().SetTitle(y_title)
    multigraph.SetTitle(name)

    canvas.Draw()
    canvas.Write()
    ensure_directory_exists(out_dir)
    canvas.SaveAs(out_dir+"/"+canvas_name+".png")
    # Optional: return canvas and multigraph if you want to interact with them outside the function
    return canvas, multigraph

############################ IMPORT DATA ############################
filemeas=args.measured
# Open the ROOT file
with uproot.open(filemeas) as file:
    # Access the TTree
    measTree = file["mytree"]
    angDeg_meas=measTree["AngleDegree"].array(library="np")
    energy_meas=measTree["energy"].array(library="np")

fileintr=args.intrinsic
# Open the ROOT file
with uproot.open(fileintr) as file:
    # Access the TTree
    intrTree = file["angles"]
    angDeg_intr=intrTree["AngleDegree"].array(library="np")
    energy_intr=1000*(intrTree["EDep"].array(library="np"))
""" 
# Apply the cut with bitwise AND (&) and proper parentheses
meas_cut_indices = (energy_meas > Estarts[i]) & (energy_meas < Estops[i])
intr_cut_indices = (energy_intr > Estarts[i]) & (energy_intr < Estops[i])
# Filter the data arrays using the cut
angDeg_meas_filtered = angDeg_meas[meas_cut_indices]
angDeg_intr_filtered = angDeg_intr[intr_cut_indices]
energy_meas_filtered = energy_meas[meas_cut_indices]
energy_intr_filtered = energy_intr[intr_cut_indices]

 """
main = ROOT.TFile(args.output_file, "RECREATE")
main.mkdir("ImportedData")
main.cd("ImportedData")

hist(angDeg_meas,"Measured Angular distribution",normalize=True)
hist(energy_meas,"Measured Energy distribution",normalize=True)
hist(angDeg_intr,"Intrinsic Angular distribution",normalize=True,linecolor=2)
hist(energy_intr,"Intrinsic Energy distribution",normalize=True,linecolor=2)

m_meas, m_intr=np.mean(angDeg_meas),np.mean(angDeg_intr)
s_meas,s_intr=np.std(angDeg_meas),np.std(angDeg_intr)
theo_m,theo_s=m_meas+m_intr, np.sqrt(s_meas*s_meas-s_intr*s_intr)

# Generating some data for testing
if args.test is True:
    # Parameters
    n_meas = len(angDeg_meas)  # Number of points in measured distribution
    n_intr = len(angDeg_intr)  # Number of points in intrinsic distribution

    # Generate data
    angDeg_meas = np.random.normal(loc=m_meas, scale=s_meas, size=n_meas)
    angDeg_intr = np.random.normal(loc=m_intr, scale=s_intr, size=n_intr)

### TRY TO DECONVOLVE f=g*h I need to get g, f is the measured h is the intrinsic
# Define common bins
bin_edges = np.linspace(-100, 100, 151)  # 200 bins from -180 to 180

main.cd()
mean_dec, std_dev_dec,hist_dec=deconvolution(angDeg_meas, angDeg_intr, bin_edges)
write_numpy_to_root(hist_dec,bin_edges,"Deconvoluted Distribution")

########## BOOTSTRAP for Errors ############
# Assuming the rest of your setup and bin_edges are defined as before
num_bootstrap = 1000  # Number of bootstrap samples

main.mkdir("Bootstrap")
main.cd("Bootstrap")
bootstrap_means, bootstrap_stds = bootstrap_deconvolution(angDeg_meas, angDeg_intr, bin_edges, num_bootstrap)
hist(bootstrap_means,"Bootstrap Means")
hist(bootstrap_stds,"Bootstrap Stds")

# Calculate confidence intervals
mean_confidence_interval = np.percentile(bootstrap_means, [16, 84])
std_confidence_interval = np.percentile(bootstrap_stds, [16, 84])

if args.verbose is True:
    print("Mean confidence interval:", mean_confidence_interval)
    print("Standard deviation confidence interval:", std_confidence_interval)

    print("DECONVOLVED",mean_dec, std_dev_dec)
    print("THORETICHAL",theo_m,theo_s)
    print(f"Resolution @ 68% CL:",std_dev_dec,"-",std_dev_dec-std_confidence_interval[0],"+",std_dev_dec-std_confidence_interval[1])


##### SCAN IN ENERGY #####


Emean=[5,12.5,17.5,25,35]
Eerr=[5,2.5,2.5,5,5]
Estarts=[0,10,15,20,30]
Estops=[10,15,20,30,40]
print(Estarts,Estops)

""" 
Estarts=np.arange(0, 40, 2.5)
Estops=np.arange(2.5, 42.5, 2.5)
Emean=(Estops)-1.25
Eerr=1.25*np.ones(len(Emean))
"""
angRes,errAngresPLUS,errAngresMINUS=np.empty(len(Emean)),np.empty(len(Emean)),np.empty(len(Emean))
angRes_gaus,errAngres_gaus=np.empty(len(Emean)),np.empty(len(Emean))

bin_edges = np.linspace(-100, 100, 51)  # 200 bins from -180 to 180
for i,Estart in (enumerate(Estarts)):
    main.mkdir(f"Energy Scan cuts=[{Estarts[i]},{Estops[i]}]/Imported")
    main.cd(f"Energy Scan cuts=[{Estarts[i]},{Estops[i]}]/Imported")
    # Apply the cut with bitwise AND (&) and proper parentheses
    meas_cut_indices = (energy_meas > Estarts[i]) & (energy_meas < Estops[i])
    intr_cut_indices = (energy_intr > Estarts[i]) & (energy_intr < Estops[i])
    # Filter the data arrays using the cut
    angDeg_meas_filtered = angDeg_meas[meas_cut_indices]
    angDeg_intr_filtered = angDeg_intr[intr_cut_indices]
    energy_meas_filtered = energy_meas[meas_cut_indices]
    energy_intr_filtered = energy_intr[intr_cut_indices]
    hist(angDeg_meas_filtered,f"Measured Angular distribution cuts=[{Estarts[i]},{Estops[i]}]",normalize=True,channels=50)
    hist(energy_meas_filtered,f"Measured Energy distribution cuts=[{Estarts[i]},{Estops[i]}]",normalize=True,channels=50)
    hist(angDeg_intr_filtered,f"Intrinsic Angular distribution cuts=[{Estarts[i]},{Estops[i]}]",normalize=True,linecolor=2,channels=50)
    if args.verbose is True: hist_with_canvas(angDeg_intr_filtered,f"IntrinsicAngulardistributioncuts=[{Estarts[i]},{Estops[i]}]")
    hist(energy_intr_filtered,f"Intrinsic Energy distribution cuts=[{Estarts[i]},{Estops[i]}]",normalize=True,linecolor=2,channels=50)

    main.cd(f"Energy Scan cuts=[{Estarts[i]},{Estops[i]}]")
    mean_dec, std_dev_dec,hist_dec=deconvolution(angDeg_meas_filtered, angDeg_intr_filtered, bin_edges)
    write_numpy_to_root(hist_dec,bin_edges,"Deconvoluted Distribution cuts=[{Estarts[i]},{Estops[i]}]")

    main.mkdir(f"Energy Scan cuts=[{Estarts[i]},{Estops[i]}]/Bootstrap")
    main.cd(f"Energy Scan cuts=[{Estarts[i]},{Estops[i]}]/Bootstrap")
    bootstrap_means, bootstrap_stds = bootstrap_deconvolution(angDeg_meas_filtered, angDeg_intr_filtered, bin_edges, num_bootstrap)
    hist(bootstrap_stds,f"Bootstrap Stds cuts=[{Estarts[i]},{Estops[i]}]")
    # Calculate confidence intervals

    std_confidence_interval = np.percentile(bootstrap_stds, [16, 84])
    print(std_dev_dec,std_confidence_interval)
    angRes[i]=std_dev_dec
    errAngresPLUS[i]=std_dev_dec-std_confidence_interval[0]
    errAngresMINUS[i]=std_confidence_interval[1]-std_dev_dec
    print(Emean[i],angRes[i],errAngresPLUS[i],errAngresMINUS[i])

    s_meas,s_intr=np.std(angDeg_meas_filtered),np.std(angDeg_intr_filtered)
    es_meas,es_intr=GetStdErr(angDeg_meas_filtered),GetStdErr(angDeg_intr_filtered)

    angRes_gaus[i],errAngres_gaus[i]=propagate_error_std(s_meas,s_intr,es_meas,es_intr)

main.cd()
graphs=[grapherrAsym(Emean,angRes,Eerr,errAngresMINUS,errAngresPLUS,"e^{-} Energy (keV)","Angular resolution RMS (deg)",markersize=2,color=2)
        ,  grapherrAsym(Emean,angRes_gaus,Eerr,errAngres_gaus,errAngres_gaus,"e^{-} Energy (keV)","Angular resolution RMS (deg)",markersize=2)]
draw_multigraph_with_legend("Angluar Resolution",graphs,["Deconvoluted","Supposed Gaussian"],"Angular Resolution RMS (deg)","e^{-} Energy (keV)")
