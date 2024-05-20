import ROOT
import pandas as pd
import uproot
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="Plot a certain quantity for all the available mixtures")
parser.add_argument("--folder", help="offset forlder",default="./LongRun/")
parser.add_argument("root_name", help="root file name where the plot may be find",type=str)
parser.add_argument("plot_name", help="Name of the TGraph to fetch",type=str)
parser.add_argument("--log", help="enable log y scale",action='store_true')
parser.add_argument("--histogram", help="Draw the TGraph in histogram style", action='store_true')
parser.add_argument("--legend_pos", help="Position of the legend: top_right, top_left, bot_right, bot_left", type=str, default="top_right")
parser.add_argument("--v", help="enable verbose",action='store_true')
args = parser.parse_args()

def find_root_files_with_subfolders(directory,file_name):
    root_files = []
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.startswith(file_name) and filename.endswith(".root"):
                subfolder = os.path.relpath(dirpath, directory)
                if subfolder == '.':
                    subfolder = ''  # If the file is in the root directory, subfolder is an empty string
                root_files.append((os.path.join(dirpath, filename), subfolder))
    return root_files
def fetch_tgraph_from_file(file_path, graph_name):
    try:
        root_file = ROOT.TFile.Open(file_path)
        if root_file and not root_file.IsZombie():
            tgraph = root_file.Get(graph_name)
            if tgraph:
                return tgraph
            else:
                print(f"Graph {graph_name} not found in file {file_path}")
        else:
            print(f"Failed to open file {file_path}")
    except Exception as e:
        print(f"Error accessing {file_path}: {e}")
    finally:
        if root_file:
            root_file.Close()
    return None
def draw_multigraph_with_legend(tgraphs, subfolders, title, legend_position):
    if not tgraphs:
        print("No TGraphs to plot.")
        return
    # Define a set of easily distinguishable colors and marker styles
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kBlack, ROOT.kViolet, ROOT.kPink, ROOT.kYellow]
    marker_styles = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]

    # Access the first TGraph to get axis titles
    first_tgraph = tgraphs[0]
    x_axis_title = first_tgraph.GetXaxis().GetTitle()
    y_axis_title = first_tgraph.GetYaxis().GetTitle()

    canvas = ROOT.TCanvas("canvas", "Multigraph", 1000, 1000)
    multigraph = ROOT.TMultiGraph()
    # Determine legend position
    if legend_position == "top_right":
        legend = ROOT.TLegend(0.7, 0.75, 0.95, 0.95)
    elif legend_position == "top_left":
        legend = ROOT.TLegend(0.05, 0.7, 0.3, 0.9)
    elif legend_position == "bot_right":
        legend = ROOT.TLegend(0.7, 0.1, 0.95, 0.3)
    elif legend_position == "bot_left":
        legend = ROOT.TLegend(0.05, 0.05, 0.3, 0.3)
    else:
        print("Invalid legend position. Defaulting to top_right.")
        legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.9)

    for i, (tgraph, subfolder) in enumerate(zip(tgraphs, subfolders)):
        if tgraph:
            color = colors[i % len(colors)]
            marker_style = marker_styles[i % len(marker_styles)]
            tgraph.SetMarkerColor(color)
            tgraph.SetLineColor(color)
            tgraph.SetMarkerStyle(marker_style)
            tgraph.SetMarkerSize(2)

            multigraph.Add(tgraph)  # Add the TGraph to the TMultiGraph
            legend.AddEntry(tgraph, subfolder, "lp")  # Add the TGraph to the legend with subfolder name

    if args.log is True: canvas.SetLogy()

    multigraph.SetTitle(title)
    multigraph.GetXaxis().SetTitle(x_axis_title)  # Set the X axis title
    multigraph.GetYaxis().SetTitle(y_axis_title)  # Set the Y axis title
    multigraph.Draw("AP")
    legend.Draw()

    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.Update()
    canvas.SaveAs(f"{title.replace(' ', '_')}.png")  # Save the canvas as an image file
    #canvas.SaveAs("multigraph.pdf")  # Save the canvas as a PDF file
    #canvas.SaveAs("multigraph.root")  # Save the canvas as a ROOT file

def draw_multigraph_with_legend_hist(tgraphs, subfolders, title, legend_position):
    if not tgraphs:
        print("No TGraphs to plot.")
        return

    # Define a set of easily distinguishable colors and marker styles
    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kOrange, ROOT.kBlack, ROOT.kViolet, ROOT.kPink, ROOT.kYellow]
    marker_styles = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29]

    # Access the first TGraph to get axis titles
    first_tgraph = tgraphs[0]
    x_axis_title = first_tgraph.GetXaxis().GetTitle()
    y_axis_title = first_tgraph.GetYaxis().GetTitle()

    canvas = ROOT.TCanvas("canvas", "Multigraph", 1000, 1000)
    # Determine legend position
    if legend_position == "top_right":
        legend = ROOT.TLegend(0.7, 0.75, 0.95, 0.95)
    elif legend_position == "top_left":
        legend = ROOT.TLegend(0.05, 0.7, 0.3, 0.9)
    elif legend_position == "bot_right":
        legend = ROOT.TLegend(0.7, 0.1, 0.95, 0.3)
    elif legend_position == "bot_left":
        legend = ROOT.TLegend(0.05, 0.05, 0.3, 0.3)
    else:
        print("Invalid legend position. Defaulting to top_right.")
        legend = ROOT.TLegend(0.7, 0.7, 0.95, 0.9)

    histograms=[]
    for i, (tgraph, subfolder) in enumerate(zip(tgraphs, subfolders)):
        if tgraph:
            n_points = tgraph.GetN()
            x_values = np.array(tgraph.GetX())
            y_values = np.array(tgraph.GetY())
            x_errors_low = np.array(tgraph.GetEXlow()) if tgraph.GetEXlow() else np.zeros(n_points)
            x_errors_high = np.array(tgraph.GetEXhigh()) if tgraph.GetEXhigh() else np.zeros(n_points)
            y_errors_low = np.array(tgraph.GetEYlow()) if tgraph.GetEYlow() else np.zeros(n_points)
            y_errors_high = np.array(tgraph.GetEYhigh()) if tgraph.GetEYhigh() else np.zeros(n_points)

            histogram = ROOT.TH1F(f"hist_{i}", "", n_points, x_values[0] - x_errors_low[0], x_values[-1] + x_errors_high[-1])
            histogram.SetLineColor(colors[i % len(colors)])
            histogram.SetLineWidth(2)
            histogram.SetFillColor(0)  # Set the fill color to 0, which means no fill
            histogram.SetFillStyle(0)  # Set the fill style to 0, which means no fill

            for j in range(n_points):
                bin_index = histogram.FindBin(x_values[j])
                bin_width = (x_errors_low[j] + x_errors_high[j]) / 2
                histogram.SetBinContent(bin_index, y_values[j])
                histogram.SetBinError(bin_index, (y_errors_low[j] + y_errors_high[j]) / 2)

            histograms.append(histogram)
            legend.AddEntry(histogram, subfolder, "f")  # Add the histogram to the legend with subfolder name

    if args.log:
        canvas.SetLogy()

    for histogram in histograms:
        histogram.Draw("E3 SAME")  # Draw the histogram with error bars and fill

    first_histogram = histograms[0]
    first_histogram.GetXaxis().SetTitle(x_axis_title)  # Set the X axis title
    first_histogram.GetYaxis().SetTitle(y_axis_title)  # Set the Y axis title

    legend.Draw()

    canvas.SetTitle(title)
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.05)
    canvas.SetTopMargin(0.05)
    canvas.Update()
    canvas.SaveAs(f"{title.replace(' ', '_')}.png")  # Save the canvas as an image file


# Example usage
tgraphs,subfolders=[],[]
directory_to_scan = args.folder
root_files = find_root_files_with_subfolders(directory_to_scan, args.root_name)
for file, subfolder in root_files:
    if args.v is True: print(f"File: {file}, Subfolder: {subfolder}")
    tgraph = fetch_tgraph_from_file(file, args.plot_name)
    if tgraph:
        if args.v is True: print(f"Successfully fetched {args.plot_name} from {file}")
        tgraphs.append(tgraph)
        subfolders.append(subfolder)

# Draw the multigraph with legend
if args.histogram is False: draw_multigraph_with_legend(tgraphs, subfolders, f"{args.plot_name}", args.legend_pos)
else: draw_multigraph_with_legend_hist(tgraphs, subfolders, f"{args.plot_name}", args.legend_pos)