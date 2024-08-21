#!/bin/bash
#python3 plot_multiGraph.py efficiency 'Efficiency plot' --log
python3 plot_multiGraph.py modulation 'Efficiency vs Photon Energy (keV)' --log
python3 plot_multiGraph.py modulation 'Modulation Factor vs Photon Energy (keV)' --legend_pos "bot_right"
python3 plot_multiGraph.py modulation 'figure-of-merit vs Photon Energy (keV)'
python3 plot_multiGraph.py deconvolved 'Measured Angular distribution'
python3 plot_multiGraph.py deconvolved 'Measured Energy distribution'
python3 plot_multiGraph.py deconvolved 'Intrinsic Angular distribution'
python3 plot_multiGraph.py deconvolved 'Intrinsic Energy distribution'
python3 plot_multiGraph.py deconvolved 'Angular resolution RMS (deg) vs e^{-} Energy (keV)'
python3 plot_multiGraph.py deconvolved 'Deconvolved Distribution' --histogram
