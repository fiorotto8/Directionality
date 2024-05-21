# Directionality pipeline to analyze data from INAF

- Reconstruct data
  - Data with Sr
  - data with Fe/Cd

- Run directionality Code on Sr only <https://github.com/fiorotto8/CygnoAnal>
  - Usually output a 'AfterDir'
  - It is better to run on Condor and use `splitROOT.py` to split the input file (usually 'reco')

- Run `calibration.py` on the Fe and Cd reco dataset
  - Usually output contain 'calib' in the name

- Run `DistrfromDirectionality.py` on directionality data with the calibration parameters
  - Usually takes in input an 'AfterDir' file and output a 'anal' file
  - Typical Impact point selection: `(df['X_ImpactPoint'] > 1700) & (df['X_ImpactPoint'] < 1800) & (df['Y_ImpactPoint'] > 850) & (df['Y_ImpactPoint'] < 1400)`
  - Typical containment selection: `(df['Ymin'] >500) & (df['Ymax'] <2304-500) & (df["Xmin"]>500)`

- Run `deconvolveDistr.py` with Geant4 anlayzed data and the measured directionalities to get the deconvoluted intrinsic angular resolution
  - *Needs Geant4 90Sr simulation* <https://github.com/fiorotto8/MANGO_RadioactiveSource>
  - Usually contains 'deconvolved'

- Run `modulation_factor.py` with deconvolved data and efficiency data
  - *Needs Geant4 gamma efficiency simulation* <https://github.com/fiorotto8/MANGOspaceG4>
  - Usually contains 'modulation'
  - Usually contains efficiency

- the script `createPlot.sh` is running the `plot_multiGraph.py` multiple times to create 'nice' canvas of some interesting quantities cycling over all the measurements
