# Directionality

## To compile

```g++ Analyzer.cxx Base_script.cxx -o nameprog `root-config --libs --cflags` -lSpectrum```

## To run

```./nameprog path_to_rootfile [output_directory]```

### PIPELINE

- Reconstruct data
  - Data with Sr
  - data with Fe/Cd

- Run directionality Code on Sr only
  - Usually output a 'AfterDir'
  - It is better to run on Condor and use `splitROOT.py` to split the input file (usually 'reco')

- Run `calibration.py` on the Fe and Cd reco dataset
  - Usually output contain 'calib' in the name

- Run `DistrfromDirectionality.py` on directionality data with the calibration parameters
  - Usually takes in input an 'AfterDir' file and output a 'anal' file

- Run `deconvolveDistr.py` with Geant4 anlayzed data and the measured directionalities to get the deconvolded intrisic angular resolution
  - Usually contains 'deconvolved'
