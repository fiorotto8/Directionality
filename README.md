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
- Run `calibration.py` on the Fe and Cd reco dataset
- Run `DistrfromDirectionality.py` on directionality data with the calibration parameters
