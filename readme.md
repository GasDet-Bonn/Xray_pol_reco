# Polarimetry reconstruction of Timepix and Timepix3 data
This script reconstructs photoelectron angles in a Timepix or
Timepix3 dataset via a two step approach. In the first step all
pixels are used to minimised the second moment of the charge
distribution. Based on the third moment of the distribution it
is then decided which pixels to keep for the second step, to
only use the start of the track. In the second step the second
moment of the charge distribution of the remaining pixels is again
minimised and then the angle in the xy plane with respect to the
x axis is calculated.

## Requirements
The following python packages are required to use this reconstruction:
```
numpy
tqdm
h5py
```

## Usage
The script can be used by:
```
python3 reco.py <path_to_hdf5>
```
The `path_to_hdf5` is a path to an Hdf5 file that was already processed
with the reconstruction of [TimepixAnalysis](https://github.com/Vindaar/TimepixAnalysis). It stores the results in the same HDF5 file.