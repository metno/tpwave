# tpwave
tpwave is a repository for an R program used for the Multiresolution analysis of precipitation based on discrete wavelet transformation. The repository contains also an R program for generating diagnostic figures.

## Required depdendencies

The following R libraries are required:
- waveslim
- ncdf4
- waveslim
- argparser
- raster

## Examples

Multiresolution analysis of ERA5 daily precipitation based on 2D Haar wavelets:

```
./tpwave.r
```

It produces the output file ```ERA5_tp_day_example_waveEn.RData``` based on the input data ```./data/ERA5_tp_day_1950-01-01_1950-01-02.nc```.
The format of the otput RData file is described in ```tpwave.r``` (sec. write output file).

Production of several figures based on the RData files stored in ```./data/```

```
./plot_tpwave.r --filter_len 30
```




## Copyright and license

Copyright © 2021 Norwegian Meteorological Institute. tpwave is licensed under The GNU Lesser General Public License (LGPL). See LICENSE file.

## Contact

E-mail: Cristian Lussana (cristianl@met.no)
