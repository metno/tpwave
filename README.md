# tpwave
tpwave is a repository for an R program used for the Multiresolution analysis of precipitation based on discrete wavelet transformation. The repository contains also an R program for generating diagnostic figures.

## Required depdendencies

The following R libraries are required:
- waveslim (tested with version 1.8.2)
- ncdf4 (tested with version 1.17)
- argparser (tested with version 0.4)
- raster (tested with version 2.9-5)

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

The program used to post-process the output files, such as ```ERA5_tp_day_example_waveEn.RData``` or similar, is included as:

```
./tpwave_multires_ampfact.r
```

and it can also be used to perform the linear regression of the scale-dependent wavelet energies, to compute the inflation or deflation coefficients used to re-scale single weather events over climates that are different with respect to the one they have actually occurred.


## Copyright and license

Copyright © 2021 Norwegian Meteorological Institute. tpwave is licensed under The GNU Lesser General Public License (LGPL). See LICENSE file.

## Contact

E-mail: Cristian Lussana (cristianl@met.no)
