# tpwave
tpwave is a repository for an R program used for the Multiresolution analysis of precipitation based on discrete wavelet transformation. The repository contains also an R program for generating diagnostic figures.

## Required depdendencies

The following R libraries are required:
- waveslim (tested with version 1.8.2)
- ncdf4 (tested with version 1.17)
- argparser (tested with version 0.4)
- raster (tested with version 2.9-5)

## Directory content

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

The program used to plot Figure 8 of Benestad et al. (2022) (the scientific article is currently under review). The figure shows the precipitation field of Hurricane Katrina on the 28th of August 2005 under different climatological conditions.

```
./plot_fig8_Benestad_et_al_2022.r
```

it generates three figures: "katrina_orig_2005.png" (original precipitation field, panel b) in Fig.8); "katrina_deflated_1961_1990.png" (precipitation field consistent with 1961-1990 climate, panel a) in Fig.8); "katrina_eninflated_2021_2050.png" (precipitation field consistent with an hypothetical 2021-2050 climate, panel c) in Fig.8).



## Copyright and license

Copyright © 2021 Norwegian Meteorological Institute. tpwave is licensed under The GNU Lesser General Public License (LGPL). See LICENSE file.

## Contact

E-mail: Cristian Lussana (cristianl@met.no)
