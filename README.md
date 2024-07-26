### DOWNLOAD THE DATA: `download_wget_bash.csh`

csh script to download data from FRF website; uncomment and adjust what you need

### EXTRACT WAVE DATA: `nc_to_dat.py`

Python script to extract the wave data from the NetCDF files of the FRF and generate a text file with the date, hs, tp, dir and level. 

It can be adapted to resample the hourly original data into daily data, with or without interpolation.

Change the name of the final file at the end of the script. Here it generated `forcing_duck_daily_8m_1990-2022.dat`.

### ANALYSE WAVE DATA: `waves_statistics.py`

Python script that analyse the timeseries and plot the following figures:
- Figure 1: timeseries of Hs and Tp
- Figure 2: windrose of wave direction
- Figure 3: density, probability of exceedance and return periods for Hs and Tp
- Figure 4: joint repartition of (Hs, Tp) and levels of density
- Figure 5: omegas
- Figure 6: comparison between omegas(Hs) and omegas(Hb)
- Figure 7: omegas repartition by season
- Figure 8: omegas repartition by consecutive time
- Figure 9: omegas timeseries and monthly means of omegas

Run the script and follow the questions asked, answering by yes (y) or no (n).

The line 519 allows to save the analysis, here in `stats_wavesduck_daily_8m_1990-2022.csv`.  

### ANALYSE, DETECT, SHORELINE/BARS: `morpho_evolution.py`

Python script to analyse the bathymetries, extract the shoreline position, detect bars and extract their position, height, width. Works with `duck_configuration.yaml` and functions in `libs` directory.

The following list explains where to modify the script and why:

- search 'savefig' and uncomment the lines if you want to save the plots;

- in the function 'bars\_dimensions', comment/uncomment the sections `KEEP ALL BARS BETWEEN 0 AND 180 M` and `KEEP BAR CORRESPONDING TO THORNTON AND GUZA DEPTH OF BREAKING` to use one, both or none. The value of 180 can be adapted and modified lines 755 and 782;

- still in this function, look for 'math.isclose' and adapt the 'abs\_tol' values;

- for the plots, look for 'set\_xlim' and 'set\_ylim' and adapt the dimensions to your study;

- if you use the function 'omega\_eq', change the value of phi, line 1020;

- in the function 'preprocessing', adapt the values of Xorig, Yorig, Lx andLy to your study; the rotations lines 1149, 1167 and 1176 may be also adapted;

- change the name of your configuration file line 1193;

- change the name of your wave file line 1212;

- change the date or comment the line 1216;

- adapt the directory for your bathymetries line 1243;

- if you want to save the data, comment/uncomment the corresponding lines at the end of the script;

- the script generated `evolution_bars_frag_thornton_2000-2022.csv`, `evolution_bars_height.csv` and `evolution_bars_width.csv`.
