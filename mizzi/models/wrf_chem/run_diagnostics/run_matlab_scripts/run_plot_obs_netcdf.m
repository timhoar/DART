fname         = 'obs_epoch_001.nc';
region        = [0 360 -90 90 -Inf Inf];
ObsTypeString = 'IASI_O3_RETRIEVAL';
CopyString    = 'NCEP BUFR observation';
QCString      = 'DART quality control';
maxgoodQC     = 4;
verbose       = 1;
twoup         = 1;
plot          = plot_obs_netcdf(fname, ObsTypeString, region, CopyString, ...
                      QCString, maxgoodQC, verbose, twoup);
