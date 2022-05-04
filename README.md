# ISBAS: Intermittent Small Baseline Subset InSAR Timeseries algorithm

This code implements the Intermittent SBAS method described by [Sowter et al. (2013)](https://dx.doi.org/10.1080/2150704X.2013.823673) for use with GMTSAR-generated interferograms. The algorithm implements the standard SBAS timeseries approach, but allows any pixel to be included so long as it is coherent in at least a set number of interferograms. This results in a significantly improved number of coherent pixels in the final timeseries, with the tradeoff being higher noise levels for those pixels that are only coherent sometimes.

Basic usage: 

1. Run `intf_pad.sh` - this makes sure all images have the same exterior bounds by padding them with NaN values.
2. Set your parameters in the config file - see the example
3. Run `python isbas.py isbas.config`.
