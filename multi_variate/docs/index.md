# Home

## About
mv-clasp is a tool for identifying mutated segments across multiple allele frequencies.

## Installation
In a conda environment:
```
pip install mv-clasp
```
!!! tip
    For support in configuring a conda environment see the [conda documentation website].

## Usage
### Detecting change points:
```
from mv_clasp import MultivariateClaSP

mc = MultivariateClaSP(./allele_frequencies_data.filetype, mode, ./results_directory/)

mc.analyze_time_series()
```
!!! note
    While the frequencies argument on initialization defaults to find columns in the input file named "vaf", "median_baf", and "median_dr",
    you should ensure that the list of string arguments supplied to frequencies match column names within your specific data file. Failure to do so will result in an error message.

### Plotting results:
```
mc.plot_combined_profile(title="sum", save=True)
```
![Output](../images/plot_profile.png)