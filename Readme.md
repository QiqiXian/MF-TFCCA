# MF-TFCCA

Mixed-frequency time-frequency canonical correlation analysis (MF-TFCCA) is a method for identifying causal relationships between time series of different temporal resolutions.

This repository contains the MATLAB code for the paper: *[Inferring directed spectral information flow between mixed-frequency time series](https://arxiv.org/abs/2408.06109)*

## Requirements

The following MATLAB toolboxes are required to run the code:

- [The Multivariate Granger Causality (MVGC) Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/78727-the-multivariate-granger-causality-mvgc-toolbox) for generating VAR systems and estimating Granger causality.

- MFVAR_toolbox (included in this repository) for estimating Granger causality in mixed-frequency VAR systems.


## Usage

Please refer to `demo_VAR.m` for an example of how to use the code. 

## Citation
For usage of the package and associated manuscript, please cite according to the enclosed citation.bib.
