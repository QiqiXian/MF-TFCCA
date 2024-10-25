# MF-TFCCA

Mixed-frequency time-frequency canonical correlation analysis (MF-TFCCA) is a method for identifying causal relationships between time series of different temporal resolutions.

This repository contains the MATLAB code for the paper: *[Inferring directed spectral information flow between mixed-frequency time series](https://arxiv.org/abs/2408.06109)*

## Repo Contents
- `src`: MATLAB code for MF-TFCCA and related functions.
- `demo`: MATLAB scripts for running a simulation demo.
- `Data`: Data samples。

## Requirements

### Hardware and OS

The `MF-TFCCA` package requires only a standard computer with enough RAM to support the in-memory operations. The runtimes are tested on a computer with the following specs: RAM: 16GB; CPU: 6 cores@2.38 GHz

The package is tested on Windows 10 operating system. The package should work on other operating systems as well.

### Software

Before setting up the package, please ensure that you have `MATLAB` installed on your system. The package is tested on `MATLAB R2021a`.

The following MATLAB toolboxes are required to run the code:

- [The Multivariate Granger Causality (MVGC) Toolbox](https://www.mathworks.com/matlabcentral/fileexchange/78727-the-multivariate-granger-causality-mvgc-toolbox) for generating VAR systems and estimating Granger causality. (version mvgc_v1.3)

This package can be installed using the given url and will install in within 30 seconds.

- MFVAR_toolbox (included in this repository) for estimating Granger causality in mixed-frequency VAR systems.

## Installation Guide

The package can be used directly by installing the required toolboxes and adding the `src` folder and its subfolders to the MATLAB path.

## Usage

Please refer to `demo/demo_VAR.m` for an example of how to use the code. Each block of the code should take less than 1 minute to run on a standard computer using the given data samples.

## Citation
For usage of the package and associated manuscript, please cite according to the enclosed [citation.bib](https://github.com/QiqiXian/MF-TFCCA/blob/main/citation.bib).
