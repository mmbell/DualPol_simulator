# DualPol Simulator

## A dual-polarization radar simulator

This program simulates the dual-polarization radar returns from weather targets using a single, double, or fully explicit bin dropsize distribution. It is currently configured to work with WRF output and liquid water only. The simulator has been used in:

Brown, B. R., M. M. Bell, and A. J. Frambach, 2016: "Validation of Simulated Hurricane Drop Size Distributions using Polarimetric Radar", *Geophys. Res. Lett.*, **42**, doi:10.1002/2015GL067278.

Inputs:
 + Type of DSD (single, double, or full)
 + Name of WRF output file

Outputs:
 + NetCDF file with polarimetric radar variables

Current output is radar reflectivity factor at horizontal and vertical polarization, and differential reflectivity for raindrops only. Additional radar variables and hydrometeor species will be added in due course.

## Copyright

Copyright (c) 2016 Michael M. Bell

