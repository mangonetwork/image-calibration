# calibration
This repository contains the scripts required to generate calibration files for MANGO imagers.

This includes a program for generating the starcal files (lists of stars and their relative positions in an image) as well as a program that performs fitting from these files to find the rotation angle and lens funtion that will be applied in raw image processing.

## mango-starcal
This program generates starcal files (lists of stars and their relative positions in an image).  Running it will open an image window.  Users should click on stars in the windown and then enter the corresponding HIP at the prompt.  A good calibration will require 10-20 stars identified.
