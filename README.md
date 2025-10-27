# calibration
This repository contains the scripts required to generate calibration files for MANGO imagers.  This is a somewhat labor intenseive process that requires manually identifying stars.

This includes a program for generating the starcal files (lists of stars and their relative positions in an image) as well as a program that performs fitting from these files to find the rotation angle and lens funtion that will be applied in raw image processing.

## Installation
This package can be installed from GitHub with pip.  Either clone and install
```
git clone https://github.com/mangonetwork/image-calibration.git
cd image-calibration
pip install .
```
or install directly from GitHub.
```
pip install git+https://github.com/mangonetwork/image-calibration.git
```

## Command Line Programs
This package primarily functions through command line programs, described below.  Before any of these functions will work, you MUST set the environment variable `MANGONETWORK_CONFIGS` to the path that is the default location for all your configuration and starcal files.
```
export MANGONETWORK_CONFIGS=/path/to/default/config/dir
```
You will have to run this export command every time you open a new terminal session.

### mango-starcal
This program generates starcal files (lists of stars and their relative positions in an image).  Running it will open an image window.  Users should click on stars in the windown and then enter the corresponding HIP at the prompt.  A good calibration will require 10-20 stars identified.

Basic usage
```
mango-starcal <station> <instrument>
```
This will load the default starcal file for this station and instrument and allow you to add additional stars.  It will raise an error if there is no default starcal file for this station and instrument.

Create a new file from scratch
```
mango=starcal <station> <instrument> -n -t <YYYY-MM-DDTHH:MM:SS>
```
This will create a new starcal file for this station and instrument using the time (ISO format) specifed in the `-t` flag.

Specify output file
The `-o` flag lets you explicitly name the output file.  The default is `mango-starcal.txt`.  Note that the output will always be saved to a new file, regardless of if the `-n` flag is used or not.  This program does not append to the existing or provided startcal files, however the stars listed in the starting starcal fill will be included in the new output file.

Specify a starting starcal file
The `-s` flag lets to explicity specify an input starcal file instead of using the default.  Note that this file sill will not be overwritten and the output will be saved as described above.

### mango-calibrate
This produces the calibration prameters required in the config file by finding the best fit of these parameters for the stars in the starcal file.  It will open a figure window so you can verify the fit is reasonable.

Basic usage
```
mango-calibrate <station> <instrument>
```
This will perform a calibration with the default starcal file.

Specify output file
The `-o` flag lets you explicitly name the output file.  The default is `mango-config.ini`.  Note that the output will always be saved to a new file.

Specify the starcal file
The `-s` flag lets you specify a starcal file to use other than the default.

Specify the configuration file
The `-c` flag lets you specify an existing configuration file.  This file will not be modified, but its contents (with the `CALIBRATATION` section updated) will be copied to the output file.

### mango-checkcal
Check that the calculated calibration makes sense with the original image used to generate the starcal file.

Basic usage
```
mango-checkcal <station> <instrument>
```
This will open the original image with elevation angle and north marked.

Specify the starcal file
The `-s` flag lets you specify a custom starcal file.

Specify the configuration file
The `-c` flag lets you specify a custom configuration file.


## Step-by-Step Guide to Calibrating a New MANGO Camera
Before you begin, make sure you have installed both this package and [Stellarium](https://stellarium.org).

1. Find a clear image from the camera in question.  This can usually be done by flipping through the quicklook movies on the [MANGONETWORK](https://www.mangonetwork.org/mango/v1/database/sites) website.  You ideally want to find a time where you can see the stars in the image clearly without clouds or haze and with moon down.

2. Run `mango-starcal` for this camera and time.
```
mango-starcal <site> <instrument> -n -t <YYYY-MM-DDTHH:MM:SS> -o starcal-<site>-<instrument>.txt
```
This will open a matplotlib window with the raw image you selected.  DO NOT close this window.

3. Open Stellarium and set to the date and time of the image.  The terminal window running the `mango-starcal` command should have the time and coordinates printed.
    a. In the Configuration window (F2), under the "Extras" tab in "Additional information and settings" check "Use decimal degrees".
    b. Pause time progression by hitting the "Play/Pause" button in the lower pop up bar or by pressing "K".
    c. In the Location window (F6), set the latitude and longitude of the site.  Check "Use custom time zone" and set the time zone in the drop down menu to "UTC+00:00". If you want to make it easier to pull up this site in the future, enter the site name/abbreviation in the "Name/City" field and click "Add to list".
    d. In the Date/time window (F5), set the date and time.
    e. In the Sky and viewing options window (F4),  under the "Sky" tab in "Projection" select "Fish-eye". Also under the "SSO" tab uncheck "Solar System objects".
    f. In the lower pop up bar, turn off all constellation highlighting, any grids, deep sky objects, planet labels, exoplanets, meteroid showers, and artificial satellites.  It is usually helpful to have "Ground" on.  "Atmosphere" may or may not be helpful.
    g. Rotate the view and zoom in or out so you can see the entire circle of the sky on your screen.
4. Find some identifiable stars or constellations in the MANGO image and rotate the Stellarium view so the two are approximately oriented in the same direction.
5. Select a star in Stellarium and click on the same star in the matplotlib window (it will look like nothing happens, but only click ONCE). The terminal window should now additionally list the the x, y coordinates of the star you just selected and have a promp for the "HIP #". Enter the HIP number that Stellarium lists in the header of information about that star and hit enter. The matplotlib window should now have a red circle around that star.
6. Repete step 5 until you've selected 10-20 stars.  Aim to get a good azimuthal spread of stars across the entire image.
7. When done selecting stars, close the matplotlib window.  The file listed as the output file in step 2 should have been created.  If you open this file, it contains a header with camera and time information and a table of every star you selected with its name (if available), HIP number, aximuth, elevation, x coordinate, and y coordiinate.
8. Run `mango-calibrate` with this starcal file.
```
mango-calibrate <site> <instrument> -s starcal-<site>-<instrument>.txt -o <site>-<instrument>.ini
```
This should pop up a window showing how the stars match the fitted rotation and lens function.  As long as nothing looks wierd, you can close this window and the output file will be created.  This is an \*.ini file that only includes the `CALIBRATION_PARAMS` section.
9. Run `mango-checkcal` to confirm the calibration parameters transforms make sense on the original image.
10. Move the newly createcd starcal and config files to your default config file directory.  You should add the `SITE_INFO`, `PROCESSING`, and `QUICKLOOK` fields to the config file so it can be used for processing data.  Alternatively, if you have a prior config file with these sections, you can specify it with the `-c` flag in step 8 and they will automatically be copied.

**Notes:**
- This program ONLY works with stars, not planets, moons or other objects.  Based on the HIP catelog.
- Stellarium appends the HIP number of some stars with alphebetical characters, usually indicating a binary star system.  In these cases, just enter the numerical digits.
- If the figure that `mango-calibrate` pops up has a star that is wildly out of place, manually edit the starcal file and delete that line.  These can offset the fit and result in bad calibration.
- You can create a new starcal file starting with the image and stars in an existing starcal file by specifying the existing file with the `-s` flag.
```
mango-starcal <site> <instrument> -s <existing_starcal_file>.txt -o starcal-<site>-<instrument>.txt
```
This is useful if you saved a starcal file midway through working on it, or realize you need to add more stars to improve the fit.
- After config and starcal files are added to the default location, it is not nessisary to explicity list them with the `-c` and `-s` flages to use them.
