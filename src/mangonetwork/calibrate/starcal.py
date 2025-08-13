#!/usr/bin/env python
"""Plot Star Calibration files"""

##########################################################################
#
#   Plot starcal file
#
#   2022-xx-xx  Leslie Lamarche and Asti Bhatt
#               Initial implementation
#
#   2023-03-08  Todd Valentic
#               PEP8 compliance
#
##########################################################################

import datetime as dt
import argparse
import io
import logging
import os
import sys
import requests

import h5py
import numpy as np

# Workaround due to bug in matplotlib event handling interface
# https://github.com/matplotlib/matplotlib/issues/30419
import matplotlib
if matplotlib.get_backend() == 'macosx':
    matplotlib.use('tkagg')
import matplotlib.pyplot as plt

from skyfield.api import Star, load, wgs84
from skyfield.data import hipparcos

if sys.version_info < (3, 9):
    import importlib_resources as resources
else:
    from importlib import resources


class StarCal:
    """Star calibration"""

    def __init__(self, image_file, output_file):

        image = self.load_image(image_file)
        self.find_stars(image)
        self.save_starcal_file(output_file)

    def add_star(self, click):
        """Add user selected star and az,el based on HIP."""
        # This makes use of the Hipparcos Catolog
        # https://rhodesmill.org/skyfield/stars.html

        # Star location in figure from click event
        x = click.xdata
        y = click.ydata
        print(f"Star at {x=:02f}, {y=:02f}")

        # User entered HIP number
        hip = input('HIP #: ')
        print(hip)

        # Look up star based on HIP and calculate az/el
        try:
            s = Star.from_dataframe(self.df.loc[float(hip)])
        except KeyError:
            print(f'Entered Hipparcos designation {hip} is not in database!')
            return

        elev, azmt, _ = self.site_ref.observe(s).apparent().altaz()

        # Append star information
        self.star_azel.append([azmt.degrees, elev.degrees])
        self.star_pos.append([x, y])
        self.star_hip.append(hip)

        # Mark star on plot
        self.ax.scatter(x, y, facecolors='none', edgecolors='r')
        self.fig.canvas.draw()

    def load_image(self, raw_file):
        """Load image and metadata from raw file"""

        image = h5py.File(raw_file, 'r')['image']
        cooked_image = self.prep_image(image)

        self.time = dt.datetime.utcfromtimestamp(image.attrs['start_time'])
        self.site_lat = image.attrs['latitude']
        self.site_lon = image.attrs['longitude']
        self.site_station = image.attrs['station']
        self.site_instrument = image.attrs['instrument']

        return cooked_image

    def prep_image(self, image, contrast=99.95, rotation_angle=0.):
        """Prepare image to display"""

        cooked_image = np.array(image)
        cooked_image = equalize(cooked_image, contrast)
        #cooked_image = imageops.rotate(cooked_image, rotation_angle)

        return cooked_image

    def find_stars(self, image):
        """Display image and track manual selection of stars"""

        #az, el, i, j, raw_file = self.parse_file(star_cal_file)

        self.prep_star_lookup()

        self.star_hip = list()
        self.star_azel = list()
        self.star_pos = list()

        # Display image with stars
        self.fig, self.ax = plt.subplots()
        # Set up button press event trigger
        self.fig.canvas.mpl_connect('button_press_event', self.add_star)
        # Display image
        self.ax.imshow(image, cmap='gray')
        plt.show()


    def prep_star_lookup(self):
        """Prepare skyfield for star lookups"""
        ts = load.timescale()
        t = ts.utc(self.time.year,self.time.month,self.time.day,self.time.hour,self.time.minute,self.time.second)
        planets = load('de421.bsp')
        earth = planets['earth']
        site = earth + wgs84.latlon(self.site_lat, self.site_lon, elevation_m=0)
        self.site_ref = site.at(t)

        with load.open(hipparcos.URL) as f:
            self.df = hipparcos.load_dataframe(f)


    def save_starcal_file(self, output):
        """ Save output starcal file"""

        with open(output, 'w') as f:
            # write header
            f.write(f'# {self.site_station.upper()}    {self.site_instrument}\n')
            f.write(f'# {self.time.isoformat()}\n')
            f.write(f'# GLAT={self.site_lat:10.6f}    GLON={self.site_lon:10.6f}\n')
            f.write(80*'#'+'\n\n')
            f.write(f'# {"HIP":8}{"Azimuth":>20}{"Elevation":>20}{"X":>15}{"Y":>15}\n')

            # add new stars
            for hip, azel, pos in zip(self.star_hip, self.star_azel, self.star_pos):
                f.write(f'{hip:10}{azel[0]:20.10f}{azel[1]:20.10f}{pos[0]:15.5f}{pos[1]:15.5f}\n')

    def parse_file(self, star_cal_file):
        """Read starcal file"""

        raw_filename = star_cal_file.split('\n')[0].split()[-1]

        az, el, i, j = np.loadtxt(io.StringIO(star_cal_file), usecols=(1,2,3,4), unpack=True)

        return az, el, i, j, raw_filename



def equalize(image, contrast, num_bins=10000):
    """Histogram Equalization to adjust contrast [1%-99%]"""
    # copied function from imageops.py
    # needed to make the image visable - there may be more efficient ways of doing this

    image_array_1d = image.flatten()

    image_histogram, bins = np.histogram(image_array_1d, num_bins)
    image_histogram = image_histogram[1:]
    bins = bins[1:]
    cdf = np.cumsum(image_histogram)

    # spliced to cut off non-image area
    # any way to determine this dynamically?  How periminant is it?
    cdf = cdf[:9996]

    max_cdf = max(cdf)
    max_index = np.argmin(abs(cdf - contrast / 100 * max_cdf))
    min_index = np.argmin(abs(cdf - (100 - contrast) / 100 * max_cdf))
    vmax = float(bins[max_index])
    vmin = float(bins[min_index])
    low_value_indices = image_array_1d < vmin
    image_array_1d[low_value_indices] = vmin
    high_value_indices = image_array_1d > vmax
    image_array_1d[high_value_indices] = vmax

    return image_array_1d.reshape(image.shape)





# ------------------------------------------------------------------------
# Main application
# ------------------------------------------------------------------------


def parse_args():
    """Command line parsing"""

    parser = argparse.ArgumentParser(
        description="Check the star identification for calibration"
    )

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-sc", "--starcal", metavar="FILE", help="Alternate starcal file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="starcal-mango.txt",
        help="Output starcal filename (default is starcal-mango.txt)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    return parser.parse_args()


def find_starcal(station, instrument):
    """Find starcal file in package data"""

    starcal_file = f"starcal-{station}-{instrument}.txt"

    logging.debug("Using package starcal file: %s", starcal_file)

    #return resources.files("mangonetwork.raw.data").joinpath(starcal_file).read_text()
    return starcal_file

def download_image(station, instrument):
    """Download image for star matching"""

    dtstr = input('Enter date-time (YYYY-MM-DD HH:MM:SS): ')
    time = dt.datetime.fromisoformat(dtstr)
    print(time)

    url = f'https://data.mangonetwork.org/data/transport/mango/archive/{station.lower()}/{instrument}/raw/{time:%Y}/{time:%j}/{time:%H}/mango-{station.lower()}-{instrument}-{time:%Y%m%d-%H%M%S}.hdf5'
    print(url)
    r=requests.get(url)
    open('mango_image.hdf5', 'wb').write(r.content)
    return 'mango_image.hdf5'


def main():
    """Main application"""

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.starcal:
        logging.debug("Alternate starcal file: %s", args.starcal)
        if not os.path.exists(args.starcal):
            logging.error("Config file not found")
            sys.exit(1)
        with open(args.starcal, encoding="utf-8") as f:
            contents = f.read()
    else:
        contents = find_starcal(args.station, args.instrument)

    image_filename = download_image(args.station, args.instrument)
    #StarCal(contents, args.output)
    StarCal(image_filename, args.output)

    sys.exit(0)


if __name__ == "__main__":
    main()
