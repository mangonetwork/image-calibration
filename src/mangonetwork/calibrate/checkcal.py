# checkcal.py

import argparse
import configparser
import io
import logging
import os
import sys
import datetime as dt
import requests
import h5py
import numpy as np
import matplotlib.pyplot as plt


class Check:
    """check calibration on starcal image"""

    def __init__(self, starcal_file, config_file):

        self.starcal_file = starcal_file
        self.config_file = config_file

        self.config = configparser.ConfigParser()
        self.config.read(config_file)

        self.station, self.instrument, self.time = self.read_header(starcal_file)

        self.run()

    def run(self):

        image_file = self.download_image(self.station, self.instrument, self.time)
        image = self.load_image(image_file)
        cooked_image = self.prep_image(image)
        self.load_calibration_params()
        self.display(cooked_image)

    def read_header(self, starcal_file):
        """Read header from starcal file"""
    
        with open(starcal_file, 'r') as f:
            line1 = f.readline()
            _, station, instrument = line1.split()
            line2 = f.readline()
            time = dt.datetime.fromisoformat(line2.split()[1])
    
        return station, instrument, time

    def download_image(self, station, instrument, time):
        """Download image for star matching"""
    
        url = f'https://data.mangonetwork.org/data/transport/mango/archive/{station.lower()}/{instrument}/raw/{time:%Y}/{time:%j}/{time:%H}/mango-{station.lower()}-{instrument}-{time:%Y%m%d-%H%M%S}.hdf5'
        logging.debug("Downloading raw image file: %s", url)
        r=requests.get(url)
        open('mango_image.hdf5', 'wb').write(r.content)
    
        return 'mango_image.hdf5'

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

        return cooked_image

    def load_calibration_params(self):

        self.x0 = self.config.getfloat("CALIBRATION_PARAMS", "X0")
        self.y0 = self.config.getfloat("CALIBRATION_PARAMS", "Y0")
        self.rl = self.config.getfloat("CALIBRATION_PARAMS", "RL")
        self.theta = self.config.getfloat("CALIBRATION_PARAMS", "THETA")

        self.A = self.config.getfloat("CALIBRATION_PARAMS", "A")
        self.B = self.config.getfloat("CALIBRATION_PARAMS", "B")
        self.C = self.config.getfloat("CALIBRATION_PARAMS", "C")
        self.D = self.config.getfloat("CALIBRATION_PARAMS", "D")


    def display(self, image):
        # Display image with stars
        self.fig, self.ax = plt.subplots()
        # Display image
        self.ax.imshow(image, cmap='gray')

        self.ax.scatter(self.x0, self.y0, edgecolor='gold', facecolor='none', label='zenith')
        t = np.linspace(0., 2*np.pi, 100)
        #self.ax.scatter(self.rl*np.cos(t)+self.x0, self.rl*np.sin(t)+self.y0, edgecolor='gold', facecolor='none', label='zenith')

        r = np.linspace(0., 1., 100)
        lam = np.rad2deg(self.A + self.B * r + self.C * r**2 + self.D * r**3)
        for el in [0., 15., 30., 60.]:
            r0 = np.interp(el, lam[::-1], r[::-1])
            self.ax.plot(r0*self.rl*np.cos(t)+self.x0, r0*self.rl*np.sin(t)+self.y0, label=f'el={el}')

        self.ax.plot([self.x0, self.x0+self.rl*np.sin(np.deg2rad(self.theta))], [self.y0, self.y0+self.rl*np.cos(np.deg2rad(self.theta))], color='k', linestyle=':', label='N')

        r0 = np.interp(self.site_lat, lam[::-1], r[::-1])
        r0 = r0*self.rl
        self.ax.scatter(self.x0+r0*np.sin(np.deg2rad(self.theta)), self.y0+r0*np.cos(np.deg2rad(self.theta)), edgecolor='magenta', facecolor='none', label='Polaris')

        self.ax.legend()

        plt.show()



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


####################################################################################


def parse_args():
    """Command line options"""

    parser = argparse.ArgumentParser(description="Calculate camera calibration")

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-c", "--config", metavar="FILE", help="Alternate configuration file"
    )
    parser.add_argument(
        "-sc", "--starcal", metavar="FILE", help="Alternate starcal file"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    return parser.parse_args()


def find_config(station, instrument):
    """Find configuration file from pacakge data"""

    # Placeholder for default config file location
    #   This function can be rewritten later
    config_dir = os.environ['MANGONETWORK_CONFIGS']

    config_file = os.path.join(config_dir, f"{station}-{instrument}.ini")

    logging.debug("Using package configuration file: %s", config_file)

    #return resources.files('mangonetwork.raw.data').joinpath(config_file).read_text()
    return config_file


def find_starcal(station, instrument):
    """Find starcal file in package data"""

    # Placeholder for default config file location
    #   This function can be rewritten later
    config_dir = os.environ['MANGONETWORK_CONFIGS']

    starcal_file = os.path.join(config_dir, f"starcal-{station}-{instrument}.txt")

    logging.debug("Using package starcal file: %s", starcal_file)

    #return resources.files("mangonetwork.raw.data").joinpath(starcal_file).read_text()
    return starcal_file


def main():
    """Main application"""

    args = parse_args()

    fmt = "[%(asctime)s] %(levelname)s %(message)s"

    if args.verbose:
        logging.basicConfig(format=fmt, level=logging.DEBUG)
    else:
        logging.basicConfig(format=fmt, level=logging.INFO)

#    if args.config:
#        logging.debug("Alternate configuration file: %s", args.config)
#        if not os.path.exists(args.config):
#            logging.error("Config file not found")
#            sys.exit(1)
#        with open(args.config, encoding="utf-8") as f:
#            config_contents = f.read()
#    else:
#        config_contents = find_config(args.station, args.instrument)

#    if args.starcal:
#        logging.debug("Alternate starcal file: %s", args.starcal)
#        if not os.path.exists(args.starcal):
#            logging.error("StarCal file not found")
#            sys.exit(1)
#        with open(args.starcal, encoding="utf-8") as f:
#            starcal_contents = f.read()
#    else:
#        starcal_contents = find_starcal(args.station, args.instrument)
#
#
#    if args.new:
#        # If new flag set, generate a fresh starcal file
#        logging.debug("Generating new starcal file")
#        station = args.station
#        instrument = args.instrument
#        time = dt.datetime.fromisoformat(args.time)
#        starcal_file = None

#    # Determine configuration filename
#    if args.new:
#        # If new flag set, generate a fresh starcal file
#        logging.debug("Generating new config file")
#        config_file = None
#    elif args.config:
#        logging.debug("Alternate configuration file: %s", args.config)
#        # If configuration file provided, check that it exists
#        if not os.path.exists(args.config):
#            logging.error("Configurationl file not found")
#            sys.exit(1)
#        logging.debug("Using provided configuration file: %s", args.config)
#        config_file = args.config
#    else:
    # If no configuration file provided, find the default
    config_file = find_config(args.station, args.instrument)
    if not os.path.exists(config_file):
        logging.error("No default configuration file found for %s %s!", args.station, args.instrument)
        sys.exit(1)
    logging.debug("Using defalt configuration file: %s", config_file)


    # Determine starcal filename
    if args.starcal:
        # If starcal file specified, check that it exists
        logging.debug("Alternate starcal file: %s", args.starcal)
        if not os.path.exists(args.starcal):
            logging.error("Starcal file not found")
            sys.exit(1)
        logging.debug("Using provided starcal file: %s", args.starcal)
        starcal_file = args.starcal
    else:
        # If no starcal file provided, find the default
        starcal_file = find_starcal(args.station, args.instrument)
        if not os.path.exists(starcal_file):
            logging.error("No default starcal file found for %s %s!", args.station, args.instrument)
            sys.exit(1)
        logging.debug("Using defalt starcal file: %s", starcal_file)


    Check(starcal_file, config_file)

    sys.exit(0)


if __name__ == "__main__":
    main()
