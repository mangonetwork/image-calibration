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

from asistarcalibration.starcal import StarCal
from asistarcalibration.wizard import equalize


def run_checkcal(starcal_file, config_file):

    # Load config file
    config = configparser.ConfigParser()
    config.read(config_file)

    # Read header information from starcal file
    station, instrument, time = read_header(starcal_file)

    # prepare image
    image_file = download_image(station, instrument, time)
    image, site_lat = load_image(image_file)
    cooked_image = prep_image(image)


    sc = StarCal(starcal_file)

    load_calibration_params(sc, config)

    sc.checkcal(image, site_lat)


def read_header(starcal_file):
    """Read header from starcal file"""

    with open(starcal_file, 'r') as f:
        line1 = f.readline()
        _, station, instrument = line1.split()
        line2 = f.readline()
        time = dt.datetime.fromisoformat(line2.split()[1])

    return station, instrument, time

def download_image(station, instrument, time):
    """Download image for star matching"""

    url = f'https://data.mangonetwork.org/data/transport/mango/archive/{station.lower()}/{instrument}/raw/{time:%Y}/{time:%j}/{time:%H}/mango-{station.lower()}-{instrument}-{time:%Y%m%d-%H%M%S}.hdf5'
    logging.debug("Downloading raw image file: %s", url)
    r=requests.get(url)
    open('mango_image.hdf5', 'wb').write(r.content)

    return 'mango_image.hdf5'

def load_image(raw_file):
    """Load image and metadata from raw file"""

    image = h5py.File(raw_file, 'r')['image']
    cooked_image = prep_image(image)

    time = dt.datetime.utcfromtimestamp(image.attrs['start_time'])
    site_lat = image.attrs['latitude']
    site_lon = image.attrs['longitude']
    site_station = image.attrs['station']
    site_instrument = image.attrs['instrument']

    return cooked_image, site_lat

def prep_image(image, contrast=99, rotation_angle=0.):
    """Prepare image to display"""

    cooked_image = np.array(image)
    cooked_image = equalize(cooked_image, contrast)

    return cooked_image

def load_calibration_params(sc, config):

    sc.x0 = config.getfloat("CALIBRATION_PARAMS", "X0")
    sc.y0 = config.getfloat("CALIBRATION_PARAMS", "Y0")
    sc.rl = config.getfloat("CALIBRATION_PARAMS", "RL")
    sc.theta = config.getfloat("CALIBRATION_PARAMS", "THETA")

    sc.A = config.getfloat("CALIBRATION_PARAMS", "A")
    sc.B = config.getfloat("CALIBRATION_PARAMS", "B")
    sc.C = config.getfloat("CALIBRATION_PARAMS", "C")
    sc.D = config.getfloat("CALIBRATION_PARAMS", "D")



####################################################################################


def parse_args():
    """Command line options"""

    parser = argparse.ArgumentParser(description="Check calibration against original starcal image")

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-c", "--config", metavar="FILE", help="Alternate configuration file"
    )
    parser.add_argument(
        "-s", "--starcal", metavar="FILE", help="Alternate starcal file"
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

    # Determine config filename
    if args.config:
        logging.debug("Alternate configuration file: %s", args.config)
        # If configuration file provided, check that it exists
        if not os.path.exists(args.config):
            logging.error("Configurationl file not found")
            sys.exit(1)
        logging.debug("Using provided configuration file: %s", args.config)
        config_file = args.config
    else:
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


    run_checkcal(starcal_file, config_file)

    sys.exit(0)


if __name__ == "__main__":
    main()
