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
#   2026-09-09 Leslie Lamarche
#              Transition to using external asistarcalibration library
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

from asistarcalibration.starfinder import StarFinder
from asistarcalibration.wizard import equalize


def load_image(raw_file):
    """Load image and metadata from raw file"""

    image = h5py.File(raw_file, 'r')['image']
    cooked_image = prep_image(image)

    time = dt.datetime.utcfromtimestamp(image.attrs['start_time'])
    site_lat = image.attrs['latitude']
    site_lon = image.attrs['longitude']
    #site_station = image.attrs['station']
    #site_instrument = image.attrs['instrument']

    return cooked_image, time, site_lat, site_lon


def prep_image(image, contrast=99., rotation_angle=0.):
    """Prepare image to display"""

    cooked_image = np.array(image)
    cooked_image = equalize(cooked_image, contrast)

    return cooked_image





# ------------------------------------------------------------------------
# Main application
# ------------------------------------------------------------------------


def parse_args():
    """Command line parsing"""

    parser = argparse.ArgumentParser(
        description="Manually identify stars for calibration"
    )

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")
    parser.add_argument("-n", "--new", action="store_true", default=False, help="Generate new file from scratch")
    parser.add_argument("-t", "--time", help="Time for star idenfication")

    parser.add_argument(
        "-s", "--starcal", metavar="FILE", help="Existing starcal file (for appending stars)"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="mango-starcal.txt",
        help="Output starcal filename (default is mango-starcal.txt)",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")

    return parser.parse_args()


def find_starcal(station, instrument):
    """Find starcal file in package data"""

    # Placeholder for default config file location
    #   This function can be rewritten later
    config_dir = os.environ['MANGONETWORK_CONFIGS']

    starcal_file = os.path.join(config_dir, f"starcal-{station}-{instrument}.txt")

    logging.debug("Using package starcal file: %s", starcal_file)

    #return resources.files("mangonetwork.raw.data").joinpath(starcal_file).read_text()
    return starcal_file


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


def main():
    """Main application"""

    args = parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    if args.new:
        # If new flag set, generate a fresh starcal file
        logging.debug("Generating new starcal file")
        station = args.station
        instrument = args.instrument
        time = dt.datetime.fromisoformat(args.time)
        starcal_file = None

    elif args.starcal:
        # If starcal file provided, read in header
        if not os.path.exists(args.starcal):
            logging.error("Starcal file not found")
            sys.exit(1)
        logging.debug("Using provided starcal file: %s", args.starcal)
        starcal_file = args.starcal
        station, instrument, time = read_header(starcal_file)

    else:
        # If no starcal file provided, find the default and read in header
        starcal_file = find_starcal(args.station, args.instrument)
        if not os.path.exists(starcal_file):
            logging.error("No default starcal file found for %s %s!", args.station, args.instrument)
            sys.exit(1)
        logging.debug("Using defalt starcal file: %s", starcal_file)
        station, instrument, time = read_header(starcal_file)

    # Download image
    image_filename = download_image(station, instrument, time)

    # Load image and retrieve actual time and site coordinates
    img, truetime, site_lat, site_lon = load_image(image_filename)

    # Run star calibration
    find = StarFinder(site_lat, site_lon, truetime, station=station, instrument=instrument)
    find.load_stars(starcal_file)
    find.find_stars(img)
    find.save_starcal_file(args.output)


    sys.exit(0)


if __name__ == "__main__":
    main()
