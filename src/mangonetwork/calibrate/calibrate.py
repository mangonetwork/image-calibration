#!/usr/bin/env python
"""Calibration"""

##########################################################################
#
#   Calibration
#
#   2022-xx-xx  Leslie Lamarche and Asti Bhatt
#               Initial implementation
#
#   2023-03-08  Todd Valentic
#               Make PEP8 compliant
#
#   2026-09-09 Leslie Lamarche
#              Transition to using external asistarcalibration library
#
##########################################################################

import argparse
import configparser
import io
import logging
import os
import sys

from asistarcalibration.starcal import StarCal


def save_calibration_params(cal, output, config_file=None):
    """Save results"""

    config = configparser.ConfigParser()

    # If a real filename is provide for a config file, read it in
    if config_file:
        config.read(config_file)

    config["CALIBRATION_PARAMS"] = dict(
            x0 = str(cal.x0),
            y0 = str(cal.y0),
            rl = str(cal.rl),
            theta = str(cal.theta),
            a = str(cal.A),
            b = str(cal.B),
            c = str(cal.C),
            d = str(cal.D))

    with open(output, "w", encoding="utf-8") as cf:
        config.write(cf)


# -------------------------------------------------------------------------
# Main application
# -------------------------------------------------------------------------


def parse_args():
    """Command line options"""

    parser = argparse.ArgumentParser(description="Calculate camera calibration")

    parser.add_argument("station", help="Station code")
    parser.add_argument("instrument", help="redline or greenline")

    parser.add_argument(
        "-c", "--config", metavar="FILE", help="Existing configuration file"
    )
    parser.add_argument(
        "-s", "--starcal", metavar="FILE", help="Starcal file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="mango-config.ini",
        help="Output configuration filename (default is mango-config.ini)",
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

    if args.config:
        logging.debug("Alternate configuration file: %s", args.config)
        # If configuration file provided, check that it exists
        if not os.path.exists(args.config):
            logging.error("Configurationl file not found")
            sys.exit(1)
        logging.debug("Using provided configuration file: %s", args.config)
        config_file = args.config
    else:
        config_file = None

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


    cal = StarCal(starcal_file)
    cal.calculate_calibration_params(695, 519)
    save_calibration_params(cal, args.output, config_file)



    sys.exit(0)


if __name__ == "__main__":
    main()
