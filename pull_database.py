#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile
from tqdm import *
import argparse as ap


GOTTCHA_DB_LATEST = "https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.tar"
FILE_NAME = 'database/BacteriaArchaeaViruses.species.tar'

def parse_params(args):
    parser = argparse.ArgumentParser(add_help=False)

    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='gottcha2 pull will automatically download the latest database for Gottcha2 taxonomic profiler')

    return p.parse_args([args])

def download_db():
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')

    os.mkdir('database')
    with requests.get(GOTTCHA_DB_LATEST, stream=True) as r:
        r.raise_for_status()
        with open(FILE_NAME, 'wb') as f:
            pbar = tqdm(total=int(r.headers['Content-Length']),unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    tar = tarfile.open('database/BacteriaArchaeaViruses.species.tar','r:gz')
    tar.extractall()
    tar.close()


def main(args):
    if args != None:
        argvs = parse_params(args)
        if argvs.help:
            sys.exit()
    else:
        download_db()

if __name__ == '__main__':
    main(sys.argv[1:])
