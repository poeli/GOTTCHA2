#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile
from tqdm import *
import argparse
import pkg_resources
import json
from gottcha import GOTTCHA_DB_LATEST, FILE_NAME



def parse_params(args):
    parser = argparse.ArgumentParser(prog='pull_database.py', description="""This script will pull the latest version of the Gottcha2 database.""")

    parser.add_argument('-u', '--url', default=GOTTCHA_DB_LATEST,
                    help='specify a URL to pull from (will override the default)')
    parser.add_argument('-r', '--rank', default='species',
                    help='taxonomic rank of the database (superkingdom, phylum, class, order, famiily, genus, species)')
    return parser.parse_args(args)

def download_db(argvs):
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')

    os.mkdir('database')

    if argvs.url:
        GOTTCHA_DB_LATEST = argvs.url
    with requests.get(GOTTCHA_DB_LATEST.replace('species',argvs.rank,1), stream=True) as r:
        r.raise_for_status()
        with open(FILE_NAME, 'wb') as f:
            pbar = tqdm(total=int(r.headers['Content-Length']),unit='B', unit_scale=True, unit_divisor=1024)
            for chunk in r.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
                    pbar.update(len(chunk))
    tar = tarfile.open(FILE_NAME.replace('species',argvs.rank,1))
    tar.extractall('database')
    tar.close()


def main(args):
    argvs = parse_params(args)
    download_db(argvs)

if __name__ == '__main__':
    main(sys.argv[1:])
