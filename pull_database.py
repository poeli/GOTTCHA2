#!/usr/bin/env python3

#This script allows the user to pull the existing gottcha database
import requests
import sys
import os
import tarfile


GOTTCHA_DB_LATEST = "https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-r90.cg.BacteriaArchaeaViruses.species.tar"

def download_db():
    if os.path.isdir('database'):
        sys.exit('Please make sure a database directory does not exist.')
    os.mkdir('database')
    r = requests.get(GOTTCHA_DB_LATEST)
    open('database/BacteriaArchaeaViruses.species.tar', 'wb').write(r.content)
    tar = tarfile.open('database/BacteriaArchaeaViruses.species.tar','r:gz')
    tar.extractall()
    tar.close()

if __name__ == '__main__':
    download_db()
