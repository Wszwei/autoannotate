# -*- coding: utf-8 -*-
#!/usr/bin/python3
import os
from generate_glabratadb import load_cbs138_data


def main():
    """Generate cbs138 database with online database and fosmid information."""
    # Load the online database
    dirn = "/home/zhuwei/data/processdata/anno/180226_auto_annotation/"
    os.chdir(dirn)
    load_cbs138_data(
        gff="cbs138.seq.gff",
        cdna='cbs138.online.cdna.fa')


if __name__ == '__main__':
    main()
