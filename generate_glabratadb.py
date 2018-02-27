# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""Program to generate sql database for cbs138 strain.
Author: Zhuwei Xu
e-mail: zxu19@jhmi.edu
"""
import sqlite3
import pandas
from biofile import simple_gff3_load, simple_fasta_write, load_fasta
import time
import re


def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print("%s: %s" % (time_str, mess))


def load_cbs138_data(fasta: str, gff: str):
    """Load fasta and gff file of cbs138 genome.
    The published cbs138 genome is loaded, and the raw
    database is built.
    """
    genome = load_fasta(fasta)
    gff_entries = simple_gff3_load(gff)

    # Load the genome into


def replace_fosmid(fasta: str, gff: str):
    """Replace the data of cbs138 online genome with fosmid sequence.
    The publisehd cbs138 genome is known to have problems.  The subtelomeric
    regions are cloned into fosmids, and re-sequenced by Sanger sequencing.
    This program updated the fosmid regions.
    """
    pass


def main():
    """Generate cbs138 database with online database and fosmid information."""
    pass


if __name__ == '__main__':
    main()
