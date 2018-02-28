# -*- coding: utf-8 -*-
#!/usr/bin/python3

"""Program to generate sql database for cbs138 strain.
Author: Zhuwei Xu
e-mail: zxu19@jhmi.edu
"""
import sqlite3
import pandas
from lib.biofile import simple_gff3_load, simple_fasta_write, load_fasta
from lib.sequence_lib import rc_seq
import time
import re
import os


def _message(mess=""):
    """Print a message with time.
    Args:
        mess: message to be printed

    """
    time_str = time.strftime("%Y-%m-%d %H:%M:%S")
    print("%s: %s" % (time_str, mess))


def convert_percent_encoding(string: str):
    """Convert the percent encoding characters in the string"""
    encoding = {
        "%20": " ",
        "%21": "!",
        "%23": "#",
        "%24": "$",
        "%26": "&",
        "%27": "'",
        "%28": "(",
        "%29": ")",
        "%2A": "*",
        "%2B": "+",
        "%2C": ",",
        "%2F": "/",
        "%3A": ":",
        "%3B": ";",
        "%3D": "=",
        "%3F": "?",
        "%40": "@",
        "%5B": "[",
        "%5D": "]",
        "'": "''"
    }
    new_str = string
    for code in encoding:
        new_str = re.sub(code, encoding[code], new_str)
    return new_str


def load_cbs138_data(gff: str, cdna: str):
    """Load fasta and gff file of cbs138 genome.
    The published cbs138 genome is loaded, and the raw
    database is built.
    """
    gff_entries, chroms, seqs = simple_gff3_load(gff, return_fasta=True)
    # Remove the previous test version of the database
    if os.path.isfile("cglabrata.db"):
        os.system("rm cglabrata.db")

    database = sqlite3.connect("cglabrata.db")
    _message("Cglabrata Database Generated in cglabrata.db")
    # Load the genome info
    # Load the chromosome information
    chroms = [chrom.split()[0] for chrom in chroms]
    genome = {chroms[idx]: seqs[idx] for idx in range(len(chroms))}
    cdna = load_fasta(cdna)
    cur = database.cursor()
    cur.execute('''CREATE TABLE genome (
        strain  TEXT    NOT NULL,
        source  TEXT    NOT NULL,
        chrom   TEXT    NOT NULL,
        sequence LONGTEXT   NOT NULL,
        note    TEXT
        );''')
    for chrom, seq in genome.items():
        cur.execute('''INSERT INTO genome
            VALUES('CBS138', 'CGD', '%s', '%s', 'From online Candida database')
            ''' % (chrom, seq))
    _message("Online Cbs138 genome loaded into C.glabrata database")
    database.commit()

    # Load the gene information
    cur.execute('''CREATE TABLE gene (
        strain TEXT NOT NULL,
        source TEXT NOT NULL,
        name TEXT NOT NULL,
        gene TEXT,
        chrom TEXT NOT NULL,
        start INT NOT NULL,
        end INT NOT NULL,
        note LONGTEXT,
        direction CHAR(1) NOT NULL,
        sequence LONGTEXT NOT NULL,
        cDNA LONGTEXT,
        utr5 TEXT,
        utr3 TEXT,
        cds_start INT,
        cds_end INT,
        is_multi CHAR(1)
        )''')
    cur.execute('''CREATE TABLE exon(
        strain TEXT NOT NULL,
        source TEXT NOT NULL,
        gene TEXT NOT NULL,
        chrom TEXT NOT NULL,
        direction CHAR(1) NOT NULL,
        no TEXT NOT NULL,
        start INT NOT NULL,
        end INT NOT NULL,
        seq LONGTEXT NOT NULL
        )''')
    gene_exon_list = {}
    for entry in gff_entries:
        if entry[2] == "gene" or entry[2] == "pseudogene":
            chrom = entry[0]
            start = int(entry[3])
            end = int(entry[4])
            seq = genome[chrom][start - 1: end]
            direction = entry[6]
            if direction == "-":
                seq = rc_seq(seq)
            name = entry[8][3:]
            note = ""
            for ent in entry[9:]:
                if ent.startswith('Note='):
                    note = ent[5:]
                    note = convert_percent_encoding(note)
                    break
            if entry[2] == "pseudogene":
                note = "Pseudogene; " + note
            cur.execute('''INSERT INTO gene(
                strain, source, chrom, name, start, end, direction, sequence, note)
                VALUES('CBS138', 'CGD', '%s','%s', %d, %d, '%s', '%s', '%s')
                ''' % (chrom, name, start, end, direction, seq, note))
            for ent in entry[9:]:
                if ent.startswith('Gene='):
                    gene = ent[5:]
                    gene = convert_percent_encoding(gene)
                    cur.execute('''UPDATE gene
                        SET gene = '%s'
                        WHERE strain = 'CBS138' and name = '%s'
                        ''' % (gene, name))
                    break
        if entry[2] == "exon":
            gene_name = entry[8][3: entry[8].index('-T-E')]
            no = int(entry[8][-1])
            start = int(entry[3])
            end = int(entry[4])
            seq = genome[chrom][start - 1: end]
            direction = entry[6]
            if direction == "-":
                seq = rc_seq(seq)
            gene_exon_list.setdefault(gene_name, {})
            gene_exon_list[gene_name][no] = (seq, start, end, direction)

            if direction == "-":
                seq = rc_seq(seq)
            cur.execute('''INSERT INTO exon
                VALUES('CBS138', 'CGD', '%s', '%s', '%s', %d, %d, %d, '%s')
                ''' % (gene_name, chrom, direction, no, start, end, seq))
    _message("Genes and exons loaded into database")

    # Load the cDNA information
    for gene, seq in cdna.items():
        gene_name = gene.split()[0]

        cur.execute('''UPDATE gene SET cdna = '%s'
            WHERE name = '%s' and strain = 'CBS138'
            ''' % (seq, gene_name))
        exons = gene_exon_list[gene_name]
        start_seq = seq[:20].upper()
        end_seq = seq[-20:].upper()
        if len(exons) == 1:
            cur.execute('''UPDATE gene SET is_multi = 'N'
                WHERE name = '%s' and strain = 'CBS138'
                ''' % (gene_name))
        else:
            cur.execute('''UPDATE gene SET is_multi = 'Y'
                WHERE name = '%s' and strain = 'CBS138'
                ''' % (gene_name))
        # Locate the 5'UTR
        exon_seqs = [exons[idx + 1][0] for idx in range(len(exons))]
        mrna = "".join(exon_seqs)
        mrna = mrna.upper()
        try:
            start_id = mrna.index(start_seq)
        except ValueError:
            _message("Error Working on utr for gene %s" % gene_name)
            print(mrna)
            print(start_seq)
            raise ValueError
        end_id = mrna.rindex(end_seq)
        cur_left = 0
        pos_utr5 = 0
        pos_utr3 = 0
        utr5 = ""
        utr3 = ""
        for start_no in range(1, len(exons) + 1):
            if start_id <= cur_left + exons[start_no][2] - exons[start_no][1] + 1:
                break
            cur_left += exons[start_no][2] - exons[start_no][1] + 1
        cur_right = 0
        for end_no in range(1, len(exons) + 1):
            if end_id <= cur_left + exons[end_no][2] - exons[end_no][1] + 1:
                break
            cur_right += exons[end_no][2] - exons[end_no][1] + 1

        rem_pos = start_id
        for no in range(1, start_no):
            rem_pos = rem_pos - (exons[no][2] - exons[no][1] + 1)
            utr5 += exons[no][0]
        utr5 += exons[start_no][0][: rem_pos]
        if exons[1][-1] == "+":
            pos_utr5 = exons[start_no][1] + rem_pos
        else:
            pos_utr5 = exons[start_no][2] - rem_pos
        rem_pos = len(mrna) - end_id - 1
        for no in range(end_no + 1, len(exons) + 1):
            rem_pos = rem_pos - (exons[no][2] - exons[no][1] + 1)
            utr3 += exons[no][0]
        if rem_pos != 19:
            utr3 = exons[end_no][0][-rem_pos + 19:] + utr3
        if exons[1][-1] == "+":
            pos_utr3 = exons[end_no][2] - rem_pos + 19
        else:
            pos_utr3 = exons[end_no][1] + rem_pos - 19
        cur.execute('''UPDATE gene
            SET
            utr5 = '%s',
            utr3 = '%s',
            cds_start = %d,
            cds_end = %d
            WHERE name = '%s' and strain = 'CBS138'
            ''' % (utr5, utr3, pos_utr5, pos_utr3, gene_name))
    database.commit()

    cur.close()
    database.close()


def replace_fosmid(fasta: str, gff: str):
    """Replace the data of cbs138 online genome with fosmid sequence.
    The publisehd cbs138 genome is known to have problems.  The subtelomeric
    regions are cloned into fosmids, and re-sequenced by Sanger sequencing.
    This program updated the fosmid regions.
    """
    pass
