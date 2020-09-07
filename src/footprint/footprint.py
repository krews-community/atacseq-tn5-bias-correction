#!/usr/bin/env python

import sys
import argparse

from typing import List, Tuple
from pysam import Fastafile, Samfile

from rgt.Util import ErrorHandler, HmmData, GenomeData, OverlapType
from rgt.HINT.signalProcessing import GenomicSignal
from rgt.HINT.biasTable import BiasTable

from .constants import *

def expandRegion(chromosome, start, end, name = None, w = 500):
    m = int((int(start) + int(end)) / 2)
    return chromosome, m - w, m + w, name

def regionDict(k, forward, reverse):
    chromosome, start, end, name = k
    return {
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "forward": forward,
        "reverse": reverse,
        "name": name
    }

def footprint(bam: str, bed: str, assembly: str = "hg38", w: int = 500, dnase: bool = False, bias_type = "SH"):

    # load HMM and bias parameters for ATAC-seq
    g = GenomeData(organism = assembly)
    hmm_data = HmmData()
    if dnase:
        hmm_file = hmm_data.get_default_hmm_dnase_bc()
        if bias_type == 'SH':
            table_F = hmm_data.get_default_bias_table_F_SH()
            table_R = hmm_data.get_default_bias_table_R_SH()
            bias_table = BiasTable().load_table(table_file_name_F = table_F, table_file_name_R = table_R)
        elif bias_type == 'DH':
            table_F = hmm_data.get_default_bias_table_F_DH()
            table_R = hmm_data.get_default_bias_table_R_DH()
            bias_table = BiasTable().load_table(table_file_name_F = table_F, table_file_name_R = table_R)
    else:
        hmm_file = hmm_data.get_default_hmm_atac_paired()
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table = BiasTable().load_table(table_file_name_F = table_F, table_file_name_R = table_R)

    # load reads from BAM
    reads_file = GenomicSignal(bam)
    reads_file.load_sg_coefs(SG_WINDOW_SIZE)

    # open data and sequence
    bam = Samfile(bam, "rb")
    fasta = Fastafile(g.get_genome())

    # load and expand regions
    with open(bed, 'r') as f:
        regions = [ expandRegion(*tuple(line.strip().split()[:3]), line.strip().split()[3] if len(line.strip().split()) >= 4 else None, w) for line in f ]
    
    # load signal
    forward = []; reverse = []
    for i, x in enumerate(regions):
        chromosome, start, end, _ = x
        atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r = reads_file.get_signal_atac(
            chromosome, start, end, 0, 0, FORWARD_SHIFT, REVERSE_SHIFT,
            50, 98, 98, bias_table, g.get_genome()
        )
        forward.append(atac_norm_f)
        reverse.append(atac_norm_r)
        if i % 500 == 0: print("INFO: aggregating region %d of %d" % (i, len(regions)), file = sys.stderr)

    return [ regionDict(regions[i], forward[i], reverse[i]) for i in range(len(regions)) ]
