#!/usr/bin/env python

import sys
import argparse
import math

from typing import List, Tuple
from pysam import Fastafile, Samfile

from rgt.Util import ErrorHandler, HmmData, GenomeData, OverlapType
from rgt.HINT.signalProcessing import GenomicSignal
from rgt.HINT.biasTable import BiasTable

from .constants import *

def expandRegion(chromosome, start, end, name = None, w = 500, strand = '.'):
    m = (math.ceil if strand == '-' else math.floor)((int(start) + int(end)) / 2)
    return chromosome, m - w, m + w, name, strand

def regionDict(k, forward, reverse):
    chromosome, start, end, name, strand = k
    return {
        "chromosome": chromosome,
        "start": start,
        "end": end,
        "forward": forward,
        "reverse": reverse,
        "name": name,
        "strand": strand
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
        regions = [ expandRegion(
            *tuple(line.strip().split()[:3]), line.strip().split()[3] if len(line.strip().split()) >= 4 else None, w,
            line.strip().split()[4] if len(line.strip().split()) >= 5 else '.'
        ) for line in f ]
    
    # load signal
    forward = []; reverse = []; failed = 0
    get_reads = reads_file.get_signal_atac if not dnase else reads_file.get_signal
    for i, x in enumerate(regions):
        try:
            chromosome, start, end, _, strand = x
            atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r = get_reads(
                chromosome, start, end, 0, 0, FORWARD_SHIFT if not dnase else 0, REVERSE_SHIFT if not dnase else 0,
                1000 if dnase else 150, 98, 98, bias_table, g.get_genome()
            )
            if strand == '-':
                atac_norm_f.reverse()
                atac_norm_r.reverse()
            forward.append(atac_norm_f if strand != '-' else atac_norm_r)
            reverse.append(atac_norm_r if strand != '-' else atac_norm_f)
            if i % 500 == 0: print("INFO: aggregating region %d of %d" % (i, len(regions)), file = sys.stderr)
        except:
            if len(forward) <= i: forward.append(None)
            if len(reverse) <= i: reverse.append(None)
            failed += 1
    if failed > 0:
        print("WARNING: failed to generate bias-corrected signal profiles for %d regions" % failed, file = sys.stderr)

    return [ regionDict(regions[i], forward[i], reverse[i]) for i in range(len(regions)) if forward[i] is not None and reverse[i] is not None ]
