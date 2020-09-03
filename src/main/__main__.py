#!/usr/bin/env python

import os
import sys
import json

from argparse import ArgumentParser

from ..footprint.footprint import footprint
from ..regions.filter import FilteredRegions
from ..plot.plot import plot

def args():
    parser = ArgumentParser()
    parser.add_argument("--bam", type = str, help = "path to alignments in BAM format", required = True)
    parser.add_argument("--bed", type = str, help = "path to regions across which to compute signal", required = True)
    parser.add_argument("--assembly", type = str, help = "genomic assembly to use", default = "hg38")
    parser.add_argument("--plot-output", type = str, help = "if provided, saves an aggregate plot to this path", default = None)
    parser.add_argument("--font", type = str, help = "if set, path to a font to use during plotting", default = None)
    parser.add_argument("--ext-size", type = int, help = "expands regions by the given number of basepairs around their centers", default = 500)
    parser.add_argument("--aggregate", action = "store_true", help = "if set, outputs aggregate signal rather than profiles for each region", default = False)
    parser.add_argument("--occurrence-threshold", type = float, help = "specificies that the given BED file contains FIMO occurrences which should be filtered at this q-value.", default = None)
    return parser.parse_args()

def main():

    cArgs = args()

    if not os.path.exists("/root/rgtdata/{assembly}/genome_{assembly}.fa".format(assembly = cArgs.assembly)):
        print(
            "WARNING: genomic data is not present for {assembly}. We will attempt to download it.".format(assembly = cArgs.assembly),
            file = sys.stderr
        )
        print(
            "If you are running many jobs, they might run faster if you mount the appropriate data at /root/rgtdata/{assembly}.".format(assembly = cArgs.assembly),
            file = sys.stderr
        )
        result = os.system("python3 /reg-gen/data/setupGenomicData.py --{assembly}".format(assembly = cArgs.assembly))
        if result != 0:
            print("FATAL: Unable to load genome data for {assembly}.".format(assembly = cArgs.assembly), file = sys.stderr)
            return 1

    if cArgs.occurrence_threshold is None:
        signal = footprint(cArgs.bam, cArgs.bed, cArgs.assembly)
    else:
        with FilteredRegions(cArgs.bed, cArgs.occurrence_threshold) as b:
            signal = footprint(cArgs.bam, b.name, cArgs.assembly)

    if cArgs.aggregate or cArgs.plot_output is not None:

        forward = [ 0 for _ in range(cArgs.ext_size * 2) ]
        reverse = [ 0 for _ in range(cArgs.ext_size * 2) ]
        for x in signal:
            if x["forward"] is not None:
                for i, xx in enumerate(x["forward"]): forward[i] += xx
            if x["reverse"] is not None:
                for i, xx in enumerate(x["reverse"]): reverse[i] += xx
        if cArgs.aggregate: signal = { "forward": forward, "reverse": reverse }

        if cArgs.plot_output is not None:
            plot(forward, reverse, cArgs.font, cArgs.plot_output)

    print(json.dumps(signal))

    return 0

if __name__ == "__main__":
    sys.exit(main())
