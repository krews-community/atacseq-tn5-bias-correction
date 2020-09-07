#!/usr/bin/env python

import sys
import tempfile

class FilteredRegions:

    def __init__(self, path: str, threshold: float, fa: str, chromSizes: str):
        self.path = path
        self.threshold = threshold
        with open(chromSizes, 'r') as f:
            self.chromosomes = set([ x.strip().split()[0] for x in f ])
        with open(fa, 'r') as f:
            self.chromosomes = self.chromosomes.intersection(set([ x.strip().split('>')[1] for x in f if x[0] == '>' ]))
    
    def __enter__(self):
        self.tempfile = tempfile.NamedTemporaryFile('wt')
        with open(self.path, 'r') as f:
            f.readline()
            for line in f:
                if float(line.strip().split()[-1]) < self.threshold and line.strip().split()[1] in self.chromosomes:
                    self.tempfile.write('\t'.join(line.strip().split()[1:4]) + '\t' + line.strip().split()[0] + '\n')
        self.tempfile.flush()
        return self.tempfile
    
    def __exit__(self, *args):
        self.tempfile.close()
