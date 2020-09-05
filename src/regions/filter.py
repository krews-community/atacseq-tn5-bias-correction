#!/usr/bin/env python

import tempfile

class FilteredRegions:

    def __init__(self, path: str, threshold: float):
        self.path = path
        self.threshold = threshold
    
    def __enter__(self):
        self.tempfile = tempfile.NamedTemporaryFile('wt')
        with open(self.path, 'r') as f:
            f.readline()
            for line in f:
                if float(line.strip().split()[-1]) < self.threshold:
                    self.tempfile.write('\t'.join(line.strip().split()[1:4]) + '\t' + line.strip().split()[0] + '\n')
        self.tempfile.flush()
        return self.tempfile
    
    def __exit__(self, *args):
        self.tempfile.close()
