#!/usr/bin/env python3

import os
import tempfile
import unittest
import math
import hashlib
import json

INPUTS = os.path.join( os.path.dirname(os.path.realpath(__file__)), "resources" ) + ":/input"
GENOME = os.path.join( os.path.dirname(os.path.realpath(__file__)), "resources", "hg38-chrM.tar.gz" )

class TestApp(unittest.TestCase):
        
    def assertFileExists(self, f: str):
        self.assertEqual(os.path.exists(f), True)

    def test_json(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test RGTDATA=/rgtdata python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j), 7)

    def test_aggregate_json(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test RGTDATA=/rgtdata python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --aggregate > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j["all"]["forward"]), 1000)
                    self.assertEqual(len(j["all"]["reverse"]), 1000)
    
    def test_plot(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/rgtdata/hg38-chrM --volume {d}:/output test RGTDATA=/rgtdata python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --plot-output /output/test.svg > /dev/null
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.svg".format(d = d))

    def test_occurrences(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test RGTDATA=/rgtdata python3 -m app.main \
                        --bed /input/test.occ.bed --bam /input/test.bam --assembly hg38-chrM --aggregate --occurrence-threshold 1e-5 > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j), 3)
                    self.assertEqual(len(j["WTTTCTCTCWGTGYA"]["forward"]), 1000)
                    self.assertEqual(len(j["WTTTCTCTCWGTGYA"]["reverse"]), 1000)
                    self.assertEqual(len(j["ATTTCTCTCWGTGYA"]["forward"]), 1000)
                    self.assertEqual(len(j["ATTTCTCTCWGTGYA"]["reverse"]), 1000)
                    self.assertEqual(len(j["all"]["forward"]), 1000)
                    self.assertEqual(len(j["all"]["reverse"]), 1000)
