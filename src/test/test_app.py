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
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
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
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --aggregate > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j["all"]["forward"]), 1000)
                    self.assertEqual(len(j["all"]["reverse"]), 1000)

    def test_aggregate_json_extsize(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --aggregate --ext-size 20 > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j["all"]["forward"]), 40)
                    self.assertEqual(len(j["all"]["reverse"]), 40)

    def test_aggregate_dnase(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --dnase --aggregate > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j["all"]["forward"]), 1000)
                    self.assertEqual(len(j["all"]["reverse"]), 1000)

    def test_aggregate_dnase_dh(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --dnase --bias-type DH --aggregate > {d}/test.json
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
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM --volume {d}:/output test python3 -m app.main \
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
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
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

    def test_occurrences_1F(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.occ.bed --bam /input/test.bam --assembly hg38-chrM --occurrence-threshold 1.0 > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j), 9)

    def test_occurrences_extras(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.occ.extras.bed --bam /input/test.bam --assembly hg38-chrM --occurrence-threshold 1.0 > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.json".format(d = d))
                with open("{d}/test.json".format(d = d), 'r') as f:
                    j = json.load(f)
                    self.assertEqual(len(j), 1)

    def test_occurrences_extras_tsv(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --env RGTDATA=/rgtdata --volume {inputs} --volume {g}:/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.occ.extras.bed --bam /input/test.bam --assembly hg38-chrM --occurrence-threshold 1.0 \
                        --output-as-tsv --aggregate > {d}/test.tsv
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.tsv".format(d = d))
                with open("{d}/test.tsv".format(d = d), 'r') as f:
                    self.assertEqual(len([ x for x in f ]), 1)
                with open("{d}/test.tsv".format(d = d), 'r') as f:
                    self.assertEqual(len(f.readline().strip().split()), 3)               
