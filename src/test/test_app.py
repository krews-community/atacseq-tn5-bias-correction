#!/usr/bin/env python3

import os
import tempfile
import unittest
import math
import hashlib

INPUTS = os.path.join( os.path.dirname(os.path.realpath(__file__)), "resources" ) + ":/input"
GENOME = os.path.join( os.path.dirname(os.path.realpath(__file__)), "resources", "hg38-chrM.tar.gz" )

class TestApp(unittest.TestCase):
        
    def assertFileExists(self, f: str):
        self.assertEqual(os.path.exists(f), True)

    def assertMD5(self, f: str, c: str):
        self.assertFileExists(f)
        m = hashlib.md5()
        with open(f, 'rb') as f:
            m.update(f.read())
        self.assertEqual(m.hexdigest(), c)

    def test_json(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/root/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertMD5("{d}/test.json".format(d = d), "30eba49fc846259e722c6972767abd90")

    def test_aggregate_json(self):
        
        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/root/rgtdata/hg38-chrM test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --aggregate > {d}/test.json
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertMD5("{d}/test.json".format(d = d), "0321578693d15d67341148a3c68641cd")
    
    def test_plot(self):

        with tempfile.TemporaryDirectory() as g:

            if os.system("tar zfx {GENOME} --directory {g}".format(GENOME = GENOME, g = g)) != 0:
                raise Exception("unable to extract required genome files")

            with tempfile.TemporaryDirectory() as d:
                if os.system("""
                    docker run --volume {inputs} --volume {g}:/root/rgtdata/hg38-chrM --volume {d}:/output test python3 -m app.main \
                        --bed /input/test.bed --bam /input/test.bam --assembly hg38-chrM --plot-output /output/test.svg > /dev/null
                """.format(inputs = INPUTS, d = d, g = g)) != 0:
                    raise Exception("unable to run tests")
                self.assertFileExists("{d}/test.svg".format(d = d))
