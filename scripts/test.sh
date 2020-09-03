#!/bin/bash
set -e

# cd to project root directory
cd "$(dirname "$(dirname "$0")")"
docker build -t test .

# change to package directory, run tests
cd src
python3 -m unittest test.test_app
