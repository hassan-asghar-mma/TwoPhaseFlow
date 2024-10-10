#!/usr/bin/env bash

# This script extracts contact angle data from the solver log
# and combines it with the wetted area data into a properly
# formatted CSV file.

# Extract contact angle data from solver log
touch contactAngleData.csv
echo "dummyCol, time, caMin, caMean, caDev, caMax, nCLcells" > contactAngleData.csv

grep "+++" *growingDroplet3D.o* >> contactAngleData.csv

# Combine contact angle data and wetted area data into single
# CSV file
python3 mergeCaseData.py

rm -f contactAngleData.csv

