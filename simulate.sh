#!/bin/bash

# Rscript simulate.R -b 0.005 -m 0.005 -r 0.8 -N 500 -t 50000 -o 7 -v
Rscript simulate.R -b 0.5 -m 0.005 -r 0.8 -N 500 -t 50000 -o 7 -v
Rscript simulate.R -b 0.005 -m 0.5 -r 0.8 -N 500 -t 50000 -o 7 -v
