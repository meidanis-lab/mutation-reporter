#!/usr/bin/env bash

### capture analysis name
analysis=$1
echo 'analysis: ' $analysis

make -C CONFIG/ blastdb/$analysis.phr

