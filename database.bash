#!/usr/bin/env bash

### capture analysis name
analysis=$1
echo 'analysis: ' $analysis

make -C blastdb $analysis.phr

