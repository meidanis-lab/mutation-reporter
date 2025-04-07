#!/usr/bin/env bash

for i in `ls $1`; do
    echo $i
    ls -l $i | cut -d\  -f5
    time ./mutation.bash $i $2
done
