#!/usr/bin/env bash

DIR1=dir1
DIR2=dir2
LOGS=logs.txt
A=1
X=2
O=2
E=1
TIMEOUT=12h

rm $LOGS || true
for F in `ls -S -r $DIR1`; do
    for MEM in -m -s; do
        for MATCH in -c -t; do
            timeout -s 2 $TIMEOUT /usr/bin/time -v ./benchmark -a $A -x $X -o $O -e $E $MEM $MATCH $DIR1/$F $DIR2/$F 2>&1 | tee -a $LOGS
        done
    done
done
sudo poweroff