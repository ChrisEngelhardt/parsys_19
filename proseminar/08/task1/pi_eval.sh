#!/bin/bash
for threads in {1..8}
do
for pSize in 100000 1000000 10000000 100000000
    do
    ./pi $pSize $threads
    done
done