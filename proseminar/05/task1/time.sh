#!/bin/bash

echo "Vary the number of bodies with const number of steps=10000"
echo "10"
./2d_nBody 10
echo "100"
./2d_nBody 100
echo "500"
./2d_nBody 500
echo "1000"
./2d_nBody 1000

echo "Vary the number of steps with const bodies=100"
echo "100"
./2d_nBody 100 100
echo "1000"
./2d_nBody 100 1000
echo "10000"
./2d_nBody 100 10000
echo "100000"
./2d_nBody 100 100000
