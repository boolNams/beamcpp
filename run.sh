#!/bin/bash
echo "COMPILE..."
rm main.o func.o prog
make
echo "RUNNING..."
./prog
