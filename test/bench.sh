#!/bin/bash

for i do
JULIA_NUM_THREADS=$i julia -O3 runbenchmarks.jl
done
