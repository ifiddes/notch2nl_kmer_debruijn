#!/bin/bash


if [ -e ./jobTree ]; then
    rm -rf ./jobTree
fi

maxThreads=4
batchSystem=parasol
defaultJobMemory=8589934592

#Set the python path to just these local directories
export PYTHONPATH=:./:${PYTHONPATH}
#Preferentially put the local binaries at the front of the path
export PATH=:./jobTree/bin:${PATH}: