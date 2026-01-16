#!/bin/bash 
source /cvmfs/sbnd.opensciencegrid.org/products/sbnd/setup_sbnd.sh
source /exp/sbnd/app/users/vito/products/setup.sh

setup sbndcode v10_06_03 -q e26:prof

export BEARER_TOKEN_FILE=/tmp/bt_u$(id -u)
htgettoken -v -a htvaultprod.fnal.gov -i sbnd
