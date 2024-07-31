#!/usr/bin/env python
import os
import shutil

# Define output directory and ensemble size
DARTDIR = '/home/cchang/DART/DA-TEM/DA_system/models/dvmdostem/work' # output directory
BKGDDIR = 'exp02/bkgd'          # model ens bkgd directory
ANALDIR = 'exp02/anal'          # model analysis ens directory
ensnum = 40
infile = 'restart_input.txt'    # input file name list
outfile = 'restart_output.txt'  # output file name list

# open a text file to write on ensemble INPUT file list
with open(infile, 'w') as f:
    # loop over ensembles
    for ens in range(1, ensnum + 1):
        # define filename
        filename = f"{BKGDDIR}/model_bkgd_{ens:02d}.nc"
        # write out filename
        f.write(filename + '\n')

print(f" file {infile} has been generated!")

# open a text file to write on ensemble file names
with open(outfile, 'w') as f:
    # loop over ensembles
    for ens in range(1, ensnum + 1):
        # define filename
        filename = f"{ANALDIR}/model_anal_{ens:02d}.nc"
        # write out filename
        f.write(filename + '\n')

print(f" file {outfile} has been generated!")


# move to DARTDIR
dest= f"{DARTDIR}/{infile}"
shutil.move(infile, dest)

print(f" file {infile} move to {dest}")

dest= f"{DARTDIR}/{outfile}"
shutil.move(outfile, dest)

print(f" file {outfile} move to {dest}")

