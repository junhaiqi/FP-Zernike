import pymol
from pymol import *
import os
import sys
import argparse


def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', type=str, required=True,
                        help='input file, it is a pdb file.')
    parser.add_argument('-o', type=str, required=True,
                        help='output file, it is a *.wrl or *.obj(only surface mode).')
    parser.add_argument('-mode', type=str, required=True,
                        help='protein representation type. The following options are available: surface, mesh, sphere')
    args = parser.parse_args()
    return args


def main(input_pdb_file, output_file, mode="surface"):
    cmd.load(input_pdb_file)
    # cmd.hide("everything")
    cmd.show(mode)
    cmd.save(output_file)


if __name__ == "__main__":
    paras = initialization_parameters()
    main(paras.i, paras.o, paras.mode)
