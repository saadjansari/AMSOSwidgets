#!/usr/bin/env python3

from pathlib import Path
import os
import argparse
import glob
import pdb

from sim import Sim


'''
Name: aman.py
Description: aman ( am(sos) an(alysis) ) performs frame-by-frame analysis sequentially on simulations
Input: To see type aman.py -h
'''


def parseArgs():

    parser = argparse.ArgumentParser(prog='aman.py')
    parser.add_argument(
        '-T',
        '--tactoid',
        action='store_true',
        default=False,
        help='use tactoid options for analysis')
    parser.add_argument(
        '-C',
        '--confine',
        action='store_true',
        default=False,
        help='use confinement options for analysis')
    parser.add_argument(
        '-G',
        '--graph',
        action='store_true',
        default=False,
        help='make graphs for sims')
    parser.add_argument(
        '--overwrite',
        action='store_true',
        default=False,
        help='overwrite feather files')
    opts = parser.parse_args()

    # Specify analysis options
    if opts.tactoid:
        opts.analyze_cluster = True
        opts.analyze_global_order = True
        opts.analyze_local_order = True
        opts.analyze_xlinks = True
        opts.analyze_pairpair_separation = False
        opts.analyze_aspect_ratio = True
        opts.analyze_z_ordering = False
        opts.length_distribution = False

    elif opts.confine:
        opts.analyze_cluster = True
        opts.analyze_global_order = True
        opts.analyze_local_order = True
        opts.analyze_xlinks = False
        opts.analyze_pairpair_separation = False
        opts.analyze_aspect_ratio = False
        opts.analyze_z_ordering = True
        opts.length_distribution = False

    if opts.tactoid and opts.confine:
        raise Exception('Cannot specify both tactoid and confinement opts')
    elif not opts.tactoid and not opts.confine:
        raise Exception('Please specify either tactoid or confinement opts')

    if opts.tactoid:
        print('Configuration: Tactoid')
    elif opts.confine:
        print('Configuration: Confine')
    print('Analysis flags: ')
    print('\tCluster: {}'.format(opts.analyze_cluster))
    print('\tGlobal order: {}'.format(opts.analyze_global_order))
    print('\tXlinks: {}'.format(opts.analyze_xlinks))
    print('\tPair-pair separation: {}'.format(opts.analyze_pairpair_separation))
    print('\tAspect Ratio: {}'.format(opts.analyze_aspect_ratio))
    print('\tZ-ordering: {}'.format(opts.analyze_z_ordering))
    return opts


def analyze(opts):

    # Prompt user for relative path to folder containing sims
    prompt = '\nSpecify relative path to run folder with simulation folders: '
    relpath = Path(input(prompt))  # get path from user
    fpath = Path(relpath).resolve()
    if not fpath.exists():
        raise Exception('specified path does not exist')
    print('Run path: {}\n'.format(fpath))

    # get sim folders (folders with result dirs) recursively using ** notation
    spaths = [res_path.parent for res_path in relpath.glob('**/result')]
    snames = ['__'.join(path.name.split('/')) for path in spaths]
    print(snames)

    # For each sim, analyze it
    for idx in range(len(spaths)):
        sim = Sim(spaths[idx], snames[idx], opts)
        sim.analyze()


if __name__ == "__main__":

    # Preamble for display
    # ****************************************************************
    # ****************************************************************
    # ****************************************************************
    filename = os.path.basename(__file__) + '.py'
    hlen = 80
    hlen2 = int((hlen - len(filename) - 4) / 2)
    symbol = '*'
    print(symbol * hlen)
    print(symbol * hlen)
    print(symbol * hlen2 + '  ' + filename.upper() + '  ' + symbol * hlen2)
    print(symbol * hlen)
    print(symbol * hlen + '\n')
    # ****************************************************************
    # ****************************************************************
    # ****************************************************************

    opts = parseArgs()
    analyze(opts)
