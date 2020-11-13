#!/usr/bin/env python3

import os
import numpy as np
import shutil
from shutil import ignore_patterns
from pathlib import Path
import stat
import pdb

'''
Name: merge_sims.py
Description: merge_sims performs a merge operation on a single set of simulations 
Input: To see help type merge_sims.py -h
'''

def parse_args():

    parser = argparse.ArgumentParser(prog='merge_sims.py')
    opts = parser.parse_args()
    return opts

def get_sim_names(opts=None):
    # Prompt user for sim path and list of sim names to merge

    # Prompt user for path to parent folder containing sims
    print('Get simulations to merge\n')
    prompt1 = 'Specify absolute path to parent containing simulations to merge:\n'
    relpath = Path(input(prompt1))  # get path from user
    fpath = Path(relpath).resolve()
    if not fpath.exists():
        raise Exception('specified path does not exist')

    print('\nSpecify simulation names to merge: (type stop when no more left to add)')
    sname = 'none'
    spaths  = []
    snames = []
    while len(snames) < 7:
        sname = input('Next sim-name : ')  # get path from user
        if sname == 'stop':
            break
        else:
            spath = fpath / sname
            if not spath.exists():
                raise Exception('Specified path does not exist')
            spath.resolve()
            snames.append( sname)
            spaths.append( spath)
            print('\nCurrent sims to merge: {}'.format(snames))
    
    print('\nFinal sims to merge: {}'.format(snames))
    print('*****************************************************************')
    return spaths

def merge_single_sim(sim_folders):
    # Merge a single sim

    if len(sim_folders) == 1:
        return

    # Specify file patterns to ignore and those to ignore just for the latter
    # copy operations
    ext_ignore = [
        "*.log",
        "*.err",
        "*.csv",
        "*.png",
        "*.pdf",
        "*.jobscript",
        ".DS_Store",
        "*.mp4",
        "*pvtp.pvd"]

    # Get path of first folder. This is where the merge will be created.
    head_tail = sim_folders[0].parent

    # Begin by creating the merge folder
    prompt = 'Specify name of merge folder (default "merge") : '
    mname = input(prompt)  # get path from user
    if mname == '':
        mname = 'merge'
    merge_dir = head_tail / mname

    if merge_dir.exists():
        shutil.rmtree(merge_dir)
    merge_dir.mkdir()
    print('Created merge folder: {}'.format(merge_dir))

    for idx, src in enumerate(sim_folders):

        # For the first folder, just copy almost every thing except those
        # things to ignore
        if idx == 0:
            if not merge_dir.exists() or not any(merge_dir.iterdir()):
                copytree(
                    sim_folders[0],
                    merge_dir,
                    ignore=ignore_patterns(
                        *ext_ignore))
            else:
                raise Exception('There is already stuff inside merge folder')

        # For the latter folders
        else:

            # Find last index in results of merge
            # find result folders inside the main result folder
            res_folders = merge_dir.glob('result/result*-*')
            res_folders = sorted(res_folders, key=fold_key)
            # search inside last folder
            idx_files = res_folders[-1].glob('SylinderAscii*.dat')
            idx_files = sorted(idx_files, key=file_key)
            max_key = file_key(idx_files[-1])
            # max_key = max([file_key(i) for i in idx_files])

            # Search in src folder
            # find result folders inside the main result folder
            src_files = src.glob('result/result*-*/*')
            src_files = sorted(src_files, key=file_key)
            src_files = remove_first_two_files(src_files)

            resdest = res_folders[-1]
            # pdb.set_trace()
            for fil in src_files:

                # if key=0, skip it
                if file_key(fil) == 0:
                    continue

                # check if key exceeds folder bound. If yes, switch to new
                # folder
                if check_key_bound_folder(max_key + file_key(fil), resdest):
                    resdest = get_new_bound_folder(resdest)
                    print('Generated new destination folder: {}'.format(resdest))
                    resdest.mkdir()
                    # os.mkdir(resdest)

                dstfile = get_dstfile_with_key_offset(fil, resdest, max_key)
                # copy each file from src results to merge results after
                # updating its index
                shutil.copyfile(fil, dstfile)

                # check if file is .pvtp file
                # head_tail = os.path.split(dstfile)
                if dstfile.suffix == '.pvtp':
                    # read file
                    with dstfile.open('r') as ff:
                        lines = ff.readlines()
                        for idx, line in enumerate(lines):
                            if line.startswith('<Piece Source="'):
                                # Change source in this line
                                str_2_change = line[15:-8]
                                str_pre = '_'.join(
                                    str_2_change.split('_')[:-1])
                                newkey = int(str_2_change.split(
                                    '_')[-1]) + max_key - 1
                                str_new = str_pre + '_' + str(newkey)
                                lines[idx] = line[0:15] + str_new + line[-8::]

                    # Write lines to file
                    with open(dstfile, 'w') as ff:
                        ff.writelines(lines)

# Define an updated copytree function


def copytree(src, dst, symlinks=False, ignore=None):
    if not dst.exists():
        dst.mkdir()
        shutil.copystat(str(src), str(dst))
    lst = [l.name for l in src.iterdir()]
    if ignore:
        excl = ignore(str(src), lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = src / item
        d = dst / item
        if symlinks and s.is_symlink():
            if os.path.lexists(d):
                d.unlink()
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except BaseException:
                pass  # lchmod not available
        elif s.is_dir():
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)


def file_key(f):
    k = int(str(f).split("_")[-1].split(".")[0])
    return k


def fold_key(f):
    k = int(str(f).split("-")[-1])
    return k


def check_key_bound_folder(k, fold):
    fk = fold_key(fold)
    if k > fk:
        return True
    else:
        return False


def get_new_bound_folder(old_fold):
    upper = fold_key(old_fold)
    lower = int(str(old_fold).split("/result")[-1].split("-")[0])
    newlower = upper + 1
    newupper = upper + (newlower - lower)
    new_fold = old_fold.parent / "result{}-{}".format(newlower, newupper)
    return new_fold


def get_dstfile_with_key_offset(srcfile, dstpath, offset):

    old_key = file_key(srcfile)
    new_key = old_key + offset - 1
    tag = "_".join(str(srcfile).split("/")[-1].split("_")[:-1])
    ext = srcfile.suffix
    dstfile = dstpath / "{}_{}{}".format(tag, new_key, ext)
    return dstfile


def remove_first_two_files(srcfiles):
    new_list = []
    for fil in srcfiles:
        if file_key(fil) >= 2:
            new_list.append(fil)
    return new_list


def sort_sim_names(s):
    return len(str(s).split('/')[-1].split('_c'))


if __name__ == "__main__":

    print('*****************************************************************')
    print('*********************** merge_sims.py ***************************')
    print('*****************************************************************')
    # opts = parse_args()

    # # Get all sim folders
    # # Prompt user for relative path to folder containing sims
    # prompt = '\nSpecify relative path to run folder with simulation folders: '
    # relpath = input(prompt) # get path from user
    # fpath = os.path.abspath( relpath)
    # if not os.path.exists( fpath):
        # raise Exception('specified path does not exist')
    # print('Run path: {}\n'.format(fpath) )

    # # get sim folders
    # spaths = glob.glob( os.path.join( relpath, '*')) 
    # # for each sim folder, look for sub-directories
    # for spath in spaths:
        # sim_sub_paths = glob.glob( os.path.join( spath, '*') )

        # # remove merge sub-directory if it exists
        # sim_sub_paths = [ii for ii in sim_sub_paths if ii.split('/')[-1] != 'merge']
        # sim_sub_paths = sorted(sim_sub_paths, key=sort_sim_names)

        # # merge these dirs
        # if len(sim_sub_paths) != 1:
            # merge_single_sim( sim_sub_paths)

    spaths = get_sim_names() 
    merge_single_sim(spaths)
