import os
import numpy as np
import shutil
from shutil import ignore_patterns
import stat
import glob
import pdb

def merge_single_sim( sim_folders):

    if len(sim_folders) == 1:
        return

    # Specify file patterns to ignore and those to ignore just for the latter copy operations
    ext_ignore = ["*.log","*.err","*.csv","*.png","*.pdf","*.jobscript",".DS_Store","*.mp4","*pvtp.pvd"]

    # Get path of first folder. This is where the merge will be created.
    head_tail = os.path.split(sim_folders[0])

    # Begin by creating a new folder
    merge_fold = os.path.join( head_tail[0], 'merge')
    if os.path.exists(merge_fold):
        shutil.rmtree(merge_fold)
    os.mkdir( merge_fold)
    print('Created merge folder: {}'.format(merge_fold))

    for idx,src in enumerate( sim_folders):

        # For the first folder, just copy almost every thing except those things to ignore
        if idx == 0:
            if not os.path.exists(merge_fold) or not os.listdir(merge_fold):
                copytree(sim_folders[0], merge_fold, ignore=ignore_patterns(*ext_ignore))
            else:
                raise Exception('There is already stuff inside merge folder')

        # For the latter folders
        else:

            # Find last index in results of merge
            # find result folders inside the main result folder
            res_folders = glob.glob( os.path.join(merge_fold, 'result/result*-*'))
            res_folders = sorted( res_folders, key=fold_key)
            # search inside last folder
            idx_files = glob.glob( os.path.join(res_folders[-1], 'SylinderAscii*.dat'))
            idx_files = sorted( idx_files, key=file_key)
            max_key = max([file_key(i) for i in idx_files])

            # Search in src folder
            # find result folders inside the main result folder
            src_files = glob.glob( os.path.join(src, 'result/result*-*/*'))
            src_files = sorted( src_files, key=file_key)
            src_files = remove_first_two_files(src_files)

            resdest = res_folders[-1]
            # pdb.set_trace()
            for fil in src_files:

                # if key=0, skip it
                if file_key( fil) == 0:
                    continue

                # check if key exceeds folder bound. If yes, switch to new folder
                if check_key_bound_folder( max_key+file_key(fil), resdest):
                    resdest = get_new_bound_folder( resdest)
                    print('Generated new destination folder: {}'.format(resdest))
                    os.mkdir( resdest)

                dstfile = get_dstfile_with_key_offset( fil, resdest, max_key)
                # copy each file from src results to merge results after updating its index
                shutil.copyfile( fil, dstfile)

                # check if file is .pvtp file
                head_tail = os.path.split(dstfile)
                if os.path.splitext(head_tail[1])[1] == '.pvtp':
                    # read file
                    with open(dstfile, 'r') as ff:
                        lines = ff.readlines()
                        for idx, line in enumerate( lines):
                            if line.startswith('<Piece Source="'):
                                # Change source in this line
                                str_2_change = line[15:-8]
                                str_pre = '_'.join( str_2_change.split('_')[:-1] )
                                newkey = int(str_2_change.split('_')[-1])+max_key-1
                                str_new = str_pre + '_' + str(newkey)
                                lines[idx] = line[0:15] + str_new + line[-8::]

                    # Write lines to file
                    with open(dstfile, 'w') as ff:
                        ff.writelines( lines)

# Define an updated copytree function
def copytree(src, dst, symlinks = False, ignore = None):
    if not os.path.exists(dst):
        os.makedirs(dst)
        shutil.copystat(src, dst)
    lst = os.listdir(src)
    if ignore:
        excl = ignore(src, lst)
        lst = [x for x in lst if x not in excl]
    for item in lst:
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if symlinks and os.path.islink(s):
            if os.path.lexists(d):
                os.remove(d)
            os.symlink(os.readlink(s), d)
            try:
                st = os.lstat(s)
                mode = stat.S_IMODE(st.st_mode)
                os.lchmod(d, mode)
            except:
                pass # lchmod not available
        elif os.path.isdir(s):
            copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def file_key(f):
    k = int( f.split("_")[-1].split(".")[0])
    return k

def fold_key(f):
    k = int( f.split("-")[-1])
    return k

def replace_key(f,nk):
    fnew = "".join( idx_files[0].split("_")[:-1] ) + "_" + str(nk) + '.' + idx_files[0].split("_")[-1].split(".")[-1]
    return fnew

def check_key_bound_folder( k, fold):
    fk = fold_key( fold)
    if k > fk:
        return True
    else:
        return False

def get_new_bound_folder( old_fold):
    upper = fold_key( old_fold)
    lower = int(old_fold.split("/result")[-1].split("-")[0])
    newlower = upper+1
    newupper = upper+1+ (upper-lower)
    new_fold = os.path.join( "/".join( old_fold.split("/")[:-1]), "result") + str(newlower) + "-" + str(newupper)
    return new_fold

def get_dstfile_with_key_offset( srcfile, dstpath, offset):

    old_key = file_key( srcfile)
    new_key = old_key + offset-1
    tag = "_".join(srcfile.split("/")[-1].split("_")[:-1])
    ext = srcfile.split(".")[-1]
    dstfile = os.path.join( dstpath, tag + "_" + str(new_key) + "." + ext)
    return dstfile

def remove_first_two_files( srcfiles):
    new_list = []
    for fil in srcfiles:
        if file_key(fil) >= 2:
            new_list.append(fil)
    return new_list

def sort_sim_names(s):
    return len(s.split('/')[-1].split('_c'))

if __name__ == "__main__":
    
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

    # Alternate
    sim_sub_paths = [
            "/Users/saadjansari/Documents/Projects/Results/AMSOS/Tactoids/scan_filamin_6400/run/f5",
            "/Users/saadjansari/Documents/Projects/Results/AMSOS/Tactoids/scan_filamin_6400/run/f5_c",
            ]
    merge_single_sim( sim_sub_paths)
            

