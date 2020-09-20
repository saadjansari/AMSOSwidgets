import os, time, yaml, glob, copy
import math, shutil, random
import numpy as np
import matplotlib.pyplot as plt
import pdb
import miniball

from util_func import *
from decorators import timer

def analyzeSims( spaths):
    # analyze each sim

    for spath in spaths:
        print(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -') 
        print('Analyzing Sim:\nPath: {}'.format(spath))

        analyzeSim(spath)

        print(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -') 

@timer
def analyzeSim( tpath):

    params = preInit(tpath)

    ftlist, ptlist, params = initSim( params)
    nT = params['nT']
    nF = params['nF']
    nP = params['nP']
    AC = params['AC']
    print('Filaments: {0}\nXlinks: {1}\nFrames: {2}'.format( nF,nP,nT))

    # Analysis
    nT = params['nT']
    nF = params['nF']

    # cct = []
    cct2 = np.zeros( (nT,nF),dtype=bool)
    for it in range(nT):

        # Create a graph for filaments
        g = Graph( len(ftlist[it]) ); 

        # add edges to the graph (each edge represents a binding xlink)
        for p in ptlist[it]:
            if p.link0 != -1 and p.link1 != -1:
                g.addEdge(p.link0, p.link1) 

        # find connected components
        cc = g.connectedComponents()
        cct2[[it for i in range(len(cc[0]))], cc[0]] = True

    ftlist_trimmed = []
    for row in cct2:
        fff = []
        idx = np.where(row)[0].tolist()
        for el in idx:
            fff.append( ftlist[it][el])
        ftlist_trimmed.append( fff)

    # Find minimal bounding sphere
    # For filaments in ftlist_trimmed, get their plus and minus ends
    # diam = []
    # for flist in ftlist_trimmed:
        # plusends = np.zeros( (len(flist),3))
        # minusends = np.zeros( (len(flist),3))
        # for idx,fila in enumerate(flist):
            # minusends[idx, :] = fila.pos0
            # plusends[idx, :] = fila.pos1
        # ends = np.vstack( (minusends, plusends) ) 
        # C, r2 = miniball.get_bounding_ball(ends)
        # diam.append( 2*np.sqrt(r2))
        
    # plot Filament number
    fig,ax = plt.subplots()
    nft = [np.sum(ii) for ii in cct2]
    ax.plot( params['dt']*np.arange(nT) , nft, color='red', linewidth=2)
    ax.set(xlabel='time (s)', ylabel='Count', title='Number of filaments inside aster')
    fig.savefig( os.path.join( params['rpath'],'numberInAster.pdf') , bbox_inches='tight', dpi=600)
    plt.close('all')

    # plot bounding sphere diameter
    # fig,ax = plt.subplots()
    # ax.plot( params['dt']*np.arange(nT) , np.array( diam), color='red', linewidth=2)
    # ax.set(xlabel='time (s)', ylabel='Diameter ($\mu m$)', title='Aster diameter')
    # fig.savefig( os.path.join( params['rpath'],'diameterAster.pdf') , bbox_inches='tight', dpi=600)
    # plt.close('all')

# preInit {{{
@timer
def preInit(tpath):

    # Create plots directory
    rpath = os.path.join(tpath, 'plots')
    if os.path.exists(rpath):
        shutil.rmtree( rpath)
    os.mkdir( rpath)

    # Specify parameters
    analyze_bulk = False

    # Load parameters from RunConfig.yaml and ProteinConfig.yaml files
    with open( os.path.join(tpath,'RunConfig.yaml')) as f:
        pars1 = yaml.load( f, Loader=yaml.FullLoader)
    with open( os.path.join(tpath,'ProteinConfig.yaml')) as f:
        pars2 = yaml.load( f, Loader=yaml.FullLoader)
        
    dt = pars1['timeSnap'] # Specify frame step of simulation
    boxsize = np.array( pars1['initBoxHigh']) # Specify frame step of simulation
        
    # title of plots
    ftitle='Aster analysis'

    # Display
    AC = AnalysisBook(boxsize)

    params = {
            'tpath':tpath,
            'rpath':rpath,
            'dt':dt,
            'boxsize':boxsize,
            'ftitle':ftitle,
            'AC':AC
            }
    return params
#}}}

# initSim {{{
@timer
def initSim( params):

    # Get key to sort through file indices
    def fileKey(f):
        k = int(f[f.rfind("_")+1:f.rfind(".")])
        return k

    # find filament files
    files = glob.glob( os.path.join(params['tpath'], 'result/result*/SylinderAscii_*.dat'))
    files = sorted(files, key=fileKey)
    ftlist = []
    for fil in files:
        flist = []
        with open(fil, 'r') as file1:
            for line in file1:
                if line.startswith('C'):
                    data = line.split()
                    gid = int(data[1])
                    radius = float(data[2])
                    pos0 = np.array([float(data[3]), float(data[4]), float(data[5])])
                    pos1 = np.array([float(data[6]), float(data[7]), float(data[8])])
                    flist.append( Filament(pos0, pos1, radius,gid))
        ftlist.append(flist)

    # Read protein file
    files = glob.glob( os.path.join(params['tpath'], 'result/result*/ProteinAscii_*.dat'))
    files = sorted(files, key=fileKey)
    ptlist = []
    for fil in files:
        plist = []
        with open(fil, 'r') as file2:
            for line in file2:
                if line.startswith('P'):
                    data = line.split()
                    gid = int(data[1])
                    pos0 = np.array([float(data[3]), float(data[4]), float(data[5])])
                    pos1 = np.array([float(data[6]), float(data[7]), float(data[8])])
                    link0 = int(data[9])
                    link1 = int(data[10])
                    plist.append( Protein(pos0, pos1, link0, link1, gid))
        ptlist.append(plist)

    # Toss out zeroth timestep
    ftlist.pop(0)
    ptlist.pop(0)

    params['nT'] = len(ftlist)
    params['nF'] = len(ftlist[0])
    params['nP'] = len(ptlist[0])

    return ftlist, ptlist, params
# }}}

if __name__ == "__main__":

    paths = [
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/test_free", 
            ]
    analyzeSims( paths)













