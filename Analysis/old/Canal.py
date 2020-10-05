import os, time, yaml, glob, copy
import math, random, shutil
import numpy as np
import matplotlib.pyplot as plt
import pickle
from itertools import compress
import pdb

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
    dt = params['dt']
    print('Filaments: {0}\nXlinks: {1}\nFrames: {2}'.format( nF,nP,nT))

    # Analysis
    cct,ftlist_trimmed = getFilamentsInsideCluster( ftlist, ptlist, params, plotting=0)

    # Nematic Order
    sdat, sdat_bulk = calcNematicOrder( ftlist, ftlist_trimmed, params)

    # Local Polar Order
    pbulk_par,pbulk_apar,pcluster_par,pcluster_apar = calcLocalPolarOrder( ftlist, cct, params)

    # fig,ax = plt.subplots(figsize=(6,6))
    # ax.plot(dt*np.arange(nT), p1[0,:], linewidth=2, color='red', label='bulk')
    # ax.plot(dt*np.arange(nT), p2[0,:], linewidth=2, color='green', label='cluster')
    # ax.set(xlabel='Time(s)', ylabel='Local polar order', title='Local polar order',ylim=[-0.1,1.1])
    # plt.legend()
    # plt.show()
    # fig.savefig( os.path.join( params['rpath'],'localPolarOrder.pdf') , bbox_inches='tight', dpi=600)

    dat = {
            'S_cluster': sdat,
            'S_bulk': sdat_bulk,
            'P_local_bulk_par': pbulk_par,
            'P_local_bulk_apar': pbulk_apar,
            'P_local_cluster_par': pcluster_par,
            'P_local_cluster_apar': pcluster_apar,
            'dt': params['dt']*np.arange(nT),
            'params': params
            }
    pickle.dump( dat, open( os.path.join( params['tpath'], "data.pickle"), "wb" ) )


# preInit {{{
@timer
def preInit(tpath):

    # Create plots directory
    rpath = os.path.join(tpath, 'plots')
    if os.path.exists(rpath):
        shutil.rmtree( rpath)
    os.mkdir( rpath)

    # Load parameters from RunConfig.yaml and ProteinConfig.yaml files
    with open( os.path.join(tpath,'RunConfig.yaml')) as f:
        pars1 = yaml.load( f, Loader=yaml.FullLoader)
    with open( os.path.join(tpath,'ProteinConfig.yaml')) as f:
        pars2 = yaml.load( f, Loader=yaml.FullLoader)

    dt = pars1['timeSnap'] # Specify frame step of simulation
    boxsize = np.array( pars1['initBoxHigh']) # Specify frame step of simulation

    # get filament and protein number
    if os.path.exists( os.path.join(tpath, 'TubuleInitial.dat')):
        with open(os.path.join(tpath, 'TubuleInitial.dat'), 'r') as file1:
            nfil = len( file1.readlines() )-2
    else:
        nfil = pars1['sylinderNumber']
    if os.path.exists( os.path.join(tpath, 'ProteinInitial.dat')):
        with open(os.path.join(tpath, 'ProteinInitial.dat'), 'r') as file1:
            nprot = len( file1.readlines() )-2
    else:
        nprot = pars2['freeNumber']

    # Display
    AC = AnalysisBook(boxsize)

    params = {
            'name': os.path.basename(tpath),
            'tpath':tpath,
            'rpath':rpath,
            'dt':dt,
            'boxsize':boxsize,
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

# getFilamentsInsideCluster {{{
@timer
def getFilamentsInsideCluster( ftlist, ptlist, params, plotting=0):

    nT = params['nT']
    nF = params['nF']

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

    if plotting:
        fig,ax = plt.subplots()
        nft = [np.sum(ii)/nF for ii in cct2]
        ax.plot( params['dt']*np.arange(nT) , nft, color=params['cols']['cluster'], linewidth=params['lw'])
        ax.set(xlabel='time (s)', ylabel='Fraction inside cluster', title='Fraction of filaments inside cluster', ylim=[-0.05,1.05])
        fig.suptitle( params['ftitle'])
        fig.savefig( os.path.join( params['rpath'],'fractionInCluster.pdf') , bbox_inches='tight', dpi=600)
        plt.close('all')

    return cct2,ftlist_trimmed
# }}}

# calcNematicOrder {{{
@timer
def calcNematicOrder( ftlist, ftlist_trimmed, params):

    nT = params['nT']
    nF = params['nF']
    AC = params['AC']
    rpath = params['rpath']
    dt = params['dt']

    # bulk
    Sall_bulk = []
    for it in range(nT):
        Sall_bulk.append( AC.NematicOrder( ftlist[it]))

    fig,ax = plt.subplots(figsize=(6,6))
    ax.plot(dt*np.arange(nT), Sall_bulk, color='blue', linewidth=2,label='bulk')

    # cluster
    Sall = []
    for it in range(nT):
        Sall.append( AC.NematicOrder( ftlist_trimmed[it]))

    ax.plot(dt*np.arange(nT), Sall, linewidth=2, color='red', label='cluster')
    ax.set(xlabel='Time(s)', ylabel='Nematic Order S', title='Nematic Order evolution',ylim=[-0.1,1.1])
    plt.legend()
    fig.savefig( os.path.join( params['rpath'],'nematicOrder.pdf') , bbox_inches='tight', dpi=600)
    plt.close('all')

    return Sall, Sall_bulk

# }}}

# calcLocalPolarOrder {{{
@timer
def calcLocalPolarOrder( ftlist, cct, params):

    nT = params['nT']
    nF = params['nF']
    AC = params['AC']
    rpath = params['rpath']
    dt = params['dt']
    boxsize = params['boxsize']

    dist_max = 4*0.025

    # filament centers
    # c = np.zeros( (nT,nF,3))
    # for it in range(nT):
    #     for jf in range(nF):
    #         c[it,jf,:] = getMean( ftlist[it][jf].pos0, ftlist[it][jf].pos1, params['boxsize'])

    # filament plus-ends
    c = np.zeros( (nT,nF,3))
    for it in range(nT):
        for jf in range(nF):
            c[it,jf,:] = ftlist[it][jf].pos1

    # bulk
    Porder_par = np.zeros((2,nT))
    Porder_apar = np.zeros((2,nT))
    bool_all = np.array(1+np.arange(nF), dtype=bool)
    for it in range(nT):
        dist_within = within_range(c[it,:,:],nF,dist_max, boxsize,bool_all)
        porder_par = []
        porder_apar = []
        for ifil in range(nF):
            # fils_neatest = ftlist[it][ np.argwhere( dist_within[ifil,:]).transpose()[0].tolist()]
            fils_nearest = list( compress( ftlist[it], dist_within[ifil,:].tolist()))
            if not fils_nearest:
                continue
            orient_ref = ftlist[it][ifil].GetOrientation( boxsize)
            orientList_par = []
            orientList_apar = []
            for fil in fils_nearest:
                orient_new = fil.GetOrientation( boxsize)
                if orient_ref.dot( orient_new) / (1*1) >= 0:
                    orientList_par.append( orient_new)
                else:
                    orientList_apar.append( orient_new)
            PList_par = np.array(orientList_par)
            PList_apar = np.array(orientList_apar)
            if PList_par.any():
                porder_par.append(np.linalg.norm( np.mean(PList_par, axis=0)) )
            if PList_apar.any():
                porder_apar.append(np.linalg.norm( np.mean(PList_apar, axis=0)) )
        if porder_par:
            Porder_par[:,it] = [np.mean( porder_par), np.std(porder_par)]
        if porder_apar:
            Porder_apar[:,it] = [np.mean( porder_apar), np.std(porder_apar)]

    # cluster
    Porder_par_cluster = np.zeros((2,nT))
    Porder_apar_cluster = np.zeros((2,nT))
    for it in range(nT):
        dist_within = within_range(c[it,:,:],nF,dist_max, boxsize,cct[it])
        porder_par = []
        porder_apar = []
        for ifil in range(nF):
            fils_nearest = list( compress( ftlist[it], dist_within[ifil,:].tolist()))
            if not fils_nearest:
                continue
            orient_ref = ftlist[it][ifil].GetOrientation( boxsize)
            orientList_par = []
            orientList_apar = []
            for fil in fils_nearest:
                orient_new = fil.GetOrientation( boxsize)
                if orient_ref.dot( orient_new) / (1*1) >= 0:
                    orientList_par.append( orient_new)
                else:
                    orientList_apar.append( orient_new)
            PList_par = np.array(orientList_par)
            PList_apar = np.array(orientList_apar)
            if PList_par.any():
                porder_par.append(np.linalg.norm( np.mean(PList_par, axis=0)) )
            if PList_apar.any():
                porder_apar.append(np.linalg.norm( np.mean(PList_apar, axis=0)) )
        if porder_par:
            Porder_par_cluster[:,it] = [np.mean( porder_par), np.std(porder_par)]
        if porder_apar:
            Porder_apar_cluster[:,it] = [np.mean( porder_apar), np.std(porder_apar)]

    return Porder_par, Porder_apar, Porder_par_cluster, Porder_apar_cluster

# }}}

@jit( nopython=True)
def within_range(c,nF,dist_max,boxsize,boolRodsInside):
    # Calculate mean pair-pair filament distances for all filaments
    p2p = np.zeros((nF,nF))
    for jf in range(nF):
        if not boolRodsInside[jf]:
            continue
        dist = np.absolute( c[jf,:]-c)
        # get distance to all other filaments
        for jf2 in range(nF):
            if not boolRodsInside[jf2]:
                continue
            elif jf==jf2:
               p2p[jf,jf2] = 100*dist_max
               p2p[jf2,jf] = 100*dist_max
            elif p2p[jf2,jf] == 0:
                res2=0
                # get distance
                for idx in range(3):
                    k = np.floor( dist[jf2,idx]/(0.5*boxsize[idx]))
                    dist[jf2,idx] -= k*boxsize[idx]
                    res2 += dist[jf2,idx]**2
                p2p[jf,jf2] = np.sqrt(res2)
            else:
                p2p[jf,jf2] =p2p[jf2,jf]
    return (p2p < dist_max)

if __name__ == "__main__":

    paths = [
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf01_d10/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf01_d15/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf01_d20/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf01_d25/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf01_d30/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf02_d10/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf02_d15/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf02_d20/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf02_d25/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf02_d30/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf04_d10/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf04_d15/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf04_d20/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf04_d25/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf04_d30/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf08_d10/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf08_d15/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf08_d20/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf08_d25/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf08_d30/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf16_d10/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf16_d15/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf16_d20/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf16_d25/s0",
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf_const_num/run/pf16_d30/s0",
            ]
    analyzeSims( paths)













