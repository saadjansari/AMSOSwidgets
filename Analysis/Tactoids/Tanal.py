import os, time, yaml, glob, copy
import math, random, shutil, random
import numpy as np
import matplotlib.pyplot as plt
import pickle

from pytransform3d import rotations as pyro
from mpl_toolkits.mplot3d import Axes3D

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
    cct,ftlist_trimmed,fdat = getFilamentsInsideCluster( ftlist, ptlist, params, plotting=1)
    # ftlist2,_ = transformCM( ftlist, ptlist, ftlist_trimmed, params)

    # filament centers 
    # c = np.zeros( (nT,nF,3))
    # for it in range(nT):
        # for jf in range(nF):
            # c[it,jf,:] = getMean( ftlist[it][jf].pos0, ftlist[it][jf].pos1, params['boxsize'])

    #Pair-pair separation vs time
    # PPseparation(c,cct, params)

    # MSD
    # lags = np.arange(1,nT)
    # MSDpairs(lags, c, cct, params)

    # Nematic Order
    sdat = calcNematicOrder( ftlist, ftlist_trimmed, params)

    # Xlink analysis
    num_prot, ratio_prot = analyzeXlinks( ptlist, ftlist_trimmed, cct, params)

    dat = {
            'S_bulk': sdat['bulk'],
            'S_cluster': sdat['cluster'],
            'xlink_num': num_prot,
            'xlink_ratio': ratio_prot,
            'actin_num': fdat,
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

    # Specify parameters
    analyze_bulk = False

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

    FAratio = np.around( nprot/nfil,1) # filamin to actin ratio

    # Determine if filaments are polydisperse or monodisperse
    with open( os.path.join(tpath, 'result/result0-399/SylinderAscii_0.dat'), 'r') as file1:
        lens = []
        for line in file1:
            if line.startswith('C'):
                data = line.split()
                pos0 = np.array([float(data[3]), float(data[4]), float(data[5])])
                pos1 = np.array([float(data[6]), float(data[7]), float(data[8])])
                lens.append( np.linalg.norm(pos1-pos0))
    if max(lens) - min(lens) < 0.001:
        dispersity = 'Mono'
    else:
        dispersity = 'Poly'
        
    # title of plots
    ftitle='Filamin-Actin ratio = {0}, {1}disperse'.format(FAratio,dispersity)
    print('   {}'.format(ftitle))

    # Display
    cols = {
        'cluster': "limegreen",
        'bulk': "mediumvioletred",
        'env': "dodgerblue"
    }
    lw = 2 # linewidth

    AC = AnalysisBook(boxsize)

    params = {
            'name': os.path.basename(tpath),
            'tpath':tpath,
            'rpath':rpath,
            'analyze_bulk':analyze_bulk,
            'dt':dt,
            'boxsize':boxsize,
            'ftitle':ftitle,
            'cols':cols,
            'lw':lw,
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

    if not params['analyze_bulk']:
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
            # cct.append( cc[0])
    else:
        cct2 = np.ones( (nT,nF),dtype=bool)
        # cct = [range(nF) for ii in range(nT)]

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

    return cct2,ftlist_trimmed, nft
# }}}

# transformCM {{{
@timer
def transformCM( ptlist, ftlist, ftlist_trimmed, params):

    nT = params['nT']
    nF = params['nF']
    AC = params['AC']

    ftlist2 = copy.deepcopy( ftlist)
    ptlist2 = copy.deepcopy( ptlist)

    world_origin = np.array([0,0,0])
    world_vec = np.identity(3)

    for it in range(nT):
        obj_origin = AC.COM( ftlist_trimmed[it]).transpose()
        w,obj_vec = AC.GetBasisVectors( ftlist_trimmed[it])
        T = (world_origin - obj_origin)

        # Get vector perp to z-axis and nematic director
        vperp = pyro.perpendicular_to_vectors( world_vec[:,-1], obj_vec[:,-1])
        if np.linalg.norm( vperp) >0.00001:
            # rotation needed
            vperp = vperp / np.linalg.norm(vperp)
            theta = pyro.angle_between_vectors( world_vec[:,-1], obj_vec[:,-1])
            axis_angle = np.concatenate( (vperp, -np.array([theta]) ))
            q = pyro.quaternion_from_axis_angle( axis_angle)

            # transfrom filaments
            for fil in ftlist2[it]:
                fil.pos0 = pyro.q_prod_vector( q, fil.pos0 +T.transpose() )
                fil.pos1 = pyro.q_prod_vector( q, fil.pos1 +T.transpose() )
            # transfrom proteins
#        for prot in ptlist2[it]:
#            prot.pos0 = pyro.q_prod_vector( q, prot.pos0 +T.transpose() )
#            prot.pos1 = pyro.q_prod_vector( q, prot.pos0 +T.transpose() )

        else:
            # transfrom filaments
            for fil in ftlist2[it]:
                fil.pos0 = fil.pos0 +T.transpose()
                fil.pos1 = fil.pos1 +T.transpose()
            # transfrom proteins
#        for prot in ptlist2[it]:
#            prot.pos0 = prot.pos0 +T.transpose()
#            prot.pos1 = prot.pos0 +T.transpose()

    return ftlist2, ptlist2
# }}}

# PPseparation {{{
@timer
def PPseparation(c, cct, params, plotting=1):
    
    rho_mu, rho_sig = meanPairSep( c, params['nT'], params['nF'], params['boxsize'], cct)

    if plotting:
        fig,ax = plt.subplots()
        ax.plot( params['dt']*np.arange(params['nT']), rho_mu, color=params['cols']['cluster'], label='cluster', linewidth=params['lw'])
        ax.fill_between( params['dt']*np.arange(params['nT']), rho_mu-rho_sig, rho_mu+rho_sig, color=params['cols']['cluster'], alpha=0.1)
        ax.set(xlabel='time (s)', ylabel='separation ($\mu$m)', title='Pair-pair separation (cluster)')
        fig.suptitle( params['ftitle'])
        fig.savefig( os.path.join( params['rpath'],'mean_separation_cluster.pdf') , bbox_inches='tight', dpi=600)
        plt.close('all')
# }}}

# MSDpairs {{{
@timer
def MSDpairs( lags, c, cct, params, plotting=1):

    lags = np.arange(params['nT'])
    rho_mu, rho_sig = msd_pairs( lags, c, params['nT'],params['nF'],params['boxsize'], cct)

    if plotting:
        fig,ax = plt.subplots()
        ax.plot( params['dt']*np.arange(params['nT']), rho_mu, color=params['cols']['cluster'], label='cluster', linewidth=params['lw'])
        ax.fill_between( params['dt']*np.arange(params['nT']), rho_mu-rho_sig, rho_mu+rho_sig, color=params['cols']['cluster'], alpha=0.1)
        ax.set(xlabel='time (s)', ylabel='MSD ($\mu m^2$)', title='MSD Pairs (cluster)')
        fig.suptitle( params['ftitle'])
        fig.savefig( os.path.join( params['rpath'],'msd_pairs.pdf') , bbox_inches='tight', dpi=600)
        plt.close('all')
# }}}

# calcNematicOrder {{{
@timer
def calcNematicOrder( ftlist, ftlist_trimmed, params):

    nT = params['nT']
    nF = params['nF']
    AC = params['AC']
    cols = params['cols']
    lw = params['lw']
    rpath = params['rpath']
    ftitle = params['ftitle']
    dt = params['dt']

    # bulk
    Sall_bulk = []
    for it in range(nT):
        Sall_bulk.append( AC.NematicOrder( ftlist[it]))
        
    fig,ax = plt.subplots(figsize=(6,6))
    ax.plot(dt*np.arange(nT), Sall_bulk, color=params['cols']['bulk'], linewidth=params['lw'],label='bulk')

    # cluster
    Sall = []
    for it in range(nT):
        Sall.append( AC.NematicOrder( ftlist_trimmed[it]))

    ax.plot(dt*np.arange(nT), Sall, color=params['cols']['cluster'], linewidth=params['lw'],label='cluster')
    ax.set(xlabel='Time(s)', ylabel='Nematic Order S', title='Nematic Order evolution',ylim=[-0.1,1.1])
    fig.legend()
    plt.suptitle(params['ftitle'])
    fig.savefig( os.path.join( params['rpath'],'nematicOrder.pdf') , bbox_inches='tight', dpi=600)
    plt.close('all')

    dat = {
            'bulk': Sall_bulk,
            'cluster': Sall
            }

    return dat

    # Plot theta evolution with time
    # theta = np.zeros(nT)
    # for it in range(nT):
        # world_vec = np.identity(3)
        # w,obj_vec = AC.GetBasisVectors( ftlist_trimmed[it])    
        # theta[it] = pyro.angle_between_vectors( world_vec[:,-1], obj_vec[:,-1])

    # fig,ax = plt.subplots(figsize=(6,6))
    # ax.plot(dt*np.arange(nT), theta, color=cols['cluster'], linewidth=lw,label='cluster')
    # ax.set(xlabel='Time(s)', ylabel='theta (rad)', title='Angle b/w Z and nematic eigenvector')
    # plt.suptitle(ftitle)

# }}}

# analyzeXlinks {{{
@timer
def analyzeXlinks( ptlist, ftlist_trimmed, cct, params):

    nT = params['nT']
    nF = params['nF']
    AC = params['AC']
    cols = params['cols']
    lw = params['lw']
    rpath = params['rpath']
    ftitle = params['ftitle']
    dt = params['dt']
    boxsize = params['boxsize']

    # Energy
    k=50
    restL = 0.05
    xlinkE = np.zeros((nT))
    xlinkE_err = np.zeros((nT))
    for iT in range(nT):
        vals = [ 0.5*k*(xlink.GetLength(boxsize)-restL)**2 for xlink in ptlist[iT]]
        xlinkE[iT] = np.mean( vals)
        xlinkE_err[iT] = np.std(vals)
    fig,ax = plt.subplots(figsize=(6,6))
    ax.plot(dt*np.arange(nT), xlinkE, color=cols['cluster'], linewidth=lw)
    ax.fill_between(dt*np.arange(nT), xlinkE - xlinkE_err, xlinkE + xlinkE_err, color=cols['cluster'], alpha=0.1)
    ax.set(xlabel='time (s)', ylabel='mean energy ($pN. \mu m$)', title='Xlink energy vs time')
    fig.suptitle(ftitle)
    fig.savefig( os.path.join(rpath,'xlinkenergy.pdf') , bbox_inches='tight', dpi=600)
    plt.close('all')

    # Distance
    # restL = 0.05
    # xlinkD = np.zeros((nT))
    # xlinkD_err = np.zeros((nT))
    # for iT in range(nT):
        # com = AC.COM( ftlist_trimmed[iT])
        # vals = [ np.linalg.norm( getDistance( com, xlink.GetCenter(boxsize), boxsize)) for xlink in ptlist[iT]]
        # xlinkD[iT] = np.mean( vals)
        # xlinkD_err[iT] = np.std( vals)
    # fig,ax = plt.subplots(figsize=(6,6))
    # ax.plot(dt*np.arange(nT), xlinkD, color=cols['cluster'], linewidth=lw)
    # ax.fill_between(dt*np.arange(nT), xlinkD - xlinkD_err, xlinkD + xlinkD_err, color=cols['cluster'], alpha=0.1)
    # ax.set(xlabel='time (s)', ylabel='mean distance ($\mu m$)', title='Xlink distance from CM')
    # fig.suptitle(ftitle)
    # fig.savefig( os.path.join(rpath,'xlinkdistance_cm.pdf') , bbox_inches='tight', dpi=600)
    # plt.close('all')

    # num xlinks in biggest cluster
    nprot = np.zeros(nT)
    ratio_prot = np.zeros(nT)
    for it in range(nT):
        idxRods = np.where(cct[it]==True)[0]
        for prot in ptlist[it]:
            if prot.link0 in idxRods or prot.link1 in idxRods:
                nprot[it]+=1
        ratio_prot[it] = nprot[it]/np.sum(cct[it])

    fig,ax = plt.subplots(2,1,figsize=(6,6))
    ax[0].plot(dt*np.arange(nT), ratio_prot, color=cols['cluster'], linewidth=lw)
    ax[0].set(xlabel='time (s)', ylabel='ratio', title='Ratio of crosslinks to filaments inside cluster')
    ax[1].plot(dt*np.arange(nT), nprot, color=cols['cluster'], linewidth=lw)
    ax[1].set(xlabel='time (s)', ylabel='count', title='Num filamin inside cluster')
    fig.suptitle(ftitle)
    fig.savefig( os.path.join(rpath,'xlinkInCluster.pdf') , bbox_inches='tight', dpi=600)
    plt.close('all')

    return nprot, ratio_prot

# }}}

if __name__ == "__main__":

    paths = [
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin1.0/s0", 
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin1.5/s0", 
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin2.0/s0", 
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin2.5/s0", 
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin5/s0", 
            "/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Tactoids/scan_filamin/run/filamin10/s0", 
            ]
    analyzeSims( paths)













