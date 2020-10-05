from paraview.simple import *
import os
import glob
import subprocess

paraview.simple._DisableFirstRenderCameraReset()

maindir = '/Users/saadjansari/Documents/Projects/AMSOS/resultsSummit/Confinement/scan_d_pf'
os.chdir(maindir)
state_file = os.path.join(maindir, 'display_example.pvsm')

# find all sims to analyze
sims = glob.glob( os.path.join(maindir, '*/*/*/Result2PVD.py'))
fnames = [
        "pf02_d20",
        "pf08_d15",
        "pf08_d25",
        ]

for simdir in sims:

    simdir = '/'.join( simdir.split('/')[:-1])

    sname = os.path.basename('/'.join( simdir.split('/')[:-2]) )
    n_syl = os.path.join( simdir, 'Sylinderpvtp.pvd')
    n_pro = os.path.join( simdir, 'Proteinpvtp.pvd')
    sbox = os.path.join( simdir, 'simBox.vtk')
    animationpath = os.path.join(simdir, 'PNG/frame.png')

    print('Loading State for sim: {}'.format(sname))
    if sname not in fnames:
        print("skipped {0}".format(sname))
        continue

    LoadState( state_file, LoadStateDataFileOptions='Search files under specified directory',
          DataDirectory=simdir, 
          OnlyUseFilesInDataDirectory=1,
          simBoxvtkFileNames=sbox,
          SylinderpvtppvdFileName=n_syl,
          ProteinpvtppvdFileName=n_pro)

    LoadPalette("WhiteBackground")

    # get animation scene
    animationScene1 = GetAnimationScene()

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')

    # get layout
    layout1_1 = GetLayoutByName("Layout #1")

    # set active view
    SetActiveView(renderView1)

    # Enter preview mode
    layout1_1.PreviewMode = [1920, 1280]

    savesettings = {'SaveAllViews': 1, ''
                    'viewOrLayout': layout1_1,
                    'SeparatorWidth': 5,
                    'SuffixFormat': "_%06d"
                    }

    print('Saving animation...')
    SaveAnimation(animationpath, **savesettings)
    os.chdir( os.path.join(simdir, "PNG") )
    subprocess.call(["bash","MovieGen.sh"])
