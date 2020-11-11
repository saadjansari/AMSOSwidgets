from paraview.simple import *
import os
import argparse
import glob
import subprocess

'''
Name: paraview_snapshots.py
Description: Loads .pvd files in paraview with special filters and saves animation frames
Input: To see type paraview_snapshots.py -h
'''

maindir = '/Users/saadjansari/Documents/Projects/Results/AMSOS/Confinement/scan_d_pf_const_num/run'
sims = glob.glob( os.path.join(maindir, '*') )
snames = [ii.split('/')[-1] for ii in sims]
sdirs = [ os.path.join( ii, 's0') for ii in sims]
save_names = ['/Users/saadjansari/Desktop/conf_snaps/'+ii+'.png' for ii in snames]
# simdir = '/Users/saadjansari/Documents/Projects/Results/AMSOS/Confinement/scan_d_pf_const_num/run/pf01_d125/s0'

def sort_sim_names(s):
    return len(s.split('/')[-1].split('_c'))


print('Configuration: Confine')

for sname,simpath,save_name in zip(snames, sims, save_names):

    simdirs = glob.glob( os.path.join(simpath,'*')) 
    if os.path.join(simpath, 'merge') in simdirs:
        simdir = os.path.join(simpath, 'merge') 
        os.chdir( simdir)
    else:
        simdirs = sorted( simdirs, key=sort_sim_names)
        simdir = simdirs[-1]
        os.chdir( simdir)

    print('Parsing {}'.format(sname) )
# load .pvd files
    fil_radius = 0.0125
    fil_colorby = ('POINTS', 'endLabel')

    paraview.simple._DisableFirstRenderCameraReset()
    f_syl = os.path.join( simdir, 'result/Sylinderpvtp.pvd')
    f_pro = os.path.join( simdir, 'result/Proteinpvtp.pvd')
    f_box = os.path.join( simdir, 'result/simBox.vtk')

# Load sylinder file
    reader_sylinder = OpenDataFile( f_syl)
    Show(reader_sylinder)
    layout1 = GetLayout()

# Apply tube filter to rods
    fils = Tube(Input=reader_sylinder)
    fils.Radius = fil_radius
    rep_syl = Show(fils)
    Hide(reader_sylinder)
    ColorBy( rep_syl, fil_colorby)

# Load proteinfile
    reader_protein = OpenDataFile( f_pro)
    rep_pro = Show(reader_protein)
    rep_pro.DiffuseColor = [0,1,0]
    Hide(reader_protein)

# Set View
    view = GetRenderView()
    view.CameraPosition = [0,0,0]
    view.CameraFocalPoint = [1,0,0]
    ResetCamera()
    camera = GetActiveCamera()
    camera.SetViewUp([0,0,1])
# Render()

# find view
    LoadPalette("DefaultBackground")
# get animation scene
    animationScene1 = GetAnimationScene()
# get the time-keeper
    timeKeeper1 = GetTimeKeeper()
# animationScene1.GoToLast()
# layouta = GetLayoutByName("Layout #1")

    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
# renderView1.ViewSize = [1280, 1280]]
# get layout
    layout1_1 = GetLayoutByName("Layout #1")
# layout1_1 = GetLayout()

    SetActiveView(renderView1)
# Render()
    layout1_1.PreviewMode = [2560, 1280]

# Second view with only proteins
    renderView2 = FindViewOrCreate('RenderView2', viewtype='RenderView')
    Show(fils)
    rep_pro = Show(reader_protein)
    rep_pro.DiffuseColor = [0,1,0]
    view = GetRenderView()
    view.CameraPosition = [0,0,0]
    view.CameraFocalPoint = [1,0,0]
    ResetCamera()
    camera = GetActiveCamera()
    camera.SetViewUp([0,0,1])
    Hide(fils)
    Render()

    animationScene1.StartTime = 0
    animationScene1.NumberOfFrames = int(animationScene1.EndTime+1)
    SaveAnimation(save_name, viewOrLayout=layout1_1, scene=animationScene1, 
            FrameWindow=(animationScene1.NumberOfFrames-2, animationScene1.NumberOfFrames))

    print('Disconnecting...')
    Disconnect()
    Connect()
    print('New Connection!')


# Go to last timestep
# SaveScreenshot(screenshot_path, layout=layout1_1)

