from paraview.simple import *
import os
import argparse
import glob
import subprocess
import pdb

'''
Name: paraview_snapshots.py
Description: Loads .pvd files in paraview with special filters and saves animation frames
Input: To see type paraview_snapshots.py -h
'''

def parseArgs():

    parser = argparse.ArgumentParser(prog='paraview_snapshots.py')
    parser.add_argument('-T', '--tactoid', action='store_true', default=False, help='use tactoid options for loading')
    parser.add_argument('-C', '--confine', action='store_true', default=False, help='use confinement options for loading')
    parser.add_argument('--all', action='store_true', default=False, help='save all frames. default saves just the last frame')
    opts = parser.parse_args()

    if opts.tactoid and opts.confine:
        raise Exception('Cannot specify both tactoid and confinement opts')
    elif not opts.tactoid and not opts.confine:
        raise Exception('Please specify either tactoid or confinement opts')
    if opts.tactoid:
        print('Configuration: Tactoid')
        opts.fil_radius = 0.0035
        opts.fil_colorby = ('CELLS', 'length')
    elif opts.confine:
        print('Configuration: Confine')
        opts.fil_radius = 0.0125
        opts.fil_colorby = ('POINTS', 'endLabel')

    return opts


def loadPVD(simdir, opts):
    
    os.chdir(simdir)

    # paraview.simple._DisableFirstRenderCameraReset()

    # sname = os.path.basename('/'.join( simdir.split('/')[:-2]) )
    f_syl = os.path.join( simdir, 'result/Sylinderpvtp.pvd')
    f_pro = os.path.join( simdir, 'result/Proteinpvtp.pvd')
    f_box = os.path.join( simdir, 'result/simBox.vtk')
    screenshot_path = os.path.join(simdir, 'result/last_frame.png')

# Load sylinder file
    reader_sylinder = OpenDataFile( f_syl)
    Show(reader_sylinder)
    layout1 = GetLayout()
    pdb.set_trace()

# Apply tube filter to rods
    fils = Tube(Input=reader_sylinder)
    fils.Radius = opts.fil_radius
    rep_syl = Show(fils)
    Hide(reader_sylinder)
    ColorBy( rep_syl, opts.fil_colorby)

# Load proteinfile
    reader_protein = OpenDataFile( f_pro)
    rep_pro = Show(reader_protein)
    rep_pro.DiffuseColor = [0,1,0]

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
    animationScene1.GoToLast()
    layouta = GetLayoutByName("Layout #1")
    layoutb = GetLayout()

    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
        # renderView1.ViewSize = [1280, 1280]]
# get layout
    # layout1_1 = GetLayoutByName("Layout #1")
    layout1_1 = GetLayout()
    pdb.set_trace()

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
    # Render()

# Go to last timestep
    SaveScreenshot(screenshot_path, layout=layout1_1)


if __name__ == "__main__":

    simdir = '/Users/saadjansari/Documents/Projects/Results/AMSOS/Confinement/scan_d_pf_const_num/run/pf01_d125/s0'
    opts = parseArgs()
    loadPVD( simdir, opts)
