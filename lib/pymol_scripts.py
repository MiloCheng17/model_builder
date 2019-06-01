#!/usr/bin/env python
"""
This is a program written by Qianyi Cheng in DeYonker Research Group
at University of Memphis.
"""

###  launches PyMol for stereo viewing on a VisBox. It runs PyMol fullscreen stereo, and disables the internal gui. The environment (PYTHON_PATH and PYMOL_PATH) must already be set up

### Tell pymol to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
##import __main__
##__main__.pymol_argv = [ 'pymol', '-qei' ]
##
##import pymol
### Call the function below before using any PyMOL modules.
##pymol.finish_launching()
##
##from pymol import cmd
##cmd.stereo('walleye')
##cmd.set('stereo_shift', 0.23)
##cmd.set('stereo_angle', 1.0cmd.set('stereo_angle', 1.0)

### launches PyMOL without any GUI for scripting only

import os, sys

# autocompletion
#import readline
#import rlcompleter
#readline.parse_and_bind('tab: complete')
#
## pymol environment
#moddir='/opt/pymol-svn/modules'
#sys.path.insert(0, moddir)
#os.environ['PYMOL_PATH'] = os.path.join(moddir, 'pymol/pymol_path')
#
## pymol launching: quiet (-q), without GUI (-c) and with arguments from command line
#import pymol
#pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
#pymol.finish_launching()
#cmd = pymol.cmd

#### Example: pymol_script.py res_3.pdb 300,301,302 ####

def system_run(cmd):
    print cmd
    exit = os.system(cmd)
    if ( exit != 0 ):
        print 'failed to run:'
        print cmd
        sys.exit()

input = sys.argv[1]
name  = input.split('.')[0]
ouput = name+'_h.pdb'
logf = open('log.pml','w')
if len(sys.argv) == 3:
    logf.write('load %s\ncmd.select("sel","%s and not resi %s")\ncmd.h_add("sel")\ncmd.save("./%s")'%(input,name,sys.argv[2],ouput))
elif len(sys.argv) == 2:
    logf.write('load %s\ncmd.h_add("%s")\ncmd.save("./%s")'%(input,name,ouput))
else:
    print "Something is wrong!"
print "Please run 'pymol -qc log.pml'"
#system_run('pymol -qc log.pml')

## preset.ball_and_stick(selection='all', mode=1) or preset.ball_and_stick(selection='all', mode=2)
## mode=1
## set_bond stick_color, white, selection, selection
## set_bond stick_radius, 0.14, selection, selection
## set sphere_scale, 0.25, selection
## show sticks, selection
## show spheres, selection

## mode=2
## set_bond stick_color, white, selection, selection
## set_bond stick_radius, -0.14, selection, selection
## set stick_ball, 1
## set stick_ball_ratio, -1.0
## set stick_ball_color, atomic
## show sticks, selection
