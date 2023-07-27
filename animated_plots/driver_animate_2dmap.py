#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=animate_gc
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH --output=LOGS/anim_%A.log
'''
Script to create animated plot of surface concentrations from Geos-Chem. Submits jobs to slurm for each timestep
the creates an mp4 file from the individual plots. Ends by deleting the individual png files.

Example use:
    " sbatch driver_animate_2dmap.py -r path/to/dir/ -v "GC variable to plot" 
Author:
-------
Matt Rowlinson
email@institute.ac.uk
XX Jan 20XX
'''
import os
import sys 
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import imageio
import glob

def submit_batch_plot(rundir, variable, nt, version='12.9.3', plot_pressure=False, plot_strmfunc=False):
    print('%s plotting jobs to submit' %nt)
    if plot_pressure:
        outname='/users/mjr583/scratch/GC/%s/rundirs/%s/plots/%s_%s_ps.mp4' % (version, rundir, variable, rundir)
        os.system("sbatch --wait  --array=1-%s anim_with_p.py -r %s -v %s -V %s" % (nt, rundir, variable, version))
    elif plot_strmfunc:
        outname='/users/mjr583/scratch/GC/%s/rundirs/%s/plots/%s_%s_strmfunc.mp4' % (version, rundir, variable, rundir)
        os.system("sbatch --wait  --array=1-%s anim_with_strmfunc.py -r %s -v %s -V %s" % (nt, rundir, variable, version))
    else:
        outname='/users/mjr583/scratch/GC/%s/rundirs/%s/plots/%s_%s.mp4' % (version, rundir, variable, rundir)
        os.system("sbatch --wait  --array=1-%s anim_map.py -r %s -v %s -V %s" % (nt, rundir, variable, version))
    
    ## animate pngs then delete individual plots
    with imageio.get_writer(outname, mode='I') as writer:
        for png in sorted(glob.glob('/users/mjr583/scratch/GC/%s/rundirs/%s/plots/pcolorm_*%s*png' % (version, rundir, variable))):
            print(png)
            image = imageio.imread(png) 
            writer.append_data(image)
    os.system("rm /users/mjr583/scratch/GC/%s/rundirs/%s/plots/pcolorm_*%s*png" % (version, rundir, variable) )
    return


def main():
    inputs=GC.get_arguments()
    rundir=inputs.rundir
    variable=inputs.var
    version=inputs.version
    plot_pressure=inputs.plot_ps
    plot_strmfunc=inputs.strmfunc
    
    ## get number of arrays to plot (timestep)
    nt=GC.get_n_timesteps(rundir, version)
    
    submit_batch_plot(rundir, variable, nt, version, plot_pressure, plot_strmfunc)
    
if __name__ == "__main__":
    main()
