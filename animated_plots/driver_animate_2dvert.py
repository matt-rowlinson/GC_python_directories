#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=animate_gc
#SBATCH --ntasks=1
#SBATCH --mem=2gb
#SBATCH --time=01:00:00
#SBATCH --output=LOGS/vert_%a_%A.log
import os
import sys 
sys.path.append('/users/mjr583/python_lib')
import GC_tools as GC
import imageio
import glob

def submit_batch_vertical(rundir, variable, nt, version='12.9.3'):
    print('%s plotting jobs to submit' %nt)
    outname='/users/mjr583/scratch/GC/%s/%s/plots/%s_vertical.mp4' % (version, rundir, variable)
    os.system("sbatch --wait  --array=1-%s anim_vertical.py -r %s -v %s -V %s" % (nt, rundir, variable, version))

    ## animate and delete .pngs 
    with imageio.get_writer(outname, mode='I') as writer:
        for png in sorted(glob.glob('/users/mjr583/scratch/GC/%s/%s/plots/vertical_*%s*png' % (version, rundir, variable))):
            print(png)
            image = imageio.imread(png) 
            writer.append_data(image)
    os.system("rm /users/mjr583/scratch/GC/%s/%s/plots/vertical_*%s*png" % (version, rundir, variable) )

def main():
    inputs=GC.get_arguments()
    rundir=inputs.rundir
    variable=inputs.var
    version=inputs=GC.get_arguments()

    ## get number of arrays to plot (timestep)
    nt=GC.get_n_timesteps(rundir, version)

    submit_batch_plot(rundir, variable, nt, version)
    
if __name__ == "__main__":
    main()
