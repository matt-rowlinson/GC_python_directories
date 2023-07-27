#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#SBATCH --job-name=timeseries
#SBATCH --ntasks=1
#SBATCH --mem=10gb
#SBATCH --partition=interactive
#SBATCH --time=00:10:00
#SBATCH --output=LOGS/timeseries.log
import numpy as np

from read_NMVOCs import v21, v20, v18, v21_sums, v20_sums, v18_sums
v21_voc = np.array(v21)
v20_voc = np.array(v20)
v18_voc = np.array(v18)

from read_ethane import v21, v20, v18, v21_sums, v20_sums, v18_sums
v21_ethane = np.array(v21)
v20_ethane = np.array(v20)
v18_ethane = np.array(v18)

np.save('temp_files/v21_voc', v21_voc)
np.save('temp_files/v20_voc', v20_voc)
np.save('temp_files/v18_voc', v18_voc)

np.save('temp_files/v21_ethane', v21_ethane)
np.save('temp_files/v20_ethane', v20_ethane)
np.save('temp_files/v18_ethane', v18_ethane)
