import os
import numpy as np
import json
import pandas as pd
import math
import matplotlib.pyplot as plt
import sys
import csv
#
import glob
import time
import socket
import subprocess
#
from Abaqus_Python import *
from geometry_projection import *
from optimization import *
from functions import *
from utilities import *
from plotting import *

exec(open('get_inputs.py').read())
## Start timer
tic = time.perf_counter()
## Initialization
init_geometry(FE,OPT,GEOM)
init_optimization(FE,OPT,GEOM)
## Analysis
performance_analysis(FE,OPT,GEOM)
if OPT['make_fd_check']:
    run_finite_difference_check(FE,OPT,GEOM)

## Optimization
OPT['history'] = runopt(FE,OPT,GEOM,OPT['dv'],obj,nonlcon)

# Report time
toc = time.perf_counter()
print( "Elapsed time: " + str( toc - tic ))

# hold graph
plt.ioff()

# ## Plot History
if True == OPT['options']['plot']:
    plot_history(2)
### Make Video
if True == OPT['options']['video']:
    make_video()




