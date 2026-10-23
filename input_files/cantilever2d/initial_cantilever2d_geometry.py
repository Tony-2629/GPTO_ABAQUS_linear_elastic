# def initial_cantilever2d_geometry():
## Initial design input file 
#
# *** THIS SCRIPT HAS TO BE CUSTOMIZED BY THE USER ***
#
# In this file, you must create two matrices that describe the initial
# design of bars.
#
# The first matrix contains the IDs (integer) and coordinates of the
# endpoints of the bars (point_matrix).
#
# The second matrix defines the IDs of the points that make up each bar.
# This matrix also sets the initial value of each bar's size variable, and
# the initial bar radius (half-width of the bar in 2-d).
#
# Note that this way of defining the bars allows for bars to be 'floating'
# (if the endpoints of a bar are not shared by any other bar) or
# 'connected' (if two or more bars share the same endpoint).
#

# *** Do not modify the line below ***
import numpy as np
# global FE, GEOM 

# Format of point_matrix is [ point_id, x, y] for 2-d problems, and 
# [ point_id, x, y, z] for 3-d problems)

point_matrix = np.array( (
    ( 1 , 2.4, 2.5 ) ,
    ( 2 , 2.6, 2.5 ) ,
    ( 3 , 7.4, 2.5 ) ,
    ( 4 , 7.6, 2.5 ) ,
    ( 5 , 12.4, 2.5 ) ,
    ( 6 , 12.6, 2.5 ) ,
    ( 7 , 17.4, 2.5 ) ,
    ( 8 , 17.6, 2.5 ) ,
    ( 9 ,  2.4, 7.5 ) ,
    ( 10 ,  2.6, 7.5 ) ,
    ( 11 ,  7.4, 7.5 ) ,
    ( 12 ,  7.6, 7.5 ) ,
    ( 13 , 12.4, 7.5 ),
    ( 14 , 12.6, 7.5 ),
    ( 15 , 17.4, 7.5 ) ,
    ( 16 , 17.6, 7.5 ) ) )

# Format of bar_matrix is [ bar_id, pt1, pt2, alpha, w/2 ], where alpha is
# the initial value of the bar's size variable, and w/2 the initial radius
# of the bar.
#
bar_matrix = np.array((
    ( 1 , 1 , 2 , 0.5, 0.5 ) ,
    ( 2 , 3 , 4 , 0.5, 0.5 ) ,
    ( 3 , 5 , 6 , 0.5, 0.5 ) ,
    ( 4 , 7 , 8 , 0.5, 0.5 ) ,
    ( 5 , 9 ,10 , 0.5, 0.5 ) ,
    ( 6 ,11 ,12 , 0.5, 0.5 ) ,
    ( 7 ,13 ,14 , 0.5, 0.5 ) ,
    ( 8 ,15, 16 , 0.5, 0.5 ) ) )

# *** Do not modify the code below ***
GEOM['initial_design']['point_matrix'] = point_matrix
GEOM['initial_design']['bar_matrix'] = bar_matrix

NOP = GEOM['initial_design']['point_matrix'].shape[0]
NOB = GEOM['initial_design']['bar_matrix'].shape[0]

GEOM['n_point'] = GEOM['initial_design']['point_matrix'].shape[0]
GEOM['n_bar'] = GEOM['initial_design']['bar_matrix'].shape[0]

print('initialized ' + str(FE['dim']) + 'd initial design with ' + 
    str( GEOM['initial_design']['point_matrix'].shape[0] ) + ' points and ' +
    str( GEOM['initial_design']['bar_matrix'].shape[0] ) + ' bars\n' )

with open('ZDATA.txt', 'w') as f:
    f.write('TEST INPUT DATA OF BARS\n')
    f.write('{:5d}{:5d}\n'.format(GEOM['initial_design']['point_matrix'].shape[0],
        GEOM['initial_design']['bar_matrix'].shape[0]))
    for i in range(NOP): 
        f.write('{:20.15f}{:20.15f}{:20.15f}\n'.format(point_matrix[i][0],
            point_matrix[i][1], point_matrix[i][2]))
    for i in range(NOB): 
        f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(bar_matrix[i][0],
            bar_matrix[i][1], bar_matrix[i][2],bar_matrix[i][3],bar_matrix[i][4]))        