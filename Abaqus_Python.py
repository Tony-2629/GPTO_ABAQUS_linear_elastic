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
from geometry_projection import *
from optimization import *
from functions import *
from utilities import *
from plotting import *

""" Batch Settings """
#nCores = 2
#abaqusVersion = "abaqus"
#inputFile = 'Opt_MMA'
#userFile = 'UMAT_SUBROUTINE'



def process(FE,OPT,GEOM):
    print('RUN ABAQUS')
    command = "%s job=%s user=%s interactive" % (GEOM['abaqusVersion'], GEOM['inputFile']
    	,GEOM['userFile'])
    startTime = time.perf_counter()
    subprocess.call(command, shell=True)
	#
    while(os.path.exists(GEOM['inputFile']+'.lck')):
    	time.sleep(10)
	#
    endTime = time.perf_counter()		 
    JobTime = startTime-endTime
    report = 'FINISH SIMULATION RAW DATA in %d seconds' % JobTime
    print(report) 
    return True

def export_output_raw(FE,OPT,GEOM):
	abaqus_script = 'ZScriptToReadResults.py'
	print('EXPORTING RAW DATA')
	command = "%s cae noGUI=%s" % (GEOM['abaqusVersion'], abaqus_script)
	os.system(command)
	print('FINISH EXPORTING RAW DATA')

def convert_output_raw2csv():
	current_dir = os.getcwd()
	csv_data_path = os.path.join(current_dir,'csv')
	if os.path.isdir(csv_data_path):
		print('CONTINUE')
	else:
		os.makedirs(csv_data_path)
	# remove unnecessary information in text result data of abaqus file
	txt_file_path = os.path.join(current_dir,'Zcompliance_sensitivities.txt')
	txt_file = open(txt_file_path, 'r')
	lines = txt_file.readlines()
	txt_file.close()
	new_txt_file_path = os.path.join(current_dir,'Zcompliance_sensitivities_filter.txt')
	new_txt_file = open(new_txt_file_path, 'w')
	for idx, line in enumerate(lines):
		if idx <19:
			continue
		else:
			new_txt_file.write(line)
	new_txt_file.close()
	# Create a blank frame of full data
	output_file_path = os.path.join(csv_data_path, 'compliance_sensitivities.csv')
	csv_file = pd.read_csv(new_txt_file_path, sep='\s+', header=None)
	first_column = csv_file.columns[0]
	# Delete first
	csv_file = csv_file.drop([first_column], axis=1)
	#csv_file = csv_file.dropna()
	csv_file.to_csv(output_file_path,index=None)

def assign_value():
	current_dir = os.getcwd()
	csv_data_path = os.path.join(current_dir,'csv')
	output_file_path = os.path.join(csv_data_path, 'compliance_sensitivities.csv')
	f = np.loadtxt(open(output_file_path,"rb"),delimiter=",")
	return f
def evaluate_relevant_functions(FE,OPT,GEOM):
    # Evaluate_relevant_functions() looks at OPT['functions'] and evaluates the
    # relevant functions for this problem based on the current OPT['dv']
    # global OPT
    n_dv = OPT['n_dv']
    # compute and assign the value of compliance and compliance sensitivities (analytical method)
    theta_0 = OPT['cost'][1,n_dv]
    grad_theta_0 = np.zeros((n_dv, 1))
    grad_theta_0 = OPT['cost'][1,0:n_dv]
    # This computes the volume fraction and its sensitivities
    # base on the last geometry projection
    volume_dens = OPT['cost'][1,(n_dv+100)]
    volume_sum = OPT['cost'][1,(n_dv+100+1)]
    volume_frac = volume_dens/volume_sum
    grad_volfrac = np.zeros((n_dv, 1))
    grad_volfrac = OPT['cost'][1,100:100+n_dv]/volume_sum
    #output
    #compliance and compliance sensitivities
    OPT['compliance'] = theta_0
    OPT['grad_compliance'] = grad_theta_0
    # The design sensitivity
    OPT['volume_fraction'] = volume_frac
    OPT['grad_volume_fraction'] = grad_volfrac
#
    OPT['functions']['f'][0]['value'] = OPT['compliance']
    OPT['functions']['f'][0]['grad']  = OPT['grad_compliance']
    OPT['functions']['f'][1]['value'] = OPT['volume_fraction']
    OPT['functions']['f'][1]['grad']  = OPT['grad_volume_fraction']
def performance_analysis(FE,OPT,GEOM):
	process(FE,OPT,GEOM)
	export_output_raw(FE,OPT,GEOM)
	convert_output_raw2csv()
	OPT['cost'] = assign_value()
	evaluate_relevant_functions(FE,OPT,GEOM)










'''
FE = {}
GEOM = {}
OPT = {}
exec(open('get_inputs.py').read())

process()
export_output_raw()
convert_output_raw2csv()
compute_cost = assign_value()

init_geometry(FE,OPT,GEOM)
init_optimization(FE,OPT,GEOM)
n_dv = OPT['n_dv']
# compute and assign the value of compliance and compliance sensitivities (analytical method)
theta_0 = compute_cost[1,n_dv]
grad_theta_0 = np.zeros((n_dv, 1))
grad_theta_0 = compute_cost[1,0:n_dv]
# This computes the volume fraction and its sensitivities
# base on the last geometry projection
volume_dens = compute_cost[1,(n_dv+100)]
volume_sum = compute_cost[1,(n_dv+100+1)]
volume_frac = volume_dens/volume_sum
grad_volfrac = np.zeros((n_dv, 1))
grad_volfrac = compute_cost[1,100:99+n_dv]/volume_sum
'''

#output
#compliance
#OPT['compliance'] = theta_0
#OPT['grad_compliance'] = grad_theta_0
# The design sensitivity
#OPT['volume_fraction'] = volume_frac
#OPT['grad_volume_fraction'] = grad_volfrac

#print(OPT['compliance'])
#print(OPT['grad_compliance'])
#print(OPT['volume_fraction'])
#print(OPT['grad_volume_fraction'])
'''
# compute and assign the value of compliance and compliance sensitivities (analytical method)



grad_theta_0 = np.zeros((n_dv, 1))

grad_theta_i = np.zeros((n_dv, 1))
theta_0 = compute_cost[1,n_dv]
grad_theta_0 = compute_cost[1,0:n_dv]

# compute and assign the value of compliance and compliance sensitivities (Finite Different Method)

fd_step = OPT['fd_step_size']

up_to_dv = n_dv

dv_0 = OPT['dv'].copy()
dv_i = OPT['dv'].copy()



max_error         = 0.0
max_rel_error     = 0.0
max_error_bar     = 0
max_rel_error_bar = 0

for i in range(0, up_to_dv):
# Preturb dv
	dv_i[i] = dv_0[i] + fd_step
	theta_i = obj(dv_i)
	grad_theta_i[i] = (theta_i - theta_0)/fd_step
	print(i, grad_theta_0[i], grad_theta_i[i],theta_0,theta_i)
	#print(i)
	#print(dv_0[i],dv_i[i])
#compute error relevant between analytical method and Finite Difference Method
	error = grad_theta_0[i] - grad_theta_i[i]

	if np.abs(error) > np.abs(max_error):
		max_error = error
		max_error_dv = i

	rel_error = error / theta_0
	if np.abs(rel_error) > np.abs(max_rel_error):
		max_rel_error = rel_error
		max_rel_error_dv = i

	dv_i = dv_0.copy()
#	print(GEOM['current_design']['point_matrix'])
#	print(GEOM['current_design']['bar_matrix'])
print('Finish simulation running')
#print(grad_theta_i)
print('Max. ABSOLUTE error is:')
print(max_error)
print('It occurs at:')
print('\tvariable:' + str(max_error_dv))

print('Max. RELATIVE error is:')
print(max_rel_error)
print('It occurs at:')
print('\tvariable:' + str(max_rel_error_dv))
plt.figure(3)
plt.plot(grad_theta_i, '+')
plt.plot(grad_theta_0, '.')
#plt.legend(fd,analytical)
plt.title('cost function')
plt.xlabel('design variable z')
plt.ylabel('dc/dz')
plt.show()
'''
