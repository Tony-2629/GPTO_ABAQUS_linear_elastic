import numpy as np
from geometry_projection import *
from Abaqus_Python import *
GEOM['nCores'] = 4
GEOM['abaqusVersion'] = "abq2023"
GEOM['inputFile'] = 'GPTO_3D'
GEOM['userFile'] = 'ZUMAT_1'

def compute_compliance(FE,OPT,GEOM):
    # This function computes the mean compliance and its sensitivities
    # based on the last finite element analysis
    # global FE, OPT

    # compute the compliance (Eq. (15))
    c = np.dot( FE['U'].flatten() , FE['P'].toarray().flatten() )

    # compute the design sensitivity
    Ke = FE['Ke']
    Ue = FE['U'][FE['edofMat']].repeat(FE['n_edof'],axis=2).transpose((1,2,0))
    Ue_T = Ue.transpose((1,0,2))

    Dc_Dpenalized_elem_dens = np.sum( np.sum( 
        - Ue_T * Ke * Ue , 
        0 ) , 0 )   # Eq. (24)

    Dc_Ddv = Dc_Dpenalized_elem_dens @ OPT['Dpenalized_elem_dens_Ddv'] # Eq. (25)
    grad_c = Dc_Ddv.T
    # save these values in the OPT structure
    OPT['compliance'] = c
    OPT['grad_compliance'] = grad_c

    return c , grad_c


def compute_volume_fraction(FE,OPT,GEOM):
    #
    # This function computes the volume fraction and its sensitivities
    # based on the last geometry projection
    #
    # global FE, OPT

    # compute the volume fraction
    v_e = FE['elem_vol'] # element
    V = np.sum(v_e) # full volume
    v = np.dot( v_e , OPT['elem_dens'] ) # projected volume
    volfrac =  v/V # Eq. (16)

    # compute the design sensitivity
    Dvolfrac_Ddv = (v_e @ OPT['Delem_dens_Ddv'] )/V   # Eq. (31)
    grad_volfrac = Dvolfrac_Ddv
    
    # output
    OPT['volume_fraction'] = volfrac
    OPT['grad_volume_fraction'] = grad_volfrac
    return volfrac , grad_volfrac

def process(FE,OPT,GEOM):
    print('RUN ABAQUS')
    command = "%s job=%s user=%s interactive" % (GEOM['abaqusVersion'], GEOM['inputFile']
        ,GEOM['userFile'])
    startTime = time.perf_counter()
    subprocess.call(command, shell=True)
    #
    while(os.path.exists(GEOM['inputFile']+'.lck')):
        time.sleep(1)
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

def export_image(FE,OPT,GEOM):
    abaqus_script = 'ZScriptToReadResults_2.py'
    print('EXPORTING IMAGE')  
    command = "%s cae noGUI=%s" % (GEOM['abaqusVersion'], abaqus_script)
    os.system(command)
    print('EXPORTING IMAGE')          

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
    n_dv = 1#OPT['n_dv']
    # compute and assign the value of compliance and compliance sensitivities (analytical method)
    theta_0 = OPT['cost'][1,n_dv]
    grad_theta_0 = np.zeros((n_dv, 1))
    grad_theta_0 = OPT['cost'][1,0:n_dv]
    # This computes the volume fraction and its sensitivities
    # base on the last geometry projection
    volume_dens = OPT['cost'][1,(n_dv+500)]
    volume_sum = OPT['cost'][1,(n_dv+500+1)]
    volume_frac = volume_dens/volume_sum
    grad_volfrac = np.zeros((n_dv, 1))
    grad_volfrac = OPT['cost'][1,500:500+n_dv]/volume_sum
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
    export_image(FE,OPT,GEOM)
    export_output_raw(FE,OPT,GEOM)
    convert_output_raw2csv()
    OPT['cost'] = assign_value()
    evaluate_relevant_functions(FE,OPT,GEOM)

def nonlcon(FE, OPT, GEOM,dv):
    
    OPT['dv_old'] = OPT['dv'].copy()
    OPT['dv'] = dv.copy()
    
    # Update or perform the analysis
    if ( OPT['dv'] != OPT['dv_old'] ).any():
        update_geom_from_dv(FE, OPT, GEOM)
        point_matrix = GEOM['current_design']['point_matrix']
        bar_matrix = GEOM['current_design']['bar_matrix']
        if (FE['dim'] == 2):
            with open('ZDATA.txt', 'w') as f:
                f.write('TEST INPUT DATA OF BARS\n')
                f.write('{:5d}{:5d}{:5d}\n'.format(GEOM['current_design']['point_matrix'].shape[0],
                                              GEOM['current_design']['bar_matrix'].shape[0], FE['dim']))
                for i in range(GEOM['n_point']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}\n'.format(point_matrix[i][0],
                                                                   point_matrix[i][1], point_matrix[i][2]))
                for i in range(GEOM['n_bar']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(bar_matrix[i][0],
                                                                                     bar_matrix[i][1], bar_matrix[i][2],
                                                                                     bar_matrix[i][3], bar_matrix[i][4]))
                f.write('{:10.3f}\n'.format(GEOM['min_bar_radius']))
                f.write('{:10.3f}\n'.format(GEOM['max_bar_radius']))
                f.write('{}\n'.format(OPT['parameters']['smooth_max_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['smooth_max_param']))
                f.write('{}\n'.format(OPT['parameters']['penalization_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['penalization_param']))
        else:
            with open('ZDATA.txt', 'w') as f:
                f.write('TEST INPUT DATA OF BARS\n')
                f.write('{:5d}{:5d}{:5d}\n'.format(GEOM['initial_design']['point_matrix'].shape[0],
                    GEOM['initial_design']['bar_matrix'].shape[0], FE['dim']))
                for i in range(GEOM['n_point']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(point_matrix[i][0],
                        point_matrix[i][1], point_matrix[i][2], point_matrix[i][3]))
                for i in range(GEOM['n_bar']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(bar_matrix[i][0],
                        bar_matrix[i][1], bar_matrix[i][2],bar_matrix[i][3],bar_matrix[i][4]))
                f.write('{:10.3f}\n'.format(GEOM['min_bar_radius']))
                f.write('{:10.3f}\n'.format(GEOM['max_bar_radius']))
                f.write('{}\n'.format(OPT['parameters']['smooth_max_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['smooth_max_param']))
                f.write('{}\n'.format(OPT['parameters']['penalization_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['penalization_param']))
        performance_analysis(FE, OPT, GEOM)

    n_con   = OPT['functions']['n_func']-1 # number of constraints
    g       = np.zeros((n_con,1))

    for i in range(0,n_con):
        g[i] = OPT['functions']['f'][i+1]['value']
        g = g - OPT['functions']['constraint_limit']

    return g.flatten()


def nonlcongrad(FE,OPT,GEOM,dv):

    n_con   = OPT['functions']['n_func']-1 # number of constraints
    gradg = np.zeros((OPT['n_dv'],n_con))

    for i in range(0,n_con):
        gradg[:,i] = OPT['functions']['f'][i+1]['grad']
    
    return gradg.flatten()            

def obj(FE,OPT,GEOM,dv):

    OPT['dv_old'] = OPT['dv'].copy()  # save the previous design
    OPT['dv'] = dv.copy()  # update the design
    # If different, update or perform the analysis
    if (OPT['dv'] != OPT['dv_old']).any():
        update_geom_from_dv(FE, OPT, GEOM)
        #print(GEOM['current_design']['point_matrix'])
        #print(GEOM['current_design']['bar_matrix'])
        point_matrix = GEOM['current_design']['point_matrix']
        bar_matrix = GEOM['current_design']['bar_matrix']
        if (FE['dim'] == 2):
            with open('ZDATA.txt', 'w') as f:
                f.write('TEST INPUT DATA OF BARS\n')
                f.write('{:5d}{:5d}{:5d}\n'.format(GEOM['current_design']['point_matrix'].shape[0],
                                                   GEOM['current_design']['bar_matrix'].shape[0], FE['dim']))
                for i in range(GEOM['n_point']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}\n'.format(point_matrix[i][0],
                                                                   point_matrix[i][1], point_matrix[i][2]))
                for i in range(GEOM['n_bar']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(bar_matrix[i][0],
                                                                                     bar_matrix[i][1], bar_matrix[i][2],
                                                                                     bar_matrix[i][3],
                                                                                     bar_matrix[i][4]))
                f.write('{:10.3f}\n'.format(GEOM['min_bar_radius']))
                f.write('{:10.3f}\n'.format(GEOM['max_bar_radius']))
                f.write('{}\n'.format(OPT['parameters']['smooth_max_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['smooth_max_param']))
                f.write('{}\n'.format(OPT['parameters']['penalization_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['penalization_param']))
        else:
            with open('ZDATA.txt', 'w') as f:
                f.write('TEST INPUT DATA OF BARS\n')
                f.write('{:5d}{:5d}{:5d}\n'.format(GEOM['initial_design']['point_matrix'].shape[0],
                                                   GEOM['initial_design']['bar_matrix'].shape[0], FE['dim']))
                for i in range(GEOM['n_point']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(point_matrix[i][0],
                                                                            point_matrix[i][1], point_matrix[i][2],
                                                                            point_matrix[i][3]))
                for i in range(GEOM['n_bar']):
                    f.write('{:20.15f}{:20.15f}{:20.15f}{:20.15f}{:20.15f}\n'.format(bar_matrix[i][0],
                                                                                     bar_matrix[i][1], bar_matrix[i][2],
                                                                                     bar_matrix[i][3],
                                                                                     bar_matrix[i][4]))
                f.write('{:10.3f}\n'.format(GEOM['min_bar_radius']))
                f.write('{:10.3f}\n'.format(GEOM['max_bar_radius']))
                f.write('{}\n'.format(OPT['parameters']['smooth_max_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['smooth_max_param']))
                f.write('{}\n'.format(OPT['parameters']['penalization_scheme']))
                f.write('{:10d}\n'.format(OPT['parameters']['penalization_param']))
        performance_analysis(FE, OPT, GEOM)
    f = OPT['compliance'].flatten().copy()
    g = OPT['grad_compliance'].copy()
    return f, g

#process(FE,OPT,GEOM)
#export_output_raw(FE, OPT, GEOM)
#convert_output_raw2csv()
#OPT['cost'] = assign_value()
#evaluate_relevant_functions(FE, OPT, GEOM)


