from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
import os
import json
import sys


current_path = os.getcwd()
#
working_directory = os.getcwd()
os.chdir(working_directory)

session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
#=========================== Read Input Parameter ===========================#
execfile(os.path.join(current_path,'Input_parameter.py'))

#========================== CREATE RESULTS FOLDER =============================#
img_folder = 'Exported_Images'
if os.path.isdir(os.path.join(working_directory ,img_folder)):
    results_path = os.path.join(working_directory ,img_folder)
else:
    #       create pah
    results_path = os.path.join(working_directory ,img_folder)
    os.makedirs(results_path)    


#================================== Results ==================================#

o1 = session.openOdb(name= outputFile+'.odb')
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
odb = session.odbs[outputFile+'.odb']
#==================================== VIEW ===================================#
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
session.graphicsOptions.setValues(backgroundColor='#FFFFFF', 
    backgroundBottomColor='#FFFFFF')
session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
session.viewports['Viewport: 1'].odbDisplay.commonOptions.setValues(
    visibleEdges=NONE, deformationScaling=UNIFORM, uniformScaleFactor=0.0001)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(triad=OFF, 
    legend=OFF, title=OFF, state=OFF, annotations=OFF, compass=OFF)
session.viewports['Viewport: 1'].viewportAnnotationOptions.setValues(
    legendFont='-*-verdana-medium-r-normal-*-*-180-*-*-p-*-*-*') 
#---------------------------Results of elements ------------------------------#
#-------------------------- Save von Mises stresses --------------------------#
ID_simulation = '%03d' % iteration
S_List=[variable]
for S_Name in S_List:
    aux = COMPONENT
    if S_Name == 'Mises':
        aux = INVARIANT
        #session.viewports['Viewport: 1'].setValues(displayedObject=o1)
        #session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
        #    CONTOURS_ON_DEF, ))
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel='S', outputPosition=ELEMENT_CENTROID, refinement=(aux,
                                                                             S_Name), )
    #   session.viewports['Viewport: 1'].view.setValues(session.views['Front'])
    else:
        #session.viewports['Viewport: 1'].setValues(displayedObject=o1)
        session.viewports['Viewport: 1'].odbDisplay.setPrimaryVariable(
            variableLabel=S_Name, outputPosition=ELEMENT_CENTROID, )

    session.printOptions.setValues(vpDecorations=OFF, reduceColors=False)                                       # Print von Mises stress image
    session.printToFile(fileName=img_folder+'/'+S_Name+'_'+ID_simulation, format=PNG, canvasObjects=(
        session.viewports['Viewport: 1'], ))

    