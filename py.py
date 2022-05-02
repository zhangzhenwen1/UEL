# -* - coding:UTF-8 -*-
from distutils import core
from abaqusConstants import*
from odbAccess import*
from copy import deepcopy
import numpy as np
import os
import time
import re

np.set_printoptions(precision=16)

def REPLACE(file, new_file, old_str, new_str):
    with open(file, "r") as f1, open(new_file, "w") as f2:
        i=0
        for line in f1:
            for text in old_str:
                if text in line:
                    line = line.replace(text, str(new_str[i]))
                    i=i+1
            f2.write(line)

def ReadODB(filepath,DSTRAN,BOUNDARY):
    odb = openOdb(path=filepath)
    myAssembly = odb.rootAssembly
    # Part instance determine how many instances
    # for instanceName in odb.rootAssembly.instances.keys():
    #    print instanceName
    # step1 = odb.steps.values()[0]
    # print step1.name
    #Frames  the last frame
    #print myAssembly.nodeSets.keys()
    coord=[]
    nodeLable=[]
    jacobi=[]
    volume=[]
    for frame in odb.steps['perturbation'].frames:
        if frame.frameId == 0 :
            # Reading field output data
            field=frame.fieldOutputs['EVOL'].values
            vol=[]
            for v in field:
                vol.append (v.data)
            vol=np.array(vol)
            vsum=np.sum(vol)
            volume.append(vsum)

            for p in myAssembly.nodeSets.keys():
                if re.search(r'X+\w+PAIR', p,re.I):
                    nodesets=myAssembly.nodeSets[p]
                    # Reading field output data
                    field=frame.fieldOutputs['COORD']
                    region = field.getSubset(region=nodesets)
                    regionValue=region.values
                    data_lable=[]
                    data_coord=[]
                    for q in regionValue:
                        data_lable.append (q.nodeLabel)
                        data_coord.append(q.data)
                    nodeLable.append(data_lable)
                    coord.append(data_coord)
            nodeLable=np.array(nodeLable)
            coord=np.array(coord,dtype=np.float64)                                        

        elif frame.frameId == 1 :
            # Reading field output data
            field = frame.fieldOutputs['S']
            region=field.getSubset( position=INTEGRATION_POINT )
            regionValue = region.values
            data=[]
            for p in regionValue:
                data.append(p.data)
            data=np.array(data)
            sigmaMean=np.mean(data, axis=0)
            #print sigmaMean
            jacobi.append(np.mean(data, axis=0))
    volume=np.array(volume,dtype=np.float64)
    vol=np.copy(volume[0])
    #np.savetxt('C:/repo/UEL/test.out', Vsum) 
    jacobi=np.array(jacobi)/vol
    jacobi=np.triu(jacobi.T)
    jacobi += jacobi.T - np.diag(jacobi.diagonal())
    np.savetxt('C:/repo/UEL/DDSDDE.out', jacobi)

def load(DSTRAN,BOUNDARY):
    DSTRAN=np.array(DSTRAN,dtype=np.float64)
    strainX=np.copy(DSTRAN[0])
    strainY=np.copy(DSTRAN[1])
    strainZ=np.copy(DSTRAN[2])
    strainXY=np.copy(DSTRAN[3])
    strainXZ=np.copy(DSTRAN[4])
    strainYZ=np.copy(DSTRAN[5])
    # For faces x=0 and x=a (excluding edges)
    string='nodeSet_X_minus, 1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_X_plus, 1,1,'+str(strainX)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_X_plus, 2,3'+'\n'
    BOUNDARY.append(string)
    # For faces y=0 and y=b (excluding edges)
    string='nodeSet_Y_minus,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Y_plus,1,1, '+str(strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Y_plus,2,2, '+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Y_plus,3,3'+'\n'
    BOUNDARY.append(string)
    # For faces z=0 and z=c (excluding edges
    string='nodeSet_Z_minus,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Z_plus,1,1, '+str(strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Z_plus,2,2, '+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_Z_plus,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    # For edges parallel to the x-axis (excluding vertices)
    string='nodeSet_parall_x_y0_z0,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y0_z1,1,1,'+str(strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y0_z1,2,2,'+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y0_z1,3,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z0,1,1,'+str(strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z0,2,2,'+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z0,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z1,1,1,'+str(strainXY+strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z1,2,2,'+str(strainY+strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_x_y1_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    # For edges parallel to the y-axis (excluding vertices)
    string='nodeSet_parall_y_x0_z0,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x1_z0,1,1,'+str(strainX)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x1_z0,2,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x0_z1,1,1,'+str(strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x0_z1,2,2,'+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x0_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x1_z1,1,1,'+str(strainX+strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x1_z1,2,2,'+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_y_x1_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    # For edges parallel to the z-axis (excluding vertices)
    string='nodeSet_parall_z_x0_y0,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x1_y0,1,1,'+str(strainX)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x1_y0,2,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x0_y1,1,1,'+str(strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x0_y1,2,2,'+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x0_y1,3,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x1_y1,1,1,'+str(strainX+strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x1_y1,2,2,'+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_parall_z_x1_y1,3,3'+'\n'
    BOUNDARY.append(string)
    # For the vertices
    string='nodeSet_vertice_x0_y0_z0,1,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y0_z0,1,1,'+str(strainX)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y0_z0,2,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z0,1,1,'+str(strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z0,2,2,'+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z0,3,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z0,1,1,'+str(strainX+strainXY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z0,2,2,'+str(strainY)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z0,3,3'+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y0_z1,1,1,'+str(strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y0_z1,2,2,'+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y0_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y0_z1,1,1,'+str(strainX+strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y0_z1,2,2,'+str(strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y0_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z1,1,1,'+str(strainXY+strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z1,2,2,'+str(strainY+strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x0_y1_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z1,1,1,'+str(strainX+strainXY+strainXZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z1,2,2,'+str(strainY+strainYZ)+'\n'
    BOUNDARY.append(string)
    string='nodeSet_vertice_x1_y1_z1,3,3,'+str(strainZ)+'\n'
    BOUNDARY.append(string)
 

print ('------------------------------------------------------------------------')
print ('----------------------   RVE COMPUTAION PROCESS   ----------------------')
print ('------------------------------------------------------------------------')
JSTEP = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_footer=10)
print ('JSTEP is  '), JSTEP
KINC = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_header=4,skip_footer=6)
print ('KINC is  '), KINC 
DSTRAN = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt', delimiter=28,dtype=np.float64,skip_header=5)
print ('DSTRAN is  '), DSTRAN

BOUNDARY=[]

filepath='C:/repo/UEL/RVE-PBC.log'
while 1:
    time.sleep(1)
    file=open(filepath,'r')
    lastline=file.readlines()[-1]
    if 'COMPLETED'in lastline:
        print ('Abaqus/Analysis COMPLETED')
        ReadODB(r'C:/repo/UEL/RVE-PBC-perturbation-pair.odb',DSTRAN,BOUNDARY)
        break
    elif 'errors'in lastline:
        print ('Abaqus/Analysis exited with errors')
        break
    file.close()
load(DSTRAN,BOUNDARY)
REPLACE('C:/repo/UEL/RVE_Restart-Template.inp','C:/repo/UEL/RVE_Restart.inp', ['STEP_REPLACE','INC_REPLACE','BOUNDARY_XYZ'], [JSTEP,KINC,''.join(BOUNDARY)])
job_submit='C:/repo/UEL/RVE-restart.bat'
os.system(job_submit)

print ("Hello, world!")