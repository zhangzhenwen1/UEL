# -* - coding:UTF-8 -*-
from abaqusConstants import*
from odbAccess import*
from copy import deepcopy
import numpy as np
import os
import time

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

def ReadODB(filepath):
    odb = openOdb(path=filepath)
    myAssembly = odb.rootAssembly
    # Part instance determine how many instances
    # for instanceName in odb.rootAssembly.instances.keys():
    #    print instanceName
    # step1 = odb.steps.values()[0]
    # print step1.name
    #Frames  the last frame
    jacobi=[]
    volume=[]
    for frame in odb.steps['perturbation'].frames:
        if frame.frameId == 0 :
            field=frame.fieldOutputs['EVOL'].values
            vol=[]
            for v in field:
                vol.append (v.data)
            vol=np.array(vol)
            vsum=np.sum(vol)
            volume.append(vsum)
        elif frame.frameId == 1 :
            # Reading field output data
            #print (frame.description), (':  ')
            stressField = frame.fieldOutputs['S']
            stress=stressField.getSubset( position=INTEGRATION_POINT )
            stressValues = stress.values
            sigmaRVE=[]
            for v in stressValues:
                sigmaRVE.append(v.data)
            sigmaRVE=np.array(sigmaRVE)
            sigmaMean=np.mean(sigmaRVE, axis=0)
            #print sigmaMean
            jacobi.append(sigmaMean)
    volume=np.array(volume,dtype=np.float64)
    vol=np.copy(volume[0])
    #np.savetxt('C:/repo/UEL/test.out', Vsum) 
    jacobi=np.array(jacobi)/vol
    jacobi=np.triu(jacobi.T)
    jacobi += jacobi.T - np.diag(jacobi.diagonal())
    np.savetxt('C:/repo/UEL/DDSDDE.out', jacobi)

print ('------------------------------------------------------------------------')
print ('----------------------   RVE COMPUTAION PROCESS   ----------------------')
print ('------------------------------------------------------------------------')
JSTEP = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_footer=10)
print ('JSTEP is  '), JSTEP
KINC = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_header=4,skip_footer=6)
print ('KINC is  '), KINC 
DSTRAN = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt', delimiter=28,dtype=np.float64,skip_header=5)
print ('DSTRAN is  '), DSTRAN

REPLACE('C:/repo/UEL/RVE_Restart-Template.inp','C:/repo/UEL/RVE_Restart.inp', ['STEP_REPLACE','INC_REPLACE'], [JSTEP,KINC])

job_submit='C:/repo/UEL/RVE-restart.bat'
os.system(job_submit)

filepath='C:/repo/UEL/RVE-PBC.log'
while 1:
    time.sleep(1)
    file=open(filepath,'r')
    lastline=file.readlines()[-1]
    if 'COMPLETED'in lastline:
        print ('Abaqus/Analysis COMPLETED')
        ReadODB(r'C:/repo/UEL/RVE-PBC-perturbation.odb')
        break
    elif 'errors'in lastline:
        print ('Abaqus/Analysis exited with errors')
        break
    file.close()




print ("Hello, world!")