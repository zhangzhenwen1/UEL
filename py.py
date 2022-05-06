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

def ReadODB(filepath):
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
    # obtain the RF of reference points
    frame=odb.steps['Step-1'].frames[-1]
    REF=[]
    REF.append( myAssembly.nodeSets['REF-X'])
    REF.append(myAssembly.nodeSets['REF-Y'])
    REF.append(myAssembly.nodeSets['REF-Z'])
    REF.append(myAssembly.nodeSets['REF-XY'])
    REF.append(myAssembly.nodeSets['REF-XZ'])
    REF.append(myAssembly.nodeSets['REF-YZ'])
    field = frame.fieldOutputs['RF']
    data=[]
    for r in REF:
        region=field.getSubset( region=r )
        regionValue = region.values
        for p in regionValue:
            data.append(np.array(p.data))
    data=np.array(data)
    np.savetxt('C:/repo/UEL/DSTRESS.out', data) 

    volume=np.array(volume,dtype=np.float64)
    vol=np.copy(volume[0])
    #np.savetxt('C:/repo/UEL/DSTRESS.out', Vsum) 
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
    string='Ref-X, 1, 1, '+str(strainX)+'\n'
    string=string+'Ref-Y, 2, 2, '+str(strainY)+'\n'
    string=string+'Ref-Z, 3, 3, '+str(strainZ)+'\n'
    string=string+'Ref-XY, 1, 2, '+str(strainXY)+'\n'
    string=string+'Ref-XZ, 1, 1, '+str(strainXZ)+'\n'
    string=string+'Ref-XZ, 3, 3, '+str(strainXZ)+'\n'
    string=string+'Ref-YZ, 2, 3, '+str(strainYZ)+'\n'
    BOUNDARY.append('\n*Step, name=Step-1, nlgeom=NO\n')
    BOUNDARY.append('*Static\n')
    BOUNDARY.append('** OUTPUT REQUESTS\n')
    BOUNDARY.append('*Restart, write, frequency=1\n')
    BOUNDARY.append('*Output, field, variable=ALL\n')
    BOUNDARY.append('*Output, history, frequency=1\n')
    BOUNDARY.append('*Boundary\n')
    BOUNDARY.append(string)
    BOUNDARY.append('*End Step')

def MonitorLogFile(JobName):
    logpath='C:/repo/UEL/'+JobName+'.log'
    while 1:
        try:
            file=open(logpath,'r')
        except:
            time.sleep(1)
            continue
        try:
            lastline=file.readlines()[-1]
        except:
            time.sleep(1)
            continue
        if 'COMPLETED'in lastline:
            print (' python INFO: Abaqus/Analysis '+JobName+' COMPLETED')
            ReadODB(r'C:/repo/UEL/'+JobName+'.odb')
            break
        elif 'errors'in lastline:
            print ('!!! python ERROR: Abaqus/Analysis '+JobName+' exited with errors')
            break
        time.sleep(1)
    file.close()
 

print ('------------------------------------------------------------------------')
print ('----------------------   RVE COMPUTAION PROCESS   ----------------------')
print ('------------------------------------------------------------------------')

JobName='RVE'

BOUNDARY=[]

try:
    JSTEP = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_footer=9)
except:
    print (' !!! python ERROR: no data of DSTRAN')
    exit(1)

#print ('python INFO: JSTEP is  '), JSTEP[0]
KINC = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_header=4,skip_footer=8)
#print ('python INFO: KINC is  '), KINC
NOEL = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_header=5,skip_footer=7)
#print ('python INFO: NOEL is  '), NOEL
NPT = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt',dtype=int,skip_header=6,skip_footer=6)
#print ('python INFO: NPT is  '), NPT
DSTRAN = np.genfromtxt('C:/repo/UEL/DSTRAIN.txt', delimiter=28,dtype=np.float64,skip_header=7)
print ('python INFO: DSTRAN is  '), DSTRAN

load(DSTRAN,BOUNDARY)

JobName_init=JobName+'_'+str(NOEL)+'_'+str(NPT)
JobName_res=JobName_init+'_restart'+'_'+str(KINC)
if ( KINC == 0 ):
    filepath='C:/repo/UEL/'+JobName_init+'.inp'
    
    REPLACE('C:/repo/UEL/RVE_template.inp',filepath, ['**BOUNDARY_XYZ'], [''.join(BOUNDARY)])
    REPLACE('C:/repo/UEL/RVE_submit_job_template.bat','C:/repo/UEL/RVE_submit_job.bat', ['JobName'], [JobName_init])
    
    submit_file='C:/repo/UEL/RVE_submit_job.bat'
    try:
        os.remove('C:/repo/UEL/'+JobName_init+'.log')
    except:
        print ('LOG file does not exist')
    os.system(submit_file)
    MonitorLogFile(JobName_init)

elif ( KINC == 1 ):
    filepath='C:/repo/UEL/'+JobName_res+'.inp'
    job_submit=JobName_res+' oldjob='+JobName_init

    REPLACE('C:/repo/UEL/RVE_Restart-Template.inp',filepath, ['STEP_REPLACE','INC_REPLACE','BOUNDARY_XYZ'], [KINC*2,1,''.join(BOUNDARY)])
    REPLACE('C:/repo/UEL/RVE_submit_job_template.bat','C:/repo/UEL/RVE_submit_job.bat', ['JobName'], [job_submit])
    
    submit_file='C:/repo/UEL/RVE_submit_job.bat'
    try:
        os.remove('C:/repo/UEL/'+JobName_res+'.log')
    except:
        print ('LOG file does not exist')
    os.system(submit_file)
    MonitorLogFile(JobName_res)

else:
    JobName_init=JobName_init+'_restart'+'_'+str(KINC-1)
    filepath='C:/repo/UEL/'+JobName_res+'.inp'
    job_submit=JobName_res+' oldjob='+JobName_init

    REPLACE('C:/repo/UEL/RVE_Restart-Template.inp',filepath, ['STEP_REPLACE','INC_REPLACE','BOUNDARY_XYZ'], [KINC*2,1,''.join(BOUNDARY)])
    REPLACE('C:/repo/UEL/RVE_submit_job_template.bat','C:/repo/UEL/RVE_submit_job.bat', ['JobName'], [job_submit])
    
    submit_file='C:/repo/UEL/RVE_submit_job.bat'
    try:
        os.remove('C:/repo/UEL/'+JobName_res+'.log')
    except:
        print ('LOG file does not exist')
    os.system(submit_file)
    MonitorLogFile(JobName_res)

exit(0)