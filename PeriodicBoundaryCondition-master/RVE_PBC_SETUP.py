import abaqus
import abaqusConstants
import customKernel
import customKernelSerialize
from math import sqrt
import numpy as np

fields=(('model name','RVE-PBC'),('instance name','PART-1-1'),('x min:','0.0'),('x max:','1.0'),('y min:','0.0'),('y max:','1.0'),('z min:','0.0'),('z max:','1.0'))
modelName,insName,xmin,xmax,ymin,ymax,zmin,zmax=getInputs(fields=fields,label="Input the parameters")
modelName=str(modelName);insName=str(insName);xmin=float(xmin);xmax=float(xmax);ymin=float(ymin);ymax=float(ymax);zmin=float(zmin);zmax=float(zmax)

m=mdb.models[modelName]
r=m.rootAssembly
tor=0.001

refpoints=r.referencePoints
if len(refpoints)<6:
    r.ReferencePoint(point=(0.5,0.5,1.0))
    r.ReferencePoint(point=(0.5,1.0,0.5))
    r.ReferencePoint(point=(1.0,0.5,0.5))
    r.ReferencePoint(point=(0.6,0.6,1.0))
    r.ReferencePoint(point=(0.6,1.0,0.6))
    r.ReferencePoint(point=(1.0,0.6,0.6))
for i in range(len(refpoints)):
    r.Set(name='Ref-'+str(i+1),referencePoints=(refpoints.values()[i],))

node=r.instances[insName].nodes
nodeSet_X_minus=[]
nodeSet_X_plus=[]
nodeSet_Y_minus=[]
nodeSet_Y_plus=[]
nodeSet_Z_minus=[]
nodeSet_Z_plus=[]   
for i in range(len(node)):
    x=node[i].coordinates[0]
    y=node[i].coordinates[1]
    z=node[i].coordinates[2]
    if (abs(x-xmin)<tor):
        nodeSet_X_minus.append(i+1)
    if (abs(x-xmax)<tor):
        nodeSet_X_plus.append(i+1)
    if (abs(y-ymin)<tor):
        nodeSet_Y_minus.append(i+1)
    if (abs(y-ymax)<tor):
        nodeSet_Y_plus.append(i+1) 
    if (abs(z-zmin)<tor):
        nodeSet_Z_minus.append(i+1)
    if (abs(z-zmax)<tor):
        nodeSet_Z_plus.append(i+1)
nodeSet_X_minus=set(nodeSet_X_minus)
nodeSet_X_plus=set(nodeSet_X_plus)
nodeSet_Y_minus=set(nodeSet_Y_minus)
nodeSet_Y_plus=set(nodeSet_Y_plus)
nodeSet_Z_minus=set(nodeSet_Z_minus)
nodeSet_Z_plus=set(nodeSet_Z_plus)

nodeSet_vertice_x0_y0_z0= nodeSet_X_minus & nodeSet_Y_minus & nodeSet_Z_minus
nodeSet_vertice_x1_y0_z0= nodeSet_X_plus & nodeSet_Y_minus & nodeSet_Z_minus
nodeSet_vertice_x0_y1_z0= nodeSet_X_minus & nodeSet_Y_plus & nodeSet_Z_minus
nodeSet_vertice_x0_y0_z1= nodeSet_X_minus & nodeSet_Y_minus & nodeSet_Z_plus
nodeSet_vertice_x1_y1_z0= nodeSet_X_plus & nodeSet_Y_plus & nodeSet_Z_minus
nodeSet_vertice_x0_y1_z1= nodeSet_X_minus & nodeSet_Y_plus & nodeSet_Z_plus
nodeSet_vertice_x1_y0_z1= nodeSet_X_plus & nodeSet_Y_minus & nodeSet_Z_plus
nodeSet_vertice_x1_y1_z1= nodeSet_X_plus & nodeSet_Y_plus & nodeSet_Z_plus

nodeSet_parall_x_y0_z0= nodeSet_Y_minus & nodeSet_Z_minus
nodeSet_parall_x_y0_z1= nodeSet_Y_minus & nodeSet_Z_plus
nodeSet_parall_x_y1_z0= nodeSet_Y_plus & nodeSet_Z_minus
nodeSet_parall_x_y1_z1= nodeSet_Y_plus & nodeSet_Z_plus

nodeSet_parall_y_x0_z0= nodeSet_X_minus & nodeSet_Z_minus
nodeSet_parall_y_x0_z1= nodeSet_X_minus & nodeSet_Z_plus
nodeSet_parall_y_x1_z0= nodeSet_X_plus & nodeSet_Z_minus
nodeSet_parall_y_x1_z1= nodeSet_X_plus & nodeSet_Z_plus

nodeSet_parall_z_x0_y0= nodeSet_X_minus & nodeSet_Y_minus
nodeSet_parall_z_x0_y1= nodeSet_X_minus & nodeSet_Y_plus
nodeSet_parall_z_x1_y0= nodeSet_X_plus & nodeSet_Y_minus
nodeSet_parall_z_x1_y1= nodeSet_X_plus & nodeSet_Y_plus

nodeSet_X_minus=nodeSet_X_minus-nodeSet_parall_y_x0_z0-nodeSet_parall_y_x0_z1-nodeSet_parall_z_x0_y0-nodeSet_parall_z_x0_y1
nodeSet_X_plus =nodeSet_X_plus-nodeSet_parall_y_x1_z0-nodeSet_parall_y_x1_z1-nodeSet_parall_z_x1_y0-nodeSet_parall_z_x1_y1
nodeSet_Y_minus=nodeSet_Y_minus-nodeSet_parall_x_y0_z0-nodeSet_parall_x_y0_z1-nodeSet_parall_z_x0_y0-nodeSet_parall_z_x1_y0
nodeSet_Y_plus =nodeSet_Y_plus-nodeSet_parall_x_y1_z0-nodeSet_parall_x_y1_z1-nodeSet_parall_z_x0_y1-nodeSet_parall_z_x1_y1
nodeSet_Z_minus=nodeSet_Z_minus-nodeSet_parall_x_y0_z0-nodeSet_parall_x_y1_z0-nodeSet_parall_y_x0_z0-nodeSet_parall_y_x1_z0
nodeSet_Z_plus =nodeSet_Z_plus-nodeSet_parall_x_y0_z1-nodeSet_parall_x_y1_z1-nodeSet_parall_y_x0_z1-nodeSet_parall_y_x1_z1

nodeSet_parall_x_y0_z0=nodeSet_parall_x_y0_z0 - nodeSet_vertice_x0_y0_z0 - nodeSet_vertice_x1_y0_z0
nodeSet_parall_x_y0_z1=nodeSet_parall_x_y0_z1 - nodeSet_vertice_x0_y0_z1 - nodeSet_vertice_x1_y0_z1
nodeSet_parall_x_y1_z0=nodeSet_parall_x_y1_z0 - nodeSet_vertice_x0_y1_z0 - nodeSet_vertice_x1_y1_z0
nodeSet_parall_x_y1_z1=nodeSet_parall_x_y1_z1 - nodeSet_vertice_x0_y1_z1 - nodeSet_vertice_x1_y1_z1
nodeSet_parall_y_x0_z0=nodeSet_parall_y_x0_z0 - nodeSet_vertice_x0_y0_z0 - nodeSet_vertice_x0_y1_z0
nodeSet_parall_y_x0_z1=nodeSet_parall_y_x0_z1 - nodeSet_vertice_x0_y0_z1 - nodeSet_vertice_x0_y1_z1
nodeSet_parall_y_x1_z0=nodeSet_parall_y_x1_z0 - nodeSet_vertice_x1_y0_z0 - nodeSet_vertice_x1_y1_z0
nodeSet_parall_y_x1_z1=nodeSet_parall_y_x1_z1 - nodeSet_vertice_x1_y0_z1 - nodeSet_vertice_x1_y1_z1
nodeSet_parall_z_x0_y0=nodeSet_parall_z_x0_y0 - nodeSet_vertice_x0_y0_z0 - nodeSet_vertice_x0_y0_z1
nodeSet_parall_z_x0_y1=nodeSet_parall_z_x0_y1 - nodeSet_vertice_x0_y1_z0 - nodeSet_vertice_x0_y1_z1
nodeSet_parall_z_x1_y0=nodeSet_parall_z_x1_y0 - nodeSet_vertice_x1_y0_z0 - nodeSet_vertice_x1_y0_z1
nodeSet_parall_z_x1_y1=nodeSet_parall_z_x1_y1 - nodeSet_vertice_x1_y1_z0 - nodeSet_vertice_x1_y1_z1

nodeSet_X_minus =list(nodeSet_X_minus)
nodeSet_X_plus  =list(nodeSet_X_plus)
nodeSet_Y_minus =list(nodeSet_Y_minus)
nodeSet_Y_plus  =list(nodeSet_Y_plus)
nodeSet_Z_minus =list(nodeSet_Z_minus)
nodeSet_Z_plus  =list(nodeSet_Z_plus)

nodeSet_vertice_x0_y0_z0=list(nodeSet_vertice_x0_y0_z0)
nodeSet_vertice_x1_y0_z0=list(nodeSet_vertice_x1_y0_z0)
nodeSet_vertice_x0_y1_z0=list(nodeSet_vertice_x0_y1_z0)
nodeSet_vertice_x0_y0_z1=list(nodeSet_vertice_x0_y0_z1)
nodeSet_vertice_x1_y1_z0=list(nodeSet_vertice_x1_y1_z0)
nodeSet_vertice_x0_y1_z1=list(nodeSet_vertice_x0_y1_z1)
nodeSet_vertice_x1_y0_z1=list(nodeSet_vertice_x1_y0_z1)
nodeSet_vertice_x1_y1_z1=list(nodeSet_vertice_x1_y1_z1)

nodeSet_parall_x_y0_z0=list(nodeSet_parall_x_y0_z0)
nodeSet_parall_x_y0_z1=list(nodeSet_parall_x_y0_z1)
nodeSet_parall_x_y1_z0=list(nodeSet_parall_x_y1_z0)
nodeSet_parall_x_y1_z1=list(nodeSet_parall_x_y1_z1)

nodeSet_parall_y_x0_z0=list(nodeSet_parall_y_x0_z0)
nodeSet_parall_y_x0_z1=list(nodeSet_parall_y_x0_z1)
nodeSet_parall_y_x1_z0=list(nodeSet_parall_y_x1_z0)
nodeSet_parall_y_x1_z1=list(nodeSet_parall_y_x1_z1)

nodeSet_parall_z_x0_y0=list(nodeSet_parall_z_x0_y0)
nodeSet_parall_z_x0_y1=list(nodeSet_parall_z_x0_y1)
nodeSet_parall_z_x1_y0=list(nodeSet_parall_z_x1_y0)
nodeSet_parall_z_x1_y1=list(nodeSet_parall_z_x1_y1)

for n in range(len(nodeSet_X_plus)):
    x0=node[nodeSet_X_minus[n]-1].coordinates[0]
    y0=node[nodeSet_X_minus[n]-1].coordinates[1]
    z0=node[nodeSet_X_minus[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_X_minus[n]-1:nodeSet_X_minus[n]],name='Set-X-MINUS-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(nodeSet_X_plus)):
        x1=node[nodeSet_X_plus[j]-1].coordinates[0]
        y1=node[nodeSet_X_plus[j]-1].coordinates[1]
        z1=node[nodeSet_X_plus[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j
    r.Set(nodes=node[nodeSet_X_plus[index]-1:nodeSet_X_plus[index]],name='Set-X-PLUS-'+str(n+1))
    m.Equation(name='Eq-X-U-'+str(n+1),terms=((1,'Set-X-MINUS-'+str(n+1),1),(1,'Ref-X',1),(-1,'Set-X-PLUS-'+str(n+1),1)))
    m.Equation(name='Eq-X-V-'+str(n+1),terms=((1,'Set-X-MINUS-'+str(n+1),2),(-1,'Set-X-PLUS-'+str(n+1),2)))
    m.Equation(name='Eq-X-W-'+str(n+1),terms=((1,'Set-X-MINUS-'+str(n+1),3),(-1,'Set-X-PLUS-'+str(n+1),3)))
for n in range(len(nodeSet_Y_plus)):
    x0=node[nodeSet_Y_minus[n]-1].coordinates[0]
    y0=node[nodeSet_Y_minus[n]-1].coordinates[1]
    z0=node[nodeSet_Y_minus[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_Y_minus[n]-1:nodeSet_Y_minus[n]],name='Set-Y-MINUS-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(nodeSet_Y_plus)):
        x1=node[nodeSet_Y_plus[j]-1].coordinates[0]
        y1=node[nodeSet_Y_plus[j]-1].coordinates[1]
        z1=node[nodeSet_Y_plus[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j
    r.Set(nodes=node[nodeSet_Y_plus[index]-1:nodeSet_Y_plus[index]],name='Set-Y-PLUS-'+str(n+1))
    m.Equation(name='Eq-Y-U-'+str(n+1),terms=((1,'Set-Y-MINUS-'+str(n+1),1),(1,'Ref-XY',1),(-1,'Set-Y-PLUS-'+str(n+1),1)))
    m.Equation(name='Eq-Y-V-'+str(n+1),terms=((1,'Set-Y-MINUS-'+str(n+1),2),(1,'Ref-Y',2),(-1,'Set-Y-PLUS-'+str(n+1),2)))
    m.Equation(name='Eq-Y-W-'+str(n+1),terms=((1,'Set-Y-MINUS-'+str(n+1),3),(-1,'Set-Y-PLUS-'+str(n+1),3)))
for n in range(len(nodeSet_Z_plus)):
    x0=node[nodeSet_Z_minus[n]-1].coordinates[0]
    y0=node[nodeSet_Z_minus[n]-1].coordinates[1]
    z0=node[nodeSet_Z_minus[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_Z_minus[n]-1:nodeSet_Z_minus[n]],name='Set-Z-MINUS-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(nodeSet_Z_plus)):
        x1=node[nodeSet_Z_plus[j]-1].coordinates[0]
        y1=node[nodeSet_Z_plus[j]-1].coordinates[1]
        z1=node[nodeSet_Z_plus[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j
    r.Set(nodes=node[nodeSet_Z_plus[index]-1:nodeSet_Z_plus[index]],name='Set-Z-PLUS-'+str(n+1))
    m.Equation(name='Eq-Z-U-'+str(n+1),terms=((1,'Set-Z-MINUS-'+str(n+1),1),(1,'Ref-XZ',1),(-1,'Set-Z-PLUS-'+str(n+1),1)))
    m.Equation(name='Eq-Z-V-'+str(n+1),terms=((1,'Set-Z-MINUS-'+str(n+1),2),(1,'Ref-YZ',2),(-1,'Set-Z-PLUS-'+str(n+1),2)))
    m.Equation(name='Eq-Z-W-'+str(n+1),terms=((1,'Set-Z-MINUS-'+str(n+1),3),(1,'Ref-Z',3),(-1,'Set-Z-PLUS-'+str(n+1),3)))

# For edges parallel to the x-axis
for n in range(len(nodeSet_parall_x_y1_z0)):
    x0=node[nodeSet_parall_x_y1_z0[n]-1].coordinates[0]
    y0=node[nodeSet_parall_x_y1_z0[n]-1].coordinates[1]
    z0=node[nodeSet_parall_x_y1_z0[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_parall_x_y1_z0[n]-1:nodeSet_parall_x_y1_z0[n]],name='Set-parall_x-10-'+str(n+1))
    index=0
    index2=0
    for j in range(len(nodeSet_parall_x_y0_z0)):
        x1=node[nodeSet_parall_x_y0_z0[j]-1].coordinates[0]
        y1=node[nodeSet_parall_x_y0_z0[j]-1].coordinates[1]
        z1=node[nodeSet_parall_x_y0_z0[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
            r.Set(nodes=node[nodeSet_parall_x_y0_z0[j]-1:nodeSet_parall_x_y0_z0[j]],name='Set-parall_x-0-'+str(n+1))
    m.Equation(name='Eq-P_X_10-U-'+str(n+1),terms=((1,'Set-parall_x-0-'+str(n+1),1),(1,'Ref-XY',1),(-1,'Set-parall_x-10-'+str(n+1),1)))
    m.Equation(name='Eq-P_X_10-V-'+str(n+1),terms=((1,'Set-parall_x-0-'+str(n+1),2),(1,'Ref-Y',2),(-1,'Set-parall_x-10-'+str(n+1),2)))
    m.Equation(name='Eq-P_X_10-W-'+str(n+1),terms=((1,'Set-parall_x-0-'+str(n+1),3),(-1,'Set-parall_x-10-'+str(n+1),3)))
    for j in range(len(nodeSet_parall_x_y1_z1)):
        x1=node[nodeSet_parall_x_y1_z1[j]-1].coordinates[0]
        y1=node[nodeSet_parall_x_y1_z1[j]-1].coordinates[1]
        z1=node[nodeSet_parall_x_y1_z1[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
            for k in range(len(nodeSet_parall_x_y0_z1)):
                x2=node[nodeSet_parall_x_y0_z1[k]-1].coordinates[0]
                y2=node[nodeSet_parall_x_y0_z1[k]-1].coordinates[1]
                z2=node[nodeSet_parall_x_y0_z1[k]-1].coordinates[2]
                distance=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                if (distance < 1.001):
                    index2=k
    r.Set(nodes=node[nodeSet_parall_x_y1_z1[index]-1:nodeSet_parall_x_y1_z1[index]],name='Set-parall_x-11-'+str(n+1))
    m.Equation(name='Eq-P_X_11-U-'+str(n+1),terms=((1,'Set-parall_x-10-'+str(n+1),1),(1,'Ref-XZ',1),(-1,'Set-parall_x-11-'+str(n+1),1)))
    m.Equation(name='Eq-P_X_11-V-'+str(n+1),terms=((1,'Set-parall_x-10-'+str(n+1),2),(1,'Ref-YZ',2),(-1,'Set-parall_x-11-'+str(n+1),2)))
    m.Equation(name='Eq-P_X_11-W-'+str(n+1),terms=((1,'Set-parall_x-10-'+str(n+1),3),(1,'Ref-Z',3),(-1,'Set-parall_x-11-'+str(n+1),3)))
    
    r.Set(nodes=node[nodeSet_parall_x_y0_z1[index2]-1:nodeSet_parall_x_y0_z1[index2]],name='Set-parall_x-01-'+str(n+1))
    m.Equation(name='Eq-P_X_01-U-'+str(n+1),terms=((1,'Set-parall_x-01-'+str(n+1),1),(1,'Ref-XY',1),(-1,'Set-parall_x-11-'+str(n+1),1)))
    m.Equation(name='Eq-P_X_01-V-'+str(n+1),terms=((1,'Set-parall_x-01-'+str(n+1),2),(1,'Ref-Y',2),(-1,'Set-parall_x-11-'+str(n+1),2)))
    m.Equation(name='Eq-P_X_01-W-'+str(n+1),terms=((1,'Set-parall_x-01-'+str(n+1),3),(-1,'Set-parall_x-11-'+str(n+1),3)))

# For edges parallel to the y-axis
for n in range(len(nodeSet_parall_y_x1_z0)):
    x0=node[nodeSet_parall_y_x1_z0[n]-1].coordinates[0]
    y0=node[nodeSet_parall_y_x1_z0[n]-1].coordinates[1]
    z0=node[nodeSet_parall_y_x1_z0[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_parall_y_x1_z0[n]-1:nodeSet_parall_y_x1_z0[n]],name='Set-parall_y-10-'+str(n+1))
    index=0
    index2=0
    for j in range(len(nodeSet_parall_y_x0_z0)):
        x1=node[nodeSet_parall_y_x0_z0[j]-1].coordinates[0]
        y1=node[nodeSet_parall_y_x0_z0[j]-1].coordinates[1]
        z1=node[nodeSet_parall_y_x0_z0[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_parall_y_x0_z0[index]-1:nodeSet_parall_y_x0_z0[index]],name='Set-parall_y-0-'+str(n+1))
    m.Equation(name='Eq-P_y_10-U-'+str(n+1),terms=((1,'Set-parall_y-0-'+str(n+1),1),(1,'Ref-X',1),(-1,'Set-parall_y-10-'+str(n+1),1)))
    m.Equation(name='Eq-P_y_10-V-'+str(n+1),terms=((1,'Set-parall_y-0-'+str(n+1),2),(-1,'Set-parall_y-10-'+str(n+1),2)))
    m.Equation(name='Eq-P_y_10-W-'+str(n+1),terms=((1,'Set-parall_y-0-'+str(n+1),3),(-1,'Set-parall_y-10-'+str(n+1),3)))
    for j in range(len(nodeSet_parall_y_x1_z1)):
        x1=node[nodeSet_parall_y_x1_z1[j]-1].coordinates[0]
        y1=node[nodeSet_parall_y_x1_z1[j]-1].coordinates[1]
        z1=node[nodeSet_parall_y_x1_z1[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
            for k in range(len(nodeSet_parall_y_x0_z1)):
                x2=node[nodeSet_parall_y_x0_z1[k]-1].coordinates[0]
                y2=node[nodeSet_parall_y_x0_z1[k]-1].coordinates[1]
                z2=node[nodeSet_parall_y_x0_z1[k]-1].coordinates[2]
                distance=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                if (distance < 1.001):
                    index2=k
    r.Set(nodes=node[nodeSet_parall_y_x1_z1[index]-1:nodeSet_parall_y_x1_z1[index]],name='Set-parall_y-11-'+str(n+1))
    m.Equation(name='Eq-P_y_11-U-'+str(n+1),terms=((1,'Set-parall_y-10-'+str(n+1),1),(1,'Ref-XZ',1),(-1,'Set-parall_y-11-'+str(n+1),1)))
    m.Equation(name='Eq-P_y_11-V-'+str(n+1),terms=((1,'Set-parall_y-10-'+str(n+1),2),(1,'Ref-YZ',2),(-1,'Set-parall_y-11-'+str(n+1),2)))
    m.Equation(name='Eq-P_y_11-W-'+str(n+1),terms=((1,'Set-parall_y-10-'+str(n+1),3),(1,'Ref-Z',3),(-1,'Set-parall_y-11-'+str(n+1),3)))
    
    r.Set(nodes=node[nodeSet_parall_y_x0_z1[index2]-1:nodeSet_parall_y_x0_z1[index2]],name='Set-parall_y-01-'+str(n+1))
    m.Equation(name='Eq-P_y_01-U-'+str(n+1),terms=((1,'Set-parall_y-01-'+str(n+1),1),(-1,'Ref-XZ',1),(-1,'Set-parall_y-0-'+str(n+1),1)))
    m.Equation(name='Eq-P_y_01-V-'+str(n+1),terms=((1,'Set-parall_y-01-'+str(n+1),2),(-1,'Ref-YZ',2),(-1,'Set-parall_y-0-'+str(n+1),2)))
    m.Equation(name='Eq-P_y_01-W-'+str(n+1),terms=((1,'Set-parall_y-01-'+str(n+1),3),(-1,'Ref-Z',3),(-1,'Set-parall_y-0-'+str(n+1),3)))

# For edges parallel to the z-axis
for n in range(len(nodeSet_parall_z_x1_y0)):
    x0=node[nodeSet_parall_z_x1_y0[n]-1].coordinates[0]
    y0=node[nodeSet_parall_z_x1_y0[n]-1].coordinates[1]
    z0=node[nodeSet_parall_z_x1_y0[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_parall_z_x1_y0[n]-1:nodeSet_parall_z_x1_y0[n]],name='Set-parall_z-10-'+str(n+1))
    index=0
    index2=0
    for j in range(len(nodeSet_parall_z_x0_y0)):
        x1=node[nodeSet_parall_z_x0_y0[j]-1].coordinates[0]
        y1=node[nodeSet_parall_z_x0_y0[j]-1].coordinates[1]
        z1=node[nodeSet_parall_z_x0_y0[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_parall_z_x0_y0[index]-1:nodeSet_parall_z_x0_y0[index]],name='Set-parall_z-0-'+str(n+1))
    m.Equation(name='Eq-P_z_10-U-'+str(n+1),terms=((1,'Set-parall_z-0-'+str(n+1),1),(1,'Ref-X',1),(-1,'Set-parall_z-10-'+str(n+1),1)))
    m.Equation(name='Eq-P_z_10-V-'+str(n+1),terms=((1,'Set-parall_z-0-'+str(n+1),2),(-1,'Set-parall_z-10-'+str(n+1),2)))
    m.Equation(name='Eq-P_z_10-W-'+str(n+1),terms=((1,'Set-parall_z-0-'+str(n+1),3),(-1,'Set-parall_z-10-'+str(n+1),3)))
    for j in range(len(nodeSet_parall_z_x1_y1)):
        x1=node[nodeSet_parall_z_x1_y1[j]-1].coordinates[0]
        y1=node[nodeSet_parall_z_x1_y1[j]-1].coordinates[1]
        z1=node[nodeSet_parall_z_x1_y1[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
            for k in range(len(nodeSet_parall_z_x0_y1)):
                x2=node[nodeSet_parall_z_x0_y1[k]-1].coordinates[0]
                y2=node[nodeSet_parall_z_x0_y1[k]-1].coordinates[1]
                z2=node[nodeSet_parall_z_x0_y1[k]-1].coordinates[2]
                distance=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
                if (distance < 1.001):
                    index2=k
    r.Set(nodes=node[nodeSet_parall_z_x1_y1[index]-1:nodeSet_parall_z_x1_y1[index]],name='Set-parall_z-11-'+str(n+1))
    m.Equation(name='Eq-P_z_11-U-'+str(n+1),terms=((1,'Set-parall_z-10-'+str(n+1),1),(1,'Ref-XY',1),(-1,'Set-parall_z-11-'+str(n+1),1)))
    m.Equation(name='Eq-P_z_11-V-'+str(n+1),terms=((1,'Set-parall_z-10-'+str(n+1),2),(1,'Ref-Y',2),(-1,'Set-parall_z-11-'+str(n+1),2)))
    m.Equation(name='Eq-P_z_11-W-'+str(n+1),terms=((1,'Set-parall_z-10-'+str(n+1),3),(-1,'Set-parall_z-11-'+str(n+1),3)))
    
    r.Set(nodes=node[nodeSet_parall_z_x0_y1[index2]-1:nodeSet_parall_z_x0_y1[index2]],name='Set-parall_z-01-'+str(n+1))
    m.Equation(name='Eq-P_z_01-U-'+str(n+1),terms=((1,'Set-parall_z-01-'+str(n+1),1),(-1,'Ref-XY',1),(-1,'Set-parall_z-0-'+str(n+1),1)))
    m.Equation(name='Eq-P_z_01-V-'+str(n+1),terms=((1,'Set-parall_z-01-'+str(n+1),2),(-1,'Ref-Y',2),(-1,'Set-parall_z-0-'+str(n+1),2)))
    m.Equation(name='Eq-P_z_01-W-'+str(n+1),terms=((1,'Set-parall_z-01-'+str(n+1),3),(-1,'Set-parall_z-0-'+str(n+1),3)))

# For the vertices
for n in range(len(nodeSet_vertice_x0_y0_z0)):
    x0=node[nodeSet_vertice_x0_y0_z0[n]-1].coordinates[0]
    y0=node[nodeSet_vertice_x0_y0_z0[n]-1].coordinates[1]
    z0=node[nodeSet_vertice_x0_y0_z0[n]-1].coordinates[2]
    r.Set(nodes=node[nodeSet_vertice_x0_y0_z0[n]-1:nodeSet_vertice_x0_y0_z0[n]],name='Set-vertice_0-'+str(n+1))
    index=0
    for j in range(len(nodeSet_vertice_x1_y0_z0)):
        x1=node[nodeSet_vertice_x1_y0_z0[j]-1].coordinates[0]
        y1=node[nodeSet_vertice_x1_y0_z0[j]-1].coordinates[1]
        z1=node[nodeSet_vertice_x1_y0_z0[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_vertice_x1_y0_z0[index]-1:nodeSet_vertice_x1_y0_z0[index]],name='Set-vertice_100-'+str(n+1))
    m.Equation(name='Eq-V_100-U-'+str(n+1),terms=((1,'Set-vertice_100-'+str(n+1),1),(-1,'Ref-X',1),(-1,'Set-vertice_0-'+str(n+1),1)))
    m.Equation(name='Eq-V_100-V-'+str(n+1),terms=((1,'Set-vertice_100-'+str(n+1),2),(-1,'Set-vertice_0-'+str(n+1),2)))
    m.Equation(name='Eq-V_100-W-'+str(n+1),terms=((1,'Set-vertice_100-'+str(n+1),3),(-1,'Set-vertice_0-'+str(n+1),3)))
    for j in range(len(nodeSet_vertice_x0_y0_z1)):
        x1=node[nodeSet_vertice_x0_y0_z1[j]-1].coordinates[0]
        y1=node[nodeSet_vertice_x0_y0_z1[j]-1].coordinates[1]
        z1=node[nodeSet_vertice_x0_y0_z1[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_vertice_x0_y0_z1[index]-1:nodeSet_vertice_x0_y0_z1[index]],name='Set-vertice_001-'+str(n+1))
    m.Equation(name='Eq-V_001-U-'+str(n+1),terms=((1,'Set-vertice_0-'+str(n+1),1),(1,'Ref-XZ',1),(-1,'Set-vertice_001-'+str(n+1),1)))
    m.Equation(name='Eq-V_001-V-'+str(n+1),terms=((1,'Set-vertice_0-'+str(n+1),2),(1,'Ref-YZ',2),(-1,'Set-vertice_001-'+str(n+1),2)))
    m.Equation(name='Eq-V_001-W-'+str(n+1),terms=((1,'Set-vertice_0-'+str(n+1),3),(1,'Ref-Z',3),(-1,'Set-vertice_001-'+str(n+1),3)))
    for j in range(len(nodeSet_vertice_x0_y1_z0)):
        x1=node[nodeSet_vertice_x0_y1_z0[j]-1].coordinates[0]
        y1=node[nodeSet_vertice_x0_y1_z0[j]-1].coordinates[1]
        z1=node[nodeSet_vertice_x0_y1_z0[j]-1].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_vertice_x0_y1_z0[index]-1:nodeSet_vertice_x0_y1_z0[index]],name='Set-vertice_010-'+str(n+1))
    m.Equation(name='Eq-V_010-U-'+str(n+1),terms=((1,'Set-vertice_010-'+str(n+1),1),(-1,'Ref-XY',1),(-1,'Set-vertice_0-'+str(n+1),1)))
    m.Equation(name='Eq-V_010-V-'+str(n+1),terms=((1,'Set-vertice_010-'+str(n+1),2),(-1,'Ref-Y',2),(-1,'Set-vertice_0-'+str(n+1),2)))
    m.Equation(name='Eq-V_010-W-'+str(n+1),terms=((1,'Set-vertice_010-'+str(n+1),3),(-1,'Set-vertice_0-'+str(n+1),3)))
    for k in range(len(nodeSet_vertice_x1_y1_z0)):
        x2=node[nodeSet_vertice_x1_y1_z0[k]-1].coordinates[0]
        y2=node[nodeSet_vertice_x1_y1_z0[k]-1].coordinates[1]
        z2=node[nodeSet_vertice_x1_y1_z0[k]-1].coordinates[2]
        distance=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        if (distance < 1.001):
            index=k
    r.Set(nodes=node[nodeSet_vertice_x1_y1_z0[index]-1:nodeSet_vertice_x1_y1_z0[index]],name='Set-vertice_110-'+str(n+1))
    m.Equation(name='Eq-V_110-U-'+str(n+1),terms=((1,'Set-vertice_110-'+str(n+1),1),(-1,'Ref-X',1),(-1,'Set-vertice_010-'+str(n+1),1)))
    m.Equation(name='Eq-V_110-V-'+str(n+1),terms=((1,'Set-vertice_110-'+str(n+1),2),(-1,'Set-vertice_010-'+str(n+1),2)))
    m.Equation(name='Eq-V_110-W-'+str(n+1),terms=((1,'Set-vertice_110-'+str(n+1),3),(-1,'Set-vertice_010-'+str(n+1),3)))
    for L in range(len(nodeSet_vertice_x1_y1_z1)):
        x3=node[nodeSet_vertice_x1_y1_z1[L]-1].coordinates[0]
        y3=node[nodeSet_vertice_x1_y1_z1[L]-1].coordinates[1]
        z3=node[nodeSet_vertice_x1_y1_z1[L]-1].coordinates[2]
        distance=sqrt((x3-x2)**2+(y3-y2)**2+(z3-z2)**2)
        if (distance < 1.001):
            index=L
    r.Set(nodes=node[nodeSet_vertice_x1_y1_z1[index]-1:nodeSet_vertice_x1_y1_z1[index]],name='Set-vertice_111-'+str(n+1))
    m.Equation(name='Eq-V_111-U-'+str(n+1),terms=((1,'Set-vertice_111-'+str(n+1),1),(-1,'Ref-XZ',1),(-1,'Set-vertice_110-'+str(n+1),1)))
    m.Equation(name='Eq-V_111-V-'+str(n+1),terms=((1,'Set-vertice_111-'+str(n+1),2),(-1,'Ref-YZ',2),(-1,'Set-vertice_110-'+str(n+1),2)))
    m.Equation(name='Eq-V_111-W-'+str(n+1),terms=((1,'Set-vertice_111-'+str(n+1),3),(-1,'Ref-Z',3),(-1,'Set-vertice_110-'+str(n+1),3)))
    for k in range(len(nodeSet_vertice_x1_y0_z1)):
        x2=node[nodeSet_vertice_x1_y0_z1[k]-1].coordinates[0]
        y2=node[nodeSet_vertice_x1_y0_z1[k]-1].coordinates[1]
        z2=node[nodeSet_vertice_x1_y0_z1[k]-1].coordinates[2]
        distance=sqrt((x2-x3)**2+(y2-y3)**2+(z2-z3)**2)
        if (distance < 1.001):
            index=k
    r.Set(nodes=node[nodeSet_vertice_x1_y0_z1[index]-1:nodeSet_vertice_x1_y0_z1[index]],name='Set-vertice_101-'+str(n+1))
    m.Equation(name='Eq-V_101-U-'+str(n+1),terms=((1,'Set-vertice_101-'+str(n+1),1),(1,'Ref-XY',1),(-1,'Set-vertice_111-'+str(n+1),1)))
    m.Equation(name='Eq-V_101-V-'+str(n+1),terms=((1,'Set-vertice_101-'+str(n+1),2),(1,'Ref-Y',2),(-1,'Set-vertice_111-'+str(n+1),2)))
    m.Equation(name='Eq-V_101-W-'+str(n+1),terms=((1,'Set-vertice_101-'+str(n+1),3),(-1,'Set-vertice_111-'+str(n+1),3)))
    for j in range(len(nodeSet_vertice_x0_y1_z1)):
        x1=node[nodeSet_vertice_x0_y1_z1[j]-1].coordinates[0]
        y1=node[nodeSet_vertice_x0_y1_z1[j]-1].coordinates[1]
        z1=node[nodeSet_vertice_x0_y1_z1[j]-1].coordinates[2]
        distance=sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        if (distance < 1.001):
            index=j
    r.Set(nodes=node[nodeSet_vertice_x0_y1_z1[index]-1:nodeSet_vertice_x0_y1_z1[index]],name='Set-vertice_011-'+str(n+1))
    m.Equation(name='Eq-V_011-U-'+str(n+1),terms=((1,'Set-vertice_011-'+str(n+1),1),(1,'Ref-X',1),(-1,'Set-vertice_111-'+str(n+1),1)))
    m.Equation(name='Eq-V_011-V-'+str(n+1),terms=((1,'Set-vertice_011-'+str(n+1),2),(-1,'Set-vertice_111-'+str(n+1),2)))
    m.Equation(name='Eq-V_011-W-'+str(n+1),terms=((1,'Set-vertice_011-'+str(n+1),3),(-1,'Set-vertice_111-'+str(n+1),3)))






r.SetFromNodeLabels(name='nodeSet_X_minus', nodeLabels=((insName,nodeSet_X_minus),))
r.SetFromNodeLabels(name='nodeSet_X_plus', nodeLabels=((insName,nodeSet_X_plus),))
r.SetFromNodeLabels(name='nodeSet_Y_minus', nodeLabels=((insName,nodeSet_Y_minus),))
r.SetFromNodeLabels(name='nodeSet_Y_plus', nodeLabels=((insName,nodeSet_Y_plus),))
r.SetFromNodeLabels(name='nodeSet_Z_minus', nodeLabels=((insName,nodeSet_Z_minus),))
r.SetFromNodeLabels(name='nodeSet_Z_plus', nodeLabels=((insName,nodeSet_Z_plus),))

r.SetFromNodeLabels(name='nodeSet_vertice_x0_y0_z0', nodeLabels=((insName,nodeSet_vertice_x0_y0_z0),))
r.SetFromNodeLabels(name='nodeSet_vertice_x1_y0_z0', nodeLabels=((insName,nodeSet_vertice_x1_y0_z0),))
r.SetFromNodeLabels(name='nodeSet_vertice_x0_y1_z0', nodeLabels=((insName,nodeSet_vertice_x0_y1_z0),))
r.SetFromNodeLabels(name='nodeSet_vertice_x0_y0_z1', nodeLabels=((insName,nodeSet_vertice_x0_y0_z1),))
r.SetFromNodeLabels(name='nodeSet_vertice_x1_y1_z0', nodeLabels=((insName,nodeSet_vertice_x1_y1_z0),))
r.SetFromNodeLabels(name='nodeSet_vertice_x0_y1_z1', nodeLabels=((insName,nodeSet_vertice_x0_y1_z1),))
r.SetFromNodeLabels(name='nodeSet_vertice_x1_y0_z1', nodeLabels=((insName,nodeSet_vertice_x1_y0_z1),))
r.SetFromNodeLabels(name='nodeSet_vertice_x1_y1_z1', nodeLabels=((insName,nodeSet_vertice_x1_y1_z1),))

r.SetFromNodeLabels(name='nodeSet_parall_x_y0_z0' , nodeLabels=((insName,nodeSet_parall_x_y0_z0),))
r.SetFromNodeLabels(name='nodeSet_parall_x_y0_z1' , nodeLabels=((insName,nodeSet_parall_x_y0_z1),))
r.SetFromNodeLabels(name='nodeSet_parall_x_y1_z0' , nodeLabels=((insName,nodeSet_parall_x_y1_z0),))
r.SetFromNodeLabels(name='nodeSet_parall_x_y1_z1' , nodeLabels=((insName,nodeSet_parall_x_y1_z1),))
r.SetFromNodeLabels(name='nodeSet_parall_y_x0_z0' , nodeLabels=((insName,nodeSet_parall_y_x0_z0),))
r.SetFromNodeLabels(name='nodeSet_parall_y_x0_z1' , nodeLabels=((insName,nodeSet_parall_y_x0_z1),))
r.SetFromNodeLabels(name='nodeSet_parall_y_x1_z0' , nodeLabels=((insName,nodeSet_parall_y_x1_z0),))
r.SetFromNodeLabels(name='nodeSet_parall_y_x1_z1' , nodeLabels=((insName,nodeSet_parall_y_x1_z1),))
r.SetFromNodeLabels(name='nodeSet_parall_z_x0_y0' , nodeLabels=((insName,nodeSet_parall_z_x0_y0),))
r.SetFromNodeLabels(name='nodeSet_parall_z_x0_y1' , nodeLabels=((insName,nodeSet_parall_z_x0_y1),))
r.SetFromNodeLabels(name='nodeSet_parall_z_x1_y0' , nodeLabels=((insName,nodeSet_parall_z_x1_y0),))
r.SetFromNodeLabels(name='nodeSet_parall_z_x1_y1' , nodeLabels=((insName,nodeSet_parall_z_x1_y1),))

