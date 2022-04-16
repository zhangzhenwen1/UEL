# 3D PBC  
from abaqus import *
from abaqusConstants import *
from caeModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
from interaction import * 



fields=(('x min:','0.0'),('x max:','1.0'),('y min:','0.0'),('y max:','1.0'),('z min:','0.0'),('z max:','1.0'))
xmin,xmax,ymin,ymax,zmin,zmax=getInputs(fields=fields,label="Input the parameters")
xmin=float(xmin);xmax=float(xmax);ymin=float(ymin);ymax=float(ymax);zmin=float(zmin);zmax=float(zmax)

tor=0.001

m=mdb.models['RVE-test1']
r=m.rootAssembly

refpoints=r.referencePoints
for i in range(len(refpoints)):
    r.Set(name='Ref-'+str(i+1),referencePoints=(refpoints.values()[i],))
    
    
node=r.instances['PART-1-1'].nodes
ne1=[]
ne2=[]
ne3=[]
ne4=[] 
ne5=[]
ne6=[]   
for i in range(len(node)):
    x=node[i].coordinates[0]
    y=node[i].coordinates[1]
    z=node[i].coordinates[2]
    if (abs(x-xmin)<tor):
        ne1.append(i)
    if (abs(x-xmax)<tor):
        ne2.append(i)
    if (abs(y-ymin)<tor):
        ne3.append(i)
    if (abs(y-ymax)<tor):
        ne4.append(i) 
    if (abs(z-zmin)<tor):
        ne5.append(i)
    if (abs(z-zmax)<tor):
        ne6.append(i)   


for n in range(min(len(ne1),len(ne2))):
    x0=node[ne1[n]].coordinates[0]
    y0=node[ne1[n]].coordinates[1]
    z0=node[ne1[n]].coordinates[2]
    r.Set(nodes=node[ne1[n]:ne1[n]+1],name='Set-L-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(ne2)):
        x1=node[ne2[j]].coordinates[0]
        y1=node[ne2[j]].coordinates[1]
        z1=node[ne2[j]].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j
    r.Set(nodes=node[ne2[index]:ne2[index]+1],name='Set-R-'+str(n+1))
    m.Equation(name='Eq-LR-X-'+str(n+1),terms=((1,'Set-L-'+str(n+1),1),(-1,'Ref-1',1),(-1,'Set-R-'+str(n+1),1)))
    m.Equation(name='Eq-LR-Y-'+str(n+1),terms=((1,'Set-L-'+str(n+1),2),(-1,'Ref-1',2),(-1,'Set-R-'+str(n+1),2)))
    m.Equation(name='Eq-LR-Z-'+str(n+1),terms=((1,'Set-L-'+str(n+1),3),(-1,'Ref-1',3),(-1,'Set-R-'+str(n+1),3)))


for n in range(min(len(ne3),len(ne4))):
    x0=node[ne3[n]].coordinates[0] 
    y0=node[ne3[n]].coordinates[1]
    z0=node[ne3[n]].coordinates[2]
    r.Set(nodes=node[ne3[n]:ne3[n]+1],name='Set-B-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(ne4)):
        x1=node[ne4[j]].coordinates[0]
        y1=node[ne4[j]].coordinates[1]
        z1=node[ne4[j]].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j                                                                                                                                                                                                          
    r.Set(nodes=node[ne4[index]:ne4[index]+1],name='Set-T-'+str(n+1))
    m.Equation(name='Eq-BT-X-'+str(n+1),terms=((1,'Set-B-'+str(n+1),1),(-1,'Ref-2',1),(-1,'Set-T-'+str(n+1),1)))
    m.Equation(name='Eq-BT-Y-'+str(n+1),terms=((1,'Set-B-'+str(n+1),2),(-1,'Ref-2',2),(-1,'Set-T-'+str(n+1),2)))
    m.Equation(name='Eq-BT-Z-'+str(n+1),terms=((1,'Set-B-'+str(n+1),3),(-1,'Ref-2',3),(-1,'Set-T-'+str(n+1),3)))
    

for n in range(min(len(ne5),len(ne6))):
    x0=node[ne5[n]].coordinates[0] 
    y0=node[ne5[n]].coordinates[1]
    z0=node[ne5[n]].coordinates[2]
    r.Set(nodes=node[ne5[n]:ne5[n]+1],name='Set-F-'+str(n+1))
    mindistance=sqrt((xmin-xmax)**2+(ymin-ymax)**2+(zmin-zmax)**2)*2
    index=0
    for j in range(len(ne6)):
        x1=node[ne6[j]].coordinates[0]
        y1=node[ne6[j]].coordinates[1]
        z1=node[ne6[j]].coordinates[2]
        distance=sqrt((x0-x1)**2+(y0-y1)**2+(z0-z1)**2)
        if (distance < mindistance):
            mindistance=distance
            index=j                                                                                                                                                                                                          
    r.Set(nodes=node[ne6[index]:ne6[index]+1],name='Set-N-'+str(n+1))
    m.Equation(name='Eq-FN-X-'+str(n+1),terms=((1,'Set-F-'+str(n+1),1),(-1,'Ref-3',1),(-1,'Set-N-'+str(n+1),1)))
    m.Equation(name='Eq-FN-Y-'+str(n+1),terms=((1,'Set-F-'+str(n+1),2),(-1,'Ref-3',2),(-1,'Set-N-'+str(n+1),2)))
    m.Equation(name='Eq-FN-Z-'+str(n+1),terms=((1,'Set-F-'+str(n+1),3),(-1,'Ref-3',3),(-1,'Set-N-'+str(n+1),3))) 