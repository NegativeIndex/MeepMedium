import numpy as np
import meep as mp
import argparse
import math
import cmath
import sys
import datetime
import re

##############################
# write data into a file
##############################
# now if fname is empty, it write to screen
# default to file and screen both
def matrix_output(fname,data,fmt,sig=None,onScreen=True):
    n1=data.shape[0]
    n2=data.shape[1]
    ss=''
 
    for i in range(n1):
        line=''
        if sig is not None:
            line+=sig+":"
        for j in range(n2):
            if line:
                line+=", "
            line+=fmt.format(data[i,j])
        if i<n1-1:
            ss+=line+"\n"
        else:
            ss+=line
  
    if fname:
        with open(fname, "w") as fout:
            print(ss,file=fout)
            if onScreen:
                print(ss)
    else:
        print(ss)

##############################
# organize data for saving
##############################
# def file_data(*argv):  
#     data=np.array()
#     for arg in argv: 
#         data=np.column_stack(( data,np.array(arg) ))
#     return data

##############################
# a function
##############################
def set_polarization(pol):
    myMap = {}
    myMap["Ex"] =mp.Ex
    myMap["Ey"] =mp.Ey
    myMap["Ez"] =mp.Ez
    myMap["Hx"] =mp.Hx
    myMap["Hy"] =mp.Hy
    myMap["Hz"] =mp.Hz
    return myMap[pol]

##############################
# define flux region
##############################
def FluxRegion_box(x1,x2,y1,y2,z1,z2,loc):
    # giving the coordinate of the points, return flux region
    fx=abs(x1-x2)
    fy=abs(y1-y2)
    fz=abs(z1-z2)
    
    cx=(x1+x2)/2
    cy=(y1+y2)/2
    cz=(z1+z2)/2
    
    if loc=="X+":
        center=mp.Vector3(max(x1,x2),cy,cz)
        size=mp.Vector3(0,fy,fz)
        weight=1
    elif loc=="X-":
        center=mp.Vector3(min(x1,x2),cy,cz)
        size=mp.Vector3(0,fy,fz)
        weight=-1
    elif loc=="Y+":
        center=mp.Vector3(cx,max(y1,y2),cz)
        size=mp.Vector3(fx,0,fz)
        weight=1
    elif loc=="Y-":
        center=mp.Vector3(cx,min(y1,y2),cz)
        size=mp.Vector3(fx,0,fz)
        weight=-1
    elif loc=="Z+":
        center=mp.Vector3(cx,cy,max(z1,z2))
        size=mp.Vector3(fx,fy,0)  
        weight=1     
    elif loc=="Z-":
        center=mp.Vector3(cx,cy,min(z1,z2))
        size=mp.Vector3(fx,fy,0)    
        weight=-1 
    else:
        sys.exit("Flux region setting is wrong")

    return mp.FluxRegion(center=center,size=size,weight=weight)

###############################################
# define flux region from box length and center
###############################################
def FluxRegion_box_Center(x,y,z,c,loc):
    x1=c.x-abs(x)/2
    x2=c.x+abs(x)/2
    y1=c.y-abs(y)/2
    y2=c.y+abs(y)/2
    z1=c.z-abs(z)/2
    z2=c.z+abs(z)/2
    return FluxRegion_box(x1,x2,y1,y2,z1,z2,loc)


###############################################
# return an array of 6 flux regions from corner
###############################################
def FluxRegions_Box_Corner(x1,x2,y1,y2,z1,z2):
    frx1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"X+")
    frx2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"X-")
    fry1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Y+")
    fry2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Y-")
    frz1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Z+")
    frz2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Z-")
    return [frx1,frx2,fry1,fry2,frz1,frz2]

###############################################
# return an array of 6 flux regions from center
###############################################
def FluxRegions_Box_Center(x,y,z,center):
    c=center
    x1=c.x-abs(x)/2
    x2=c.x+abs(x)/2
    y1=c.y-abs(y)/2
    y2=c.y+abs(y)/2
    z1=c.z-abs(z)/2
    z2=c.z+abs(z)/2
    frx1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"X+")
    frx2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"X-")
    fry1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Y+")
    fry2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Y-")
    frz1=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Z+")
    frz2=FluxRegion_box(x1,x2,y1,y2,z1,z2,"Z-")
    return [frx1,frx2,fry1,fry2,frz1,frz2]

def FluxRegions_Cube_Center(l,center):
    return FluxRegions_Box_Center(l,l,l,center)


###############################################
def FluxRegions_Box2D_Center(lx,ly,center):
    x1=center.x-abs(lx)/2
    x2=center.x+abs(lx)/2
    y1=center.y-abs(ly)/2
    y2=center.y+abs(ly)/2

    frx1=mp.FluxRegion(center=mp.Vector3(x1,(y1+y2)/2,0),
                      size=mp.Vector3(0,abs(y2-y1),0))
    frx2=mp.FluxRegion(center=mp.Vector3(x2,(y1+y2)/2,0),
                      size=mp.Vector3(0,abs(y2-y1),0))

    fry1=mp.FluxRegion(center=mp.Vector3((x1+x2)/2,y1,0),
                      size=mp.Vector3(abs(x2-x1),0,0))
    fry2=mp.FluxRegion(center=mp.Vector3((x1+x2)/2,y2,0),
                       size=mp.Vector3(abs(x2-x1),0,0))

    return [frx1,frx2,fry1,fry2]

###############################################
# calculate flux box center and size
###############################################
def FluxBox2D_helper(size,xsrc,xmin=None,xmax=None):
    if xmin is None:
        xmin=xsrc-10*size
    if xmax is None:
        xmax=xsrc+10*size

    if (xmax-xmin)<2*size:
        lx=(xmax-xmin)/2
        xbox=xsrc/2+xmin/4+xmax/4
    else:
        lx=size
        if (xsrc-xmin)<lx:
            xbox=(xsrc+xmin)/2+lx/2
        elif (xmax-xsrc)<lx:
            xbox=(xsrc+xmax)/2-lx/2
        else:
            xbox=xsrc
    return (lx,xbox)


def FluxBox2D(size,src,xmin=None,xmax=None,ymin=None,ymax=None):
    # size is the max size of the box
    # src is a meep vector3
    # xmin,xmax,ymin,ymax are boundaries
    # the idea is the box should be close to the src than boundary
    lx,xbox=FluxBox2D_helper(size,src.x,xmin,xmax)
    ly,ybox=FluxBox2D_helper(size,src.y,ymin,ymax)
    return (min(lx,ly),mp.Vector3(xbox,ybox,0))
    
def FluxBox3D(size,src,xmin=None,xmax=None,
              ymin=None,ymax=None,
              zmin=None,zmax=None):
    # size is the max size of the box
    # src is a meep vector3
    # xmin,xmax,ymin,ymax are boundaries
    # the idea is the box should be close to the src than boundary
    lx,xbox=FluxBox2D_helper(size,src.x,xmin,xmax)
    ly,ybox=FluxBox2D_helper(size,src.y,ymin,ymax)
    lz,zbox=FluxBox2D_helper(size,src.z,zmin,zmax)

    return (min(lx,ly,lz),mp.Vector3(xbox,ybox,zbox))


###############################################
# force meep to demonstrate information on time
###############################################
def my_flush_step(sim):
    t=datetime.datetime.now()
    print("My flush: "+str(t))
    sys.stdout.flush()


###############################################
# convert number of second into a readable string
###############################################
def nice_sec2str(x):
    sec=datetime.timedelta(seconds=x)
    d = datetime.datetime(1,1,1) + sec
    unit=["days","hours","minutes","seconds"]
    if d.day<3:
        unit[0]="day"
    if d.hour<2:
        unit[1]="hour"
    if d.minute<2:
        unit[2]="minute"
    if d.second<2:
        unit[3]="second"

    if (d.day>1):
        ss="{:d} {} {:d} {} {:d} {} {:d} {} ".format(
            d.day-1,unit[0],d.hour,unit[1], 
            d.minute,unit[2], d.second,unit[3])
    elif (d.hour>0):
        ss="{:d} {} {:d} {} {:d} {} ".format(
            d.hour,unit[1], 
            d.minute,unit[2], d.second,unit[3])
    elif (d.minute>0):
        ss="{:d} {} {:d} {} ".format(
            d.minute,unit[2], d.second,unit[3])
    else:
       ss="{:d} {} ".format(d.second,unit[3])
 
    return ss


##############################
# replace strings in a file
##############################
def file_re_sub(fname1,fname2,*argv):  
    with open(fname1, "rt") as fin:
        fin=open(fname1, 'rt')
        lines=fin.readlines()
        fin.close()

        for i,line in enumerate(lines):
            for pair in argv:
                line=re.sub(pair[0],pair[1],line)
            lines[i]=line
    
        fout=open(fname2, 'wt')
        fout.write("".join(lines))
        fout.close()



##############################
# add a line to a file
##############################
def file_add_line(fname1,fname2,newline='',
                before=None,after=None,number=None):
    # before and after are regular expression
    # number is an interger
    # if before,after and number are not None, 
    # all the condition must be satisfied
    with open(fname1, "rt") as fin:
        fin=open(fname1, 'rt')
        lines=fin.readlines()
        fin.close()
    
        # find before or after, the first match
        number1=None
        number2=None
        for i,line in enumerate(lines):
            if before is not None and re.search(before, line):
                if number1 is None:
                    number1=i
            if after is not None and re.search(after, line):
                if number2 is None:
                    number2=i
        # print(number1)
        # print(number2)
        # check condition
        idx=number
        if idx is not None:
            if number1 is not None and (number1-idx)!=0:
                idx=None
        elif number1 is not None:
            idx=number1

        if idx is not None:
            if number2 is not None and (number2-idx)!=-1:
                idx=None
        elif number2 is not None:
            idx=number2+1
            
        # write the line
        if idx!=None:
            lines.insert(idx,newline+'\n')
       
        fout=open(fname2, 'wt')
        fout.write("".join(lines))
        fout.close()


##############################
# calculate an array of centers based on thicknesses 
##############################
def center_from_thickness(ds,x0):
    ds=np.array(ds)
    return np.array([np.sum(ds[0:i+1]) for i in range(ds.size)])-ds/2+x0
