#! /usr/bin/env python3

import numpy as np
import sys, getopt
from generic import obj   
import math

atoms=['','H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn']

def replace_zxe(inputfile):

    print('\nUpdating ZXE...')
    f=open(inputfile, "r")
    
    d=np.array(f.read().split())
    
    f.seek(0)
    pptext = f.read()
    f.close()

    # HEADER
    dhead=d[0:2]
    # CHANNEL LENGTHS
    nchan=int(dhead[1])
    dchan=d[2:2+nchan]
    dchan=np.array(dchan,dtype=int)
    # POTENTIALS
    dpots=[]
    pos=0
    for i in range(len(dchan)):
        dtmp = d[2+nchan+int(pos):2+nchan+int(pos+dchan[i]*3)]
        dtmp.shape=len(dtmp)//3,3
        pos+=dchan[i]*3
        dpots.append(dtmp)
    #end for
    
    for t in dpots[2]:
        if t[0]=='1':
            pptext = pptext.replace('ZXE',str(float(t[1])*float(t[2])))
    
    f=open(inputfile, "w")
    f.write(pptext)

#end def


def add_to_vL2(pots,pwrs,exps,coeffs):

    tpots = pots.copy()

    lmax = len(tpots)-1
    for l in range(len(tpots)):
        if l<lmax:
            fctr = l*(l+1)-lmax*(lmax+1)
        else:
            fctr = lmax*(lmax+1)
        #end if
        for t in range(len(pwrs)):
            tpots[l] = np.append(tpots[l],np.array([[pwrs[t],exps[t],coeffs[t]*fctr]]),axis=0)
        #end for
    #end for
    return tpots

#end def


def transform_to_L2(pots,keep="s p",lmax=None):

    print('transform_to_L2')

    chan_labels = ['s','p','d','f','g','h','i','j']
    tpots = pots.copy() 
    
    keep_chans = keep.split()

    # Are the labels recognized?
    if keep_chans[0] not in chan_labels or keep_chans[1] not in chan_labels:
        print('Requested channel to keep is not recognized')
        quit()
    #end if

    # Does the original potential contain the requested channels?
    if chan_labels.index(keep_chans[0]) > len(tpots)-1 or chan_labels.index(keep_chans[1]) > len(tpots)-1:
        print('Cannot keep channel that is not already present')
        quit()
    #end if

    # Is lmax atleast a large as the original local
    if lmax<len(tpots)-1:
        print('Lmax is less than local of original potential')
        quit()
    #end if

    # Are the requested channels different?
    if chan_labels.index(keep_chans[0]) == chan_labels.index(keep_chans[1]):
        print('The two channels must be different from eatch other.')
        quit()
    #end if

    keep_l_vals = []
    keep_l_vals.append(chan_labels.index(keep_chans[0]))
    keep_l_vals.append(chan_labels.index(keep_chans[1]))
    keep_l_vals.sort()

    # Is one of the channels the local channel?
    if keep_l_vals[0]==len(tpots)-1 or keep_l_vals[1]==len(tpots)-1:
        keep_local=True
    else:
        keep_local=False
    #end if

    if not keep_local:
    
        vm = tpots[keep_l_vals[0]].copy()
        lm = keep_l_vals[0]
        vn = tpots[keep_l_vals[1]].copy()
        ln = keep_l_vals[1]

        vlocp = tpots[-1].copy()

        denom = lm*(lm+1)-ln*(ln+1)

        fctr1 = lm*(lm+1)/denom
        for i,p in enumerate(vn):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr1]]),axis=0)
        #end for
        fctr2 = -ln*(ln+1)/denom
        for i,p in enumerate(vm):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr2]]),axis=0)
        #end for
        fctr3 = lmax*(lmax+1)/denom
        for i,p in enumerate(vm):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr3]]),axis=0)
        #end for
        fctr4 = -lmax*(lmax+1)/denom
        for i,p in enumerate(vn):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr4]]),axis=0)
        #end for

        rpots = [] 
        for l in range(lmax):

            vtmp = vm.copy()

            fctr1 = l*(l+1)/denom
            for i,p in enumerate(vm):
                vtmp[i][2] = p[2]*fctr1
            #end for
            fctr2 = -l*(l+1)/denom
            for i,p in enumerate(vn):
                vtmp = np.append(vtmp,np.array([[p[0],p[1],p[2]*fctr2]]),axis=0)
            #end for
            fctr3 = -lmax*(lmax+1)/denom
            for i,p in enumerate(vm):
                vtmp = np.append(vtmp,np.array([[p[0],p[1],p[2]*fctr3]]),axis=0)
            #end for
            fctr4 = lmax*(lmax+1)/denom
            for i,p in enumerate(vn):
                vtmp = np.append(vtmp,np.array([[p[0],p[1],p[2]*fctr4]]),axis=0)
            #end for
            rpots.append(vtmp)

        #end for
        rpots.append(vlocp)

    else:
        
        print('KEEP LOCAL')
        # Channel selected to remain unmodified
        vm = tpots[keep_l_vals[0]].copy()
        lm = keep_l_vals[0]
      
        # Local channel. Selected to remain unmodified
        vlocp = tpots[keep_l_vals[1]].copy()
        lloc = keep_l_vals[1]


        fctr1 = (-lloc*(lloc+1))/(lm*(lm+1)-lloc*(lloc+1))
        for i,p in enumerate(vm):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr1]]),axis=0)
        #end for
        fctr2 = (lmax*(lmax+1))/(lm*(lm+1)-lloc*(lloc+1))
        for i,p in enumerate(vm):
            vlocp = np.append(vlocp,np.array([[p[0],p[1],p[2]*fctr2]]),axis=0)
        #end for

        rpots = [] 
        for l in range(lmax):

            vtmp = vm.copy()

            fctr1 = (l*(l+1))/(lm*(lm+1)-lloc*(lloc+1))
            for i,p in enumerate(vm):
                vtmp[i][2] = p[2]*fctr1
            #end for
            fctr2 = -(lmax*(lmax+1))/(lm*(lm+1)-lloc*(lloc+1))
            for i,p in enumerate(vm):
                vtmp = np.append(vtmp,np.array([[p[0],p[1],p[2]*fctr2]]),axis=0)
            #end for
            rpots.append(vtmp)

        #end for
        rpots.append(vlocp)
    #end if

    return rpots

#end def

def trans1(pots,lmax):

    tpots = pots.copy() 

    tpots[1] = tpots[0].copy()
    for i,p in enumerate(tpots[1]):
        tpots[1][i][2] = 2.*p[2]/3.

    return tpots

#end def

def trans2(pots,lmax):

    tpots = pots.copy() 

    tpots[0] = tpots[1].copy()
    for i,p in enumerate(tpots[0]):
        tpots[0][i][2] = 3.*p[2]/2.

    return tpots

#end def

def trans3(pots,lmax):

    tpots = pots.copy() 

    vs = tpots[0].copy()
    for i,term in enumerate(vs):
        vs[i][2] *= -2.

    vp = tpots[1].copy()
    for i,term in enumerate(vp):
        vp[i][2] *= 3.

    tpots[2] = np.append(tpots[2],vs,axis=0)
    tpots[2] = np.append(tpots[2],vp,axis=0)

    for i,term in enumerate(vs):
        vs[i][2] *= -1.

    for i,term in enumerate(vp):
        vp[i][2] *= -1.

    tpots[0] = np.append(tpots[0],vs,axis=0)
    tpots[0] = np.append(tpots[0],vp,axis=0)

    tpots[1] = np.append(tpots[1],vs,axis=0)
    tpots[1] = np.append(tpots[1],vp,axis=0)

    return tpots

#end def

def printchansyoon(pots,Z,outfile=None):
    if outfile is None:
        print('\nPrinting pseudopotential in yoon format...')
        print(int(Z),len(pots))
        for p in pots:
            print(len(p),end=' ')
        print()
        for i in range(len(pots)):
            for j,p in enumerate(pots[i]):
                print("{:d} {:12f} {:12f}".format(int(p[0]),p[1],p[2]))
            #end for
        #end for
    else:
        f = open(outfile, "w")
        f.write(str(int(Z))+' '+str(len(pots))+'\n')
        for p in pots:
            f.write(str(len(p))+' ')
        f.write('\n')
        for i in range(len(pots)):
            for j,p in enumerate(pots[i]):
                f.write("{:d} {} {} \n".format(int(p[0]),p[1],p[2]))
            #end for
        #end for
        f.close()
#end def

def printchansmolpro(pots,Z,zae=None,outfile=None):

    # EXAMPLE
    # ecp,Co,10,2,0
    # 4
    # 1, 25.00124115981202,  17.0
    # 3, 22.83490096546710,  425.0210997168043400
    # 2, 23.47468155934268,  -195.48211282934741
    # 2, 10.33794825313753,  -2.81572866482791
    # 2
    # 2, 23.41427030715358,  271.77708486766420
    # 2, 10.76931694175074,  54.26461121615731
    # 2
    # 2, 25.47446316631093,  201.53430745305457
    # 2, 10.68404901015315,  38.99231927238777

    if zae is None:
        print('zae cannot be none')
        quit()

    if outfile is None:
        print('\nPrinting pseudopotential in molpro format...')
        idxs = np.arange(len(pots)-1)
        idxs = np.insert(idxs,0,len(pots)-1)
        print('ecp,{},{},{},0'.format(atoms[zae],int(zae-Z),int(len(pots)-1)))
        for i in idxs:
            print('{};'.format(len(pots[i])))
            for j,p in enumerate(pots[i]):
                print("{:d}, {}, {};".format(int(p[0]),p[1],p[2]))
            #end for
        #end for
    else:
        f = open(outfile, "w")
        idxs = np.arange(len(pots)-1)
        idxs = np.insert(idxs,0,len(pots)-1)
        f.write('ecp,{},{},{},0\n'.format(atoms[zae],int(zae-Z),int(len(pots)-1)))
        for i in idxs:
            f.write('{};\n'.format(len(pots[i])))
            for j,p in enumerate(pots[i]):
                f.write("{:d}, {}, {};\n".format(int(p[0]),p[1],p[2]))
            #end for
        #end for
        f.close()
    #end if

#end def

def transform_yoon(filepath,keep,lmax,ncore=None,outfile=None,form='yoon'):

    f=open(filepath, "r")
    
    d=np.array(f.read().split(),dtype=float)
    f.close()
    
    # HEADER
    dhead=d[0:2]
    # CHANNEL LENGTHS
    nchan=int(dhead[1])
    dchan=d[2:2+nchan]
    # POTENTIALS
    dpots=[]
    pos=0
    for i in range(len(dchan)):
        dtmp = d[2+nchan+int(pos):2+nchan+int(pos+dchan[i]*3)]
        dtmp.shape=len(dtmp)//3,3
        pos+=dchan[i]*3
        dpots.append(dtmp)
    #end for

    if keep is not None:
        tpots = transform_to_L2(dpots,keep,lmax)
    #end if

    print('HERE')
    if form=='yoon':
        printchansyoon(tpots,dhead[0],outfile=outfile)
    elif form=='molpro':
        if ncore is None:
            print('For molpro format, the number of core electrons must be given via "ncore"')
            quit()
        #end if
        printchansmolpro(tpots,dhead[0],ncore,outfile=outfile)
    else:
        print('form not recognized')
    #end if

    print('HERE')
    return tpots
        
#end def


def transform_molpro(filepath,transtype,zae=None,outfile=None,maxl=None):

    if zae is None:
        print('zae cannot be none')
        quit()

    f=open(filepath, "r")
    
    d=np.array(f.read().split(),dtype=float)
    f.close()
    
    # HEADER
    dhead=d[0:2]
    # CHANNEL LENGTHS
    nchan=int(dhead[1])
    dchan=d[2:2+nchan]
    # POTENTIALS
    dpots=[]
    pos=0
    for i in range(len(dchan)):
        dtmp = d[2+nchan+int(pos):2+nchan+int(pos+dchan[i]*3)]
        dtmp.shape=len(dtmp)//3,3
        pos+=dchan[i]*3
        dpots.append(dtmp)
    #end for
    
    if transtype is not None:
    
        if transtype==0:
            tpots = dpots.copy()
        elif transtype==1:
            tpots = trans1(dpots)
        elif transtype==2:
            tpots = trans2(dpots)
        elif transtype==3:
            tpots = trans3(dpots)
        #end if

    #end if

    if maxl is not None:


        if maxl<len(tpots):
            print('MAXL must be larger than largest angular channel already present')
            quit()
        #end if

        max_pots = tpots.copy() 

        # First, construct local channel (corresponding to maxl)

        # Copy s
        vlmax = tpots[0].copy() 
        for i,p in enumerate(vlmax):
            vlmax[i][2] = (1.0-maxl*(maxl+1.0)/6.0)*p[2]
        #end for


        # Re-construct all channels. vlmax needs to be subtracted away from each non-local channel
        # new_vs = old_vs - vlmax

        for i in range(len(tpots)):
            max_pots[i]=tpots[0].copy()
            for j,p in enumerate(max_pots[i]):
                max_pots[i][j][2] = ((1.0-i*(i+1)/6.0) - (1.0-maxl*(maxl+1.0)/6.0))*p[2]
            #end for
        #end for

        # Append any additional non-local channels

        for i in range(len(tpots),maxl):
            new_chan = tpots[0].copy()
            for j,p in enumerate(new_chan):
                new_chan[j][2] = ((1.0-i*(i+1)/6.0) - (1.0-maxl*(maxl+1.0)/6.0))*p[2]
            #end for
            max_pots.append(new_chan)
        #end for

        # Now construct new local channel
        new_local = tpots[len(tpots)-1].copy()
        for p in vlmax:
            new_local = np.vstack([new_local,p])
        max_pots.append(new_local)

        tpots=max_pots.copy()
    
    printchansmolpro(tpots,dhead[0],zae=zae,outfile=outfile)
        
#end def


def get_yoon_eigs(out,valence=False):

    at_eigs=False
    f = open(out,'r')
    eigs=[]
    lvals = []
    for line in f:
        words = line.split()
        if len(words)==7:
            at_eigs=True
            eigs.append(words[4])
            lvals.append(words[1])
        else:
            if at_eigs:
                break
            #end if
        #end if
    #end for
    eigs=np.array(eigs,dtype=float)
    lvals=np.array(lvals,dtype=int)

    if valence:
        eigs_tmp=[0]*(max(lvals)+1)
        for i in range(len(eigs)):
            eigs_tmp[lvals[i]]=eigs[i]
        eigs_tmp=np.array(eigs_tmp,dtype=float)
        eigs=eigs_tmp.copy()

    return eigs

#end def


def printchans(pots):
    print('\nPrinting pseudopotential...')
    clabels=['S','P','D','F','G']
    lloc = clabels[len(pots)-1]
    for i in range(len(pots)):
        if i==len(pots)-1:
            print(lloc)
        else:
            print(str(clabels[i])+'-'+str(lloc))
        #end if
        for i,p in enumerate(pots[i]):
            print("{:3d} {:12f} {:12f}".format(int(p[0]),p[1],p[2]))
        #end for
    #end for
#end def

def ppchannel(r,coeff,alpha,powers):
    val = 0
    for ci,c in enumerate(coeff):
        val = val+c*r**(powers[ci]-2)*np.exp(-alpha[ci]*r*r)
    #end for
    return val
#end def

def p7(x,c):
    val=0
    for ci,cv in enumerate(c):
        val+=cv*x**ci
    return val
#end def

def Rs(x,dx,s,c):
    if x+1-s<-dx:
        return 0-(1-s)
    elif x+1-s>dx:
        return x
    else:
        return p7(x+1-s,c)-(1-s)
#end def

def unboundedness(inputfile,db=0.25,dbs=0.005):

    f1=open(inputfile, "r")
    
    d1=np.array(f1.read().split(),dtype=float)
    f1.close()
    
    # HEADER
    d1head=d1[0:2]
    # CHANNEL LENGTHS
    nchan=int(d1head[1])
    d1chan=d1[2:2+nchan]
    # POTENTIALS
    d1pots=[]
    pos=0
    for i in range(len(d1chan)):
        dtmp = d1[2+nchan+int(pos):2+nchan+int(pos+d1chan[i]*3)]
        dtmp.shape=len(dtmp)//3,3
        pos+=d1chan[i]*3
        d1pots.append(dtmp)
    #end for
    
    A=[]
    for i in range(8):
        row=[]
        for j in range(8):
            if i<4:
                if j-i<0:
                    dcoeff=0
                    dpower=0
                else:
                    dcoeff=math.factorial(j)/math.factorial(j-i)
                    dpower=j-i
                #end if
                row.append(dcoeff*db**dpower)
            else:
                if j-(i-4)<0:
                    dcoeff=0
                    dpower=0
                else:
                    dcoeff=math.factorial(j)/math.factorial(j-(i-4))
                    dpower=j-(i-4)
                #end if
                row.append(dcoeff*(-db)**dpower)
            #end if
        #end for
        A.append(row)
    #end for
    
    A = np.array(A)
    b = np.array([db,1]+[0]*6)
    c = np.linalg.inv(A).dot(b)
    
    ng=3000
    gmin=0.02
    gmax=0.55
    r = np.linspace( gmin, gmax, ng )
    
    # PP
    v = []
    for i in range(nchan):
        vtmp = []
        for j in range(len(r)):
            vtmp.append(ppchannel(r[j],d1pots[i][:,2],d1pots[i][:,1],d1pots[i][:,0]))
        #for
        v.append(vtmp)
    #for
    v=np.array(v)
    
    # 2*r^2*VL2 
    f = r*r*(v[1]-v[0])
    # 2*r^2*V'L2 
    fp = [Rs(fr,db,dbs,c) for fr in f]
    
    undoundedness = 0
    for fi,fx in enumerate(f):
        undoundedness+=(fp[fi]-fx)**2.*(gmax-gmin)/ng 
    #end for
    return undoundedness

#end def

def get_yoon_norm(fort26,rcuts):

    rcuts = np.array(rcuts,dtype=float)
    
    orbs = []
    orb = []
    read = True
    garbage_count = 0
    with open(fort26) as f:
        for num, line in enumerate(f, 1):
            if '  20. 0.' in line or not read:
                
                garbage_count+=1
                if garbage_count==1:
                    read=False
                    orbs.append(orb)
                #end if
    
                if garbage_count==3:
                    read=True
                    orb = []
                    garbage_count=0
                #end if
    
            else:
        
                x,y = line.split()
                # save orb data
                orb.append([x,y])
                #norm+=
    
    orbs = np.array(orbs,dtype=float)
    
    norms=[]
    for oi in range(len(orbs)):
        ovlp=0
        
        x10=orbs[oi][0][0]
        y10=orbs[oi][0][1]
    
        for vi in range(1,len(orbs[oi])): 
    
            if orbs[oi][vi][0]<rcuts[oi]:
    
                dx=orbs[oi][vi][0]-x10
    
                x1 = orbs[oi][vi][0]
                y1 = orbs[oi][vi][1]
    
                ovlp+=dx*((x10*y10)*(x10*y10)+(x1*y1)*(x1*y1))/2.0
    
                x10=x1
                y10=y1
            #end if
        #end for

        norms.append(ovlp)
    #end for
    return np.array(norms)

#end def

def get_yoon_norm_diff(fort26,ref_fort26,rcuts):

    rcuts = np.array(rcuts,dtype=float)
    
    orbs1 = []
    orb = []
    read = True
    garbage_count = 0
    with open(fort26) as f:
        for num, line in enumerate(f, 1):
            if '  20. 0.' in line or not read:
                
                garbage_count+=1
                if garbage_count==1:
                    read=False
                    orbs1.append(orb)
                #end if
    
                if garbage_count==3:
                    read=True
                    orb = []
                    garbage_count=0
                #end if
    
            else:
        
                x,y = line.split()
                # save orb data
                orb.append([x,y])
                #norm+=
    
    orbs1 = np.array(orbs1,dtype=float)
    
    
    orbs2 = []
    orb = []
    read = True
    garbage_count = 0
    with open(ref_fort26) as f:
        for num, line in enumerate(f, 1):
            if '  20. 0.' in line or not read:
                
                garbage_count+=1
                if garbage_count==1:
                    read=False
                    orbs2.append(orb)
                #end if
    
                if garbage_count==3:
                    read=True
                    orb = []
                    garbage_count=0
                #end if
    
            else:
        
                x,y = line.split()
                # save orb data
                orb.append([x,y])
                #norm+=
    
    orbs2 = np.array(orbs2,dtype=float)
   
    diffs=[]
    for oi in range(len(orbs1)):
        ovlp1=0
        ovlp2=0
        
        x10=orbs1[oi][0][0]
        y10=orbs1[oi][0][1]
        x20=orbs2[oi][0][0]
        y20=orbs2[oi][0][1]
    
        for vi in range(1,len(orbs1[oi])): 
    
            if orbs1[oi][vi][0]<rcuts[oi]:
    
                dx=orbs1[oi][vi][0]-x10
    
                x1 = orbs1[oi][vi][0]
                y1 = orbs1[oi][vi][1]
    
                x2 = orbs2[oi][vi][0]
                y2 = orbs2[oi][vi][1]
    
                ovlp1+=dx*((x10*y10)*(x10*y10)+(x1*y1)*(x1*y1))/2.0
                ovlp2+=dx*((x20*y20)*(x20*y20)+(x2*y2)*(x2*y2))/2.0
    
                x10=x1
                y10=y1
                x20=x2
                y20=y2
            #end if
        #end for

        diffs.append(ovlp1-ovlp2)
    #end for
    return np.array(diffs)

#end def

# CORRECTING UNBOUNDEDNESS

def p7(x,c):
    val=0
    for ci,cv in enumerate(c):
        val+=cv*x**ci
    return val
#end def

def Rs(x,dx,s,c):
    if x+1-s<-dx:
        return 0-(1-s)
    elif x+1-s>dx:
        return x
    else:
        return p7(x+1-s,c)-(1-s)
#end def

class fitClass:

    def __init__(self):
        pass

    def gauss_correction(self,x,c1,c2,c3):
        val = 0
        for ci,c in enumerate([c1,c2,c3]):
            val+=x**2.*c*np.exp(-self.exps[ci]*x**2.)
        #end for
        return val
    #end def

#end class

def make_bound(exps0,pots,db,dbs):

    nchan = len(pots)

    A=[]
    for i in range(8):
        row=[]
        for j in range(8):
            if i<4:
                if j-i<0:
                    dcoeff=0
                    dpower=0
                else:
                    dcoeff=math.factorial(j)/math.factorial(j-i)
                    dpower=j-i
                #end if
                row.append(dcoeff*db**dpower)
            else:
                if j-(i-4)<0:
                    dcoeff=0
                    dpower=0
                else:
                    dcoeff=math.factorial(j)/math.factorial(j-(i-4))
                    dpower=j-(i-4)
                #end if
                row.append(dcoeff*(-db)**dpower)
            #end if
        #end for
        A.append(row)
    #end for
    
    A = np.array(A)
    b = np.array([db,1]+[0]*6)
    c = np.linalg.inv(A).dot(b)

    ng=3000
    gmin=0.02
    gmax=0.85
    r = np.linspace( gmin, gmax, ng )
    
    # PP
    v = []
    for i in range(nchan):
        vtmp = []
        for j in range(len(r)):
            vtmp.append(ppchannel(r[j],pots[i][:,2],pots[i][:,1],pots[i][:,0]))
        #for
        v.append(vtmp)
    #for
    v=np.array(v)
    
    # 2*r^2*VL2 
    f = r*r*(v[1]-v[0])
    # 2*r^2*V'L2 
    fp = [Rs(fr,db,dbs,c) for fr in f]

    undoundedness = 0
    for fi,fx in enumerate(f):
        undoundedness+=(fp[fi]-fx)*(gmax-gmin)/ng 
    #end for
    print('\npseudopotential undoundedness: ',undoundedness)

    import matplotlib.pyplot as plt
    from scipy.optimize import curve_fit

    plt.plot(r, fp, 'g-', label='fp')


    fit_instance = fitClass()
    fit_instance.exps=exps0
    popt, pcov = curve_fit(fit_instance.gauss_correction, r, f-fp)

    plt.xlabel('r (bohr)')
    plt.ylabel('$2r^2v_{L^2}$')
    plt.plot(r, f-fit_instance.gauss_correction(r, *popt), 'r-',label='f-corr')
    plt.plot(r, f, 'b-',label='f')
    plt.plot(r, [-1]*len(r), 'k-',label=None)
    plt.legend()
    plt.show()

    return popt

#end def
