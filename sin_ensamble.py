#!/usr/bin/python
import sys
import numpy as np
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import cosmolopy.distance as cd
import time
#import my_auto_ssp_elines_rnd as ssp
import pyfits as pyf
import cosmolopy.distance as cd
from pyfits import getheader as ghead
from pyfits import getdata as gdata
from pyfits import writeto as wfit
from scipy.interpolate.interpolate import interp1d
#import my as my
import math
import os.path as ptt
import matplotlib
matplotlib.use('Agg')

def wfits_ext(name,hlist):
    if ptt.exists(name) == False:
        hlist.writeto(name)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        hlist.writeto(name)

def func_plot(x,ftype=0):
    if ftype == 0:
        y=10**(x)
    if ftype == 1:
        y=np.log10(x)
    if ftype == 2:
        y=x
    return y

def sycall(comand):
    import os
    os.system(comand)

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        sycall("rm "+name)
        wfit(name,data,hdr)

def map_plot(image,ages_l,rads,pdf=0,form='pdf',dir='',fname='map',title='map',minval=6,maxval=8):
    #image[np.where(np.isfinite(image) == False)]=0
    mapf=np.log10(image+0.0001)
    #print np.sum(mapf),"THIS IS THE SUM OF MASS",np.amin(mapf),np.amax(mapf)
    [nx,ny]=image.shape
    #nx=nx*.5
    #ny=ny*.5
    ft=open("ctable", "r")
    nc=256
    g=np.zeros(nc)
    r=np.zeros(nc)
    b=np.zeros(nc)
    l=np.zeros(nc)
    for cmap in ft:
        data=cmap.split(" ")
        data=filter(None,data)
        nc=int(float_(data[0]))
        r[nc-1]=(float_(data[1]))/255.
        g[nc-1]=(float_(data[2]))/255.
        b[nc-1]=(float_(data[3]))/255.
        l[nc-1]=nc/255.
    bright=1 
    contrast=0.5
    r[0]=1.0
    g[0]=1.0
    b[0]=1.0
    nc=256
    my_rgb=np.zeros([256,3])
    my_rgb[:,0]=r
    my_rgb[:,1]=g
    my_rgb[:,2]=b
    x1=np.arange(0,nx,1)-nx/2+1
    y1=np.arange(0,ny,1)-ny/2+1
    [X,Y]=np.meshgrid(x1*.5,y1*.5)
    #my_rgb[:,3]=l
    my_cmap = matplotlib.colors.ListedColormap(my_rgb, name='my_name')
    #if pdf == 1:
    #    matplotlib.use('Agg')
    from matplotlib import cm
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    font0 = FontProperties()
#    fig, ax = plt.subplots(figsize=(6,5.5))
    fig= plt.figure(figsize=(6,5.5))
#    ax=fig.add_subplot(1, 1, 1)
    ax = fig.add_axes([0.08, 0.11, 0.78, 0.82])
    plt.xlabel("RA (arcsec)",fontsize=14)
    plt.ylabel("DEC (arcsec)",fontsize=14)
    ax.yaxis.set_label_coords(-0.05, 0.5)
    plt.title(title,fontsize=15)
    font = font0.copy()
    font.set_weight('bold')
    plt.text(-nx*.5/2.*0.9,ny*.5/2*0.88,'$%3.2f' % 10.0**(ages_l-9.0) + '$ $Gyr$', fontsize=15)
    font.set_weight('normal')
    levels=[0.5,1,1.5]
    labels = ['$%3.1f' % val + " R_{50}%$" for val in levels ]
    #print labels,val
    fmt={}
    for l,s in zip( levels, labels ):
        fmt[l] = s
#    im = ax.contourf(image, 150, cmap=cm.coolwarm)#, extent=(mi,mf,np.amax(age),np.amin(age)))
    im=ax.imshow(mapf,cmap=my_cmap,vmin=minval, vmax=maxval, interpolation='nearest', origin='lower', extent=[-nx*.5/2,nx*.5/2,-ny*.5/2,ny*.5/2])
    CS = ax.contour(X,Y,rads, levels, linewidths=2, extent=(-nx*.5/2,nx*.5/2,-ny*.5/2,ny*.5/2), colors=('black','black','black'))#cmap=cm.gray
    ax.clabel(CS, levels, inline=1,fmt=fmt,fontsize=14)
    cbar_ax = fig.add_axes([0.85, 0.11, 0.03, 0.82])
    cbf=fig.colorbar(im, cax=cbar_ax)
    cbf.set_label('$log[M_{\odot}$ /arcsec$^2]$',fontsize=15)#+plus_fon_size)
#    ax.annotate(ages_l,xytext=(12, -15))
    if pdf == 1:
#        fig.tight_layout()
        plt.savefig(dir+fname+'.'+form)
    else:
        plt.show()
    plt.close()
#sys.argv=filter(None,sys.argv)
#if len(sys.argv) < 4:
#    print "USE: pack_results_name.py NAME PREFIX"
    #sys.exit(0)
#name=sys.argv[1]
#dir1=sys.argv[2]
#dir2=sys.argv[2]

def int_ensamble(name,dir1,dir3,dir4,fo,fi=0,pdf=1,fits_f=0,m_t=""):
    names=name.split("-")
    dir3=dir3#+"/"+names[1]+"/"+names[2]
    dir_map=dir3+"/"+names[1]+"/"+names[2]+m_t
    DIRS=dir3.split("/")
    DRT=""
    for DR in DIRS:
        DRT=DRT+DR+"/"
        call="mkdir -p "+DRT
        sycall(call)
    DIRS=dir_map.split("/")
    DRT=""
    for DR in DIRS:
        DRT=DRT+DR+"/"
        call="mkdir -p "+DRT
        #sycall(call)
    call="mkdir -p "+dir3+'/Plots'
    sycall(call)
    speed_of_light=299792.458
    #dir1=dir1+"/"+names[1]+"/stadi/"+names[2]+m_t
    dir1=dir1+"/"+names[1]+"/"+names[2]+m_t
    dir2=dir1
    file1=dir1+"/coeffs_auto_ssp."+name+".int.out"
    file2=dir2+"/auto_ssp."+name+".int.out"
    pdl_cube=[]
    pdl_cube_e=[]
    ages=[]
    cont=0
    f1=open(file1,"r")
    for line in f1:
        if not "#" in line:
            line=line.replace("\n","")
            data=line.split(" ")
            data=filter(None,data)
            if (cont % 4) == 0:
                ages.extend([np.log10(float_(data[1]))+9.0])
            pdl_cube.extend([float_(data[3])])
            pdl_cube_e.extend([float_(data[8])])
            cont=cont+1
    f1.close()
    f2=open(file2,"r")
    for line in f2:
        if not "#" in line:
            line=line.replace("\n","")
            data=line.split(",")
            data=filter(None,data)
            pdl_flux=float_(data[13])
            pdl_flux_e=float_(data[14])
            pdl_Av=float_(data[5])
            pdl_Av_e=float_(data[6])
            redshift=float_(data[7])
    f2.close()
    ages=np.array(ages)
    ages=sorted(ages,reverse=True)
    #np.sor
    #sys.exit()
    pdl_cube=np.array(pdl_cube)
    pdl_cube_e=np.array(pdl_cube_e)


#    f=open(dir4+"/BASE.gsd01","r")
    f=open(dir4+"/BASE.bc17_salp_Agelin_Metlin_330","r")
    yunk=f.readline()
    age_t=[]
    met_t=[]
    cor_t=[]
    for line in f:
        if not "#" in line:
            data=line.split(" ")
            data=filter(None,data)
            age_t.extend([float_(data[1])])
            met_t.extend([float_(data[2])])
            cor_t.extend([float_(data[4])])
    n_t=len(age_t)
    age_t=np.array(age_t)
    met_t=np.array(met_t)
    cor_t=np.array(cor_t)
    age_t=np.around(age_t/1e9,decimals=4)
    met_t=np.around(met_t,decimals=4)
    f.close()

    cosmo = {'omega_M_0' : 0.27, 'omega_lambda_0' : 0.73, 'h' : 0.71}
    cosmo = cd.set_omega_k_0(cosmo)
    DL1=cd.luminosity_distance(redshift,**cosmo)
    #print DL, redshift
    ratio=3.08567758e24
    modz=5.0*np.log10(DL1)+25.0
    DL=DL1*ratio
    DA=DL1/(1+redshift)**2.0*1e6*np.pi/180./3600.
    L=4.0*np.pi*(DL**2.0)#/(1+$redshift);
    Factor=(L*1e-16)/3.826e33
    filed=file1
    f2=open(filed,"r")
    n=0
    n_ini=0
    ML=np.zeros(156)
    a_age=[]
    a_met=[]
    n_age=0
    n_met=0
    AGE=np.zeros(156)
    MET=np.zeros(156)
    COR=np.zeros(156)
    for line in f2:
        if n_ini < 156:
            if not "#" in line:
                data=line.split(" ")
                data=filter(None,data)
                n=int(data[0])
                AGE[n]=float_(data[1])
                MET[n]=float_(data[2])
                ML[n]=float_(data[5])
                diff_age=1
                for i in range(0, n):
                    if AGE[n] == AGE[i]:
                        diff_age=0
                if diff_age == 1:
                    a_age.extend([AGE[n]])
                    n_age=n_age+1
                diff_met=1
                for i in range(0, n):
                    if MET[n] == MET[i]:
                        diff_met=0
                if diff_met == 1:
                    a_met.extend([MET[n]])
                    n_met=n_met+1
                n_ini=n_ini+1
                for jt in range(0, n_t):
                    if age_t[jt] == 0.02:
                        age_t[jt] = 0.0199
                    if AGE[n] == age_t[jt]:
                        if MET[n] == met_t[jt]:
                            COR[n]=cor_t[jt]
    f2.close()
    n=n+1
    MassT=0
    LighT=0
    massN=0
    massN_e=0
    lightN=0
    lightN_e=0
    #mas1=0
    #mas2=0
    #mas3=0
    #mas4=0
    mass=np.zeros([n_age])
    mass_e=np.zeros([n_age])
    light=np.zeros([n_age])
    light_e=np.zeros([n_age])
    #ages=np.zeros([n_age])
    sfrt=np.zeros([n_age])
    sfdt=np.zeros([n_age])
    nz=len(pdl_cube)
    temp_a=0
    mass_age=0
#    pdl_cube[np.isnan(pdl_cube)]=1

    #sys.exit()
    #for i in range(nt-1, 155, -1):
    #    label=hdr['FILE_'+str(i)]
    #    time=label.replace('_NORM_age.fits.gz','')
    #    time=float_(time.replace('map.CS.'+name+'_',''))
    #    ages[38-i+156]=np.log10(time)+9
    #print np.log10(time)+9, 38-i+156
    #temp=pdl_cube[i,:,:]
    #temp=temp*10.0**(ML[i-156])*pdl_flux*Factor
    #temp_a=temp+temp_a
#    MassT=MassT+np.sum(temp)
#    mas1=np.sum(temp[n1])+mas1
#    mas2=np.sum(temp[n2][n2a])+mas2
#    mas3=np.sum(temp[n3][n3a])+mas3
#    mas4=np.sum(temp[n4][n4a])+mas4
#    mass[38-i+156,0]=np.log10(mas1)
#    mass[38-i+156,1]=np.log10(mas2)
#    mass[38-i+156,2]=np.log10(mas3)
#    mass[38-i+156,3]=np.log10(mas4)
    f2=open(dir3+"/"+name+m_t+"_Ensemble_int.csv","w")
    f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    fL2=open(dir3+"/"+name+m_t+"_Ensemble_int_L.csv","w")
    fL2.write("#  LOG_AGE  N_LIGHTR_1  N_LIGHTR_2  N_LIGHTR_3  N_LIGHTR_4  LOG_LIGHTR_1  LOG_LIGHTR_2  LOG_LIGHTR_3  LOG_LIGHTR_4 \n")
    mass_age_t=np.zeros([n_age])
    mass_age_t_2=np.zeros([n_age])
    mass_age_t_e=np.zeros([n_age])
    light_age_t=np.zeros([n_age])
    light_age_t_2=np.zeros([n_age])
    light_age_t_e=np.zeros([n_age])
    for i in range(0, n_age):
        age_now=a_age[i]
        pdl_age=0
        pdl_age_2=0
        pdl_age_e=0
        pdl_ageL=0
        pdl_age_2L=0
        pdl_age_eL=0
        for j in range(0, n):
            if age_now == AGE[j]:
                #if AGE[j] <= 2:
                pdl_age=pdl_age+pdl_cube[j]*10.0**(ML[j])*pdl_flux*Factor*10.0**(0.4*pdl_Av)#*0.25/np.pi#/1.47
                pdl_age_e=pdl_age_e+((pdl_cube_e[j]/pdl_cube[j])**2.0+(pdl_flux_e/pdl_flux)**2.0+(np.log(10.0)*0.4*pdl_Av_e)**2.0)*(pdl_cube[j]*10.0**(ML[j])*pdl_flux*Factor*10.0**(0.4*pdl_Av))**2.0
                pdl_age_2=pdl_age_2+pdl_cube[j]*10.0**(ML[j])*pdl_flux*Factor*10.0**(0.4*pdl_Av)/COR[j]
                pdl_age_2L=pdl_age_2L+pdl_cube[j]*pdl_flux*Factor*10.0**(0.4*pdl_Av)/COR[j]
                pdl_ageL=pdl_ageL+pdl_cube[j]*pdl_flux*Factor*10.0**(0.4*pdl_Av)
                pdl_age_eL=pdl_age_eL+((pdl_cube_e[j]/pdl_cube[j])**2.0+(pdl_flux_e/pdl_flux)**2.0+(np.log(10.0)*0.4*pdl_Av_e)**2.0)*(pdl_cube[j]*pdl_flux*Factor*10.0**(0.4*pdl_Av))**2.0
                #print age_now    
        
        if np.isfinite(pdl_age) == False:
            pdl_age=0
        if np.isfinite(pdl_age_2) == False:
            pdl_age_2=0
        if np.isfinite(pdl_age_e) == False:
            pdl_age_e=0
        if np.isfinite(pdl_ageL) == False:
            pdl_ageL=0
        if np.isfinite(pdl_age_2L) == False:
            pdl_age_2L=0
        if np.isfinite(pdl_age_eL) == False:
            pdl_age_eL=0
          #  pdl_age_e[np.isnan(pdl_age_e)]=0
        for k in range(0, n_age):
            if np.log10(age_now)+9 == ages[k]:
                #print ages[k],np.log10(age_now)+9,age_now
                mass_age_t[k]=pdl_age
                mass_age_t_2[k]=pdl_age_2
                mass_age_t_e[k]=pdl_age_e
                light_age_t[k]=pdl_ageL
                light_age_t_2[k]=pdl_age_2L
                light_age_t_e[k]=pdl_age_eL
    #print mass_age_t

    mass_temp_total=0
    light_temp_total=0
    #sys.exit()
    for i in range(0, n_age):
        if i == 0:
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
            age_i=0.0
        elif i == n_age-1:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=2.0*10.0**(ages[i])-age_i
        else:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
        Dt_age=np.abs(age_s-age_i)
        sfdt[i]=Dt_age/1e6
        temp=mass_age_t[i]
        temp_2=mass_age_t_2[i]
        temp_e=mass_age_t_e[i]
        tempL=light_age_t[i]
        temp_2L=light_age_t_2[i]
        temp_eL=light_age_t_e[i]
        #temp[np.where(np.isfinite(temp) == False)]=0
        #temp_e[np.where(np.isfinite(temp_e) == False)]=1
        #temp_2[np.where(np.isfinite(temp_2) == False)]=0
        if i == 0:
            if fits_f == 1:
                MASS_map_cube=np.zeros([n_age])
                MGH_map_cube=np.zeros([n_age])
                SFH_map_cube=np.zeros([n_age])
                LIGHT_map_cube=np.zeros([n_age])
                LGH_map_cube=np.zeros([n_age])
            temp1=temp
            temp1L=tempL
        else:
            temp1=temp1+temp
            temp1L=temp1L+tempL
        if fits_f == 1:
            MASS_map_cube[i]=temp
            MGH_map_cube[i]=temp1
            SFH_map_cube[i]=temp_2/Dt_age
            LIGHT_map_cube[i]=tempL
            LGH_map_cube[i]=temp1L
        #if pdf==1:
        #map_plot(temp1,ages[i],pdl_rad,dir=dir_map+"/",pdf=1,title=name,form='pdf',fname=name+'_smap_'+str(i),minval=lovalue,maxval=upvalue)
        MassT=MassT+np.sum(temp)
        LighT=LighT+np.sum(tempL)
#            print temp[ind[ii]]
#            print len(temp[ind[ii]]),np.amax(inda[ii])
        Dt_mass=temp
        Dt_mass_e=temp_e
        Dt_mass_2=temp_2
        Dt_light=tempL
        Dt_light_e=temp_eL
        Dt_light_2=temp_2L
        massN=Dt_mass+massN
        massN_e=Dt_mass_e+massN_e
        lightN=Dt_light+lightN
        lightN_e=Dt_light_e+lightN_e
        #mas2=np.sum(temp[n2][n2a])+mas2
        #mas3=np.sum(temp[n3][n3a])+mas3
        #mas4=np.sum(temp[n4][n4a])+mas4
#            print temp[ind[ii]][inda[ii]]
        mass[i]=np.log10(massN)
        mass_e[i]=massN_e
        light[i]=np.log10(lightN)
        light_e[i]=lightN_e
        sfrt[i]=Dt_mass_2/Dt_age#/massN[ii]#/(np.pi*(rx[ii+1]**2.0-rx[ii]**2.0)*rad**2.0)#/(float_(len(temp[ind[ii]][inda[ii]]))*(0.5*DA)**2.0)
        #mass[i,1]=np.log10(mas2)
        #mass[i,2]=np.log10(mas3)
        #mass[i,3]=np.log10(mas4)
       # print mass[i],massN
    mass_temp_total=np.log10(np.sum(10**mass[n_age-1]))
    light_temp_total=np.log10(np.sum(10**light[n_age-1]))
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    LighT=np.log10(LighT)
    LighT=light_temp_total
    if fits_f == 1:
        #if ptt.exists() == True:
        h1=pyf.PrimaryHDU(MASS_map_cube)#.header
        h2=pyf.PrimaryHDU(SFH_map_cube)#.header
        h3=pyf.PrimaryHDU(MGH_map_cube)
        h4=pyf.PrimaryHDU(LIGHT_map_cube)
        h5=pyf.PrimaryHDU(LGH_map_cube)
        h=h1.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_age
        h["NAXIS2"]=1
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"MASS_maps_"+name+".fits",hlist)
        #if ptt.exists() == True:
        h=h2.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_age
        h["NAXIS2"]=1
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"SFH_maps_"+name+".fits",hlist) 
        #if ptt.exists() == True:   
        h=h3.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_age
        h["NAXIS2"]=1
        hlist=pyf.HDUList([h3])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"MGH_maps_"+name+".fits",hlist)
        h=h4.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_age
        h["NAXIS2"]=1
        hlist=pyf.HDUList([h4])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"LIGHT_maps_"+name+".fits",hlist)
        h=h5.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_age
        h["NAXIS2"]=1
        hlist=pyf.HDUList([h5])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"LGH_maps_"+name+".fits",hlist)
    #print MassT,name

    #print Ha,(L*1e-16)
    #sys.exit(0)
    mass_n=10**(10**(mass-mass[n_age-1]))
    mass_n_e=np.sqrt((10**(mass-mass[n_age-1]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[n_age-1]/10**(2.0*mass[n_age-1]))))
    light_n=10**(10**(light-light[n_age-1]))
    light_n_e=np.sqrt((10**(light-light[n_age-1]))**2.0*((light_e/10**(2.0*light))+(light_e[n_age-1]/10**(2.0*light[n_age-1]))))
    #print mass_n
    #mass_n=10**(mass-mass[nt-156-1,:])
    #mass_n=(mass-mass[nt-156-1,:])
    for i in range(0, n_age):
        #print ages[i],a_age[i],"test_ages"
        line=''
        line=line+str(ages[i])
        line=line+';'+str(mass_n[i])
        line=line+';'+str(mass[i])
        line=line+';'+str(sfrt[i])
        line=line+';'+str(mass_n_e[i])
        line=line+';'+str(sfdt[i])
        #print line
        line=line+' \n'
        f2.write(line)
        lineL=''
        lineL=lineL+str(ages[i])
        lineL=lineL+';'+str(light_n[i])
        lineL=lineL+';'+str(light[i])
        lineL=lineL+';'+str(sfrt[i])
        lineL=lineL+';'+str(light_n_e[i])
        lineL=lineL+';'+str(sfdt[i])
        lineL=lineL+' \n'
        fL2.write(lineL)
    
   # print mass_n.shape,ages.shape
    if not pdf == 0:
        dev=dir3+'/Plots/'+name+m_t+"_int_Relative_Mass2.pdf"
        #if pdf == 1:
        #    matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6,5.5))
        ax.set_xlabel("$log_{10}(time/yr)$",fontsize=14)
        ax.set_ylabel("$M(t)/M_{0}$",fontsize=14)
        #MassT=10.32
        ax.set_title(name+' $\log M_{tot}='+('%7.2f' % MassT)+'$',fontsize=15)
        ax.set_xlim(8.6,10.1)
            #ax.set_ylim(0,12)
        ax.set_ylim(func_plot(np.log10(1.78),ftype=1),func_plot(np.log10(12),ftype=1))
        plt.plot(ages,func_plot(np.log10(mass_n),ftype=1))#,label='$'+('%6.1f' % rx[ii])+'R_e<R<'+('%6.1f' % rx[ii+1])+'R_e$'
        plt.legend(loc=3)
        plt.plot(np.arange(0,20,.1),np.ones(200)*func_plot(0.95,ftype=1),'--',color='black')
        plt.plot(np.arange(0,20,.1),np.ones(200)*func_plot(0.50,ftype=1),'--',color='green')
        fig.canvas.draw()
        labels = [item.get_text() for item in ax.get_yticklabels()] 
        for i in range(0, len(labels)):
            labels[i]=labels[i].replace(u'\u2212','-')
        for i in range(0, len(labels)):
            if labels[i] != u'':
                if float_(labels[i]) == 0:
                    labels[i]=u'%3.2f' % 10**(0)
                else:
                    labels[i]=u'%3.2f' % 10**(float_(labels[i]))
        ax.set_yticklabels(labels)
        if pdf == 1:
            fig.tight_layout()
            plt.savefig(dev)#,dpi = 1000)
        else:
            plt.show()
        plt.close()
    f2.close()
    fL2.close()
    fo.write(name+" "+str(MassT)+" "+str(redshift)+" \n")
    

def ensamble(name,dir1,dir3,dir4,fo,fi=0,fii=0,pdf=1,rx=[0,0.5,1.0,1.5],fits_f=0):
    rad=1.0
    names=name.split("-")
    dir3=dir3#+"/"+names[1]+"/"+names[2]
    dir_map=dir3+"/"+names[1]+"/"+names[2]
    DIRS=dir3.split("/")
    DRT=""
    for DR in DIRS:
        DRT=DRT+DR+"/"
        call="mkdir -p "+DRT
        sycall(call)
    DIRS=dir_map.split("/")
    DRT=""
    for DR in DIRS:
        DRT=DRT+DR+"/"
        call="mkdir -p "+DRT
        sycall(call)
    speed_of_light=299792.458
    dir1=dir1+"/"+names[1]+"/"+names[2]
#    dir2=dir2+"/"+names[1]+"-"+names[2]
    file=dir1+"/"+name+".SFH.cube.fits.gz"
    file2=dir1+"/"+name+".p_e.pdl_r.fits"
    #file2=dir2+"/"+name+".photo.r_Lc_rad.fits"
    #file3=dir1+"/"+'mask.'+name+'.V.fits.gz'
    file3=dir1+"/"+'DMASK.'+name+'.fits.gz'
    [pdl_cube, hdr]=gdata(file,0, header=True)
    [pdl_rad, hdr2]=gdata(file2,0, header=True)
    [pdl_mask, hdr3]=gdata(file3,0, header=True)
    pdl_mask=1.0-pdl_mask
    if np.sum(pdl_mask) == 0:
        pdl_mask[:,:]=1.0
    ind=[]
    inda=[]
    nr=len(rx)
    for ii in range(0, nr-1):
        nt=np.where(pdl_rad< rx[ii+1]*rad)
        nta=np.where(pdl_rad[nt]>= rx[ii]*rad)
        ind.extend([nt])
        inda.extend([nta])
#    n2=np.where(pdl_rad< r2*rad)
#    n2a=np.where(pdl_rad[n2]>= r1*rad)
#    n3=np.where(pdl_rad< r3*rad)
#    n3a=np.where(pdl_rad[n3]>= r2*rad)
#    n4=np.where(pdl_rad< r4*rad)
#    n4a=np.where(pdl_rad[n4]>= r3*rad)
    SN_file="norm_SN_"+name+".CS.fits.gz"
    if ptt.exists(dir1+"/"+SN_file) == True:
        [pdl_SN, hdr000]=gdata(dir1+"/"+SN_file, 0, header=True)   
    Ha_file="map.CS."+name+"_flux_6562.fits.gz"
    if ptt.exists(dir1+"/"+Ha_file) == True:
        [pdl_ha, hdr001]=gdata(dir1+"/"+Ha_file,0, header=True)
        pdl_ha=pdl_ha*pdl_mask#[0,:,:]#-z_r*speed_of_light
        pdl_ha[np.isnan(pdl_ha)]=0     
    Av_file="map.CS."+name+"_Av_ssp.fits.gz"
    [pdl_Av, hdr002]=gdata(dir1+"/"+Av_file,0, header=True)
    Av_file_e="map.CS."+name+"_e_Av_ssp.fits.gz"
    if ptt.exists(dir1+"/"+Av_file_e) == True:
        [pdl_Av_e, hdr002e]=gdata(dir1+"/"+Av_file_e,0, header=True)
        pdl_Av_e[np.isnan(pdl_Av_e)]=0
    else:
        pdl_Av_e=np.zeros(Av_file.shape)
    nt=hdr['NAXIS3']-5#5#4#n_met
    flux_file="map.CS."+name+"_flux_ssp.fits.gz"
    [pdl_flux, hdr0]=gdata(dir1+"/"+flux_file, 0, header=True)
    flux_file_e="map.CS."+name+"_e_flux_ssp.fits.gz"
    if ptt.exists(dir1+"/"+flux_file_e) == True:
        [pdl_flux_e, hdr0e]=gdata(dir1+"/"+flux_file_e, 0, header=True)
        pdl_flux_e[np.isnan(pdl_flux_e)]=0
    else:
        pdl_flux_e=np.zeros(pdl_flux.shape)
    mass_file="map.CS."+name+"_Mass_dust_cor_ssp.fits.gz"#dust_cor_
    [pdl_mass, hdr00]=gdata(dir1+"/"+mass_file, 0, header=True)
    pdl_mass[np.isnan(pdl_mass)]=1
    MassT2=np.log10(np.sum(10.0**pdl_mass))
    #print MassT2, name
    
    f=open(dir4+"/BASE.gsd01","r")
    #f=open(dir4+"/BASE.bc17_salp_Agelin_Metlin_330","r")
    yunk=f.readline()
    age_t=[]
    met_t=[]
    cor_t=[]
    for line in f:
        if not "#" in line:
            data=line.split(" ")
            data=filter(None,data)
            age_t.extend([float_(data[1])])
            met_t.extend([float_(data[2])])
            cor_t.extend([float_(data[4])])
    n_t=len(age_t)
    age_t=np.array(age_t)
    met_t=np.array(met_t)
    cor_t=np.array(cor_t)
    age_t=np.around(age_t/1e9,decimals=4)
    met_t=np.around(met_t,decimals=4)
    f.close()
    a_redshift=[]
    filet="auto_ssp.CS."+name+".rss.out"
    f=open(dir1+"/"+filet, "r")
    for line in f:
        if not "#" in line:
            data=line.split(",")
            data=filter(None,data)
            #print data
            a_redshift.extend([float_(data[7])])
    f.close()
    a_redshift=np.array(a_redshift)
    redshift=np.median(a_redshift)
    cosmo = {'omega_M_0' : 0.27, 'omega_lambda_0' : 0.73, 'h' : 0.71}
    cosmo = cd.set_omega_k_0(cosmo)
    DL1=cd.luminosity_distance(redshift,**cosmo)
    #print DL, redshift
    ratio=3.08567758e24
    modz=5.0*np.log10(DL1)+25.0
    DL=DL1*ratio
    DA=DL1/(1+redshift)**2.0*1e6*np.pi/180./3600.
    L=4.0*np.pi*(DL**2.0)#/(1+$redshift);
    Factor=(L*1e-16)/3.826e33
    filed="coeffs_auto_ssp.CS."+name+".rss.out"
    f2=open(dir1+"/"+filed,"r")
    n=0
    n_ini=0
    n_ssp=156#330
    ML=np.zeros(n_ssp)
    a_age=[]
    a_met=[]
    n_age=0
    n_met=0
    AGE=np.zeros(n_ssp)
    MET=np.zeros(n_ssp)
    COR=np.zeros(n_ssp)
    for line in f2:
        if n_ini < n_ssp:
            if not "#" in line:
                data=line.split(" ")
                data=filter(None,data)
                n=int(data[0])
                AGE[n]=float_(data[1])
                MET[n]=float_(data[2])
                ML[n]=float_(data[5])
                diff_age=1
                for i in range(0, n):
                    if AGE[n] == AGE[i]:
                        diff_age=0
                if diff_age == 1:
                    a_age.extend([AGE[n]])
                    n_age=n_age+1
                diff_met=1
                for i in range(0, n):
                    if MET[n] == MET[i]:
                        diff_met=0
                if diff_met == 1:
                    a_met.extend([MET[n]])
                    n_met=n_met+1
                n_ini=n_ini+1
                for jt in range(0, n_t):
                    if age_t[jt] == 0.02:
                        age_t[jt] = 0.0199
                    if AGE[n] == age_t[jt]:
                        if MET[n] == met_t[jt]:
                            COR[n]=cor_t[jt]
    f2.close()
    n=n+1
    MassT=0
    LighT=0
    massN=np.zeros(nr)
    massN_e=np.zeros(nr)
    lightN=np.zeros(nr)
    lightN_e=np.zeros(nr)
    #mas1=0
    #mas2=0
    #mas3=0
    #mas4=0
    mass=np.zeros([n_age,nr])
    mass_e=np.zeros([n_age,nr])
    light=np.zeros([n_age,nr])
    light_e=np.zeros([n_age,nr])
    ages=np.zeros([n_age])
    sfrt=np.zeros([n_age,nr])
    sfdt=np.zeros([n_age])
    [nz,nx,ny]=pdl_cube.shape
    temp_a=np.zeros([nx,ny])
    mass_age=np.zeros([nx,ny])
    pdl_cube[np.isnan(pdl_cube)]=1
    
    pdl_cube_e=np.zeros([n,nx,ny])
    #print name
    for i in range(0, n):
        norm_file=dir1+"/"+"map.CS."+name+"_eNORM_"+str(i)+"_ssp.fits.gz"
        if ptt.exists(norm_file) == True:
            pdl_cube_e[i,:,:] = gdata(norm_file)
        else:
            pdl_cube_e[i,:,:] = np.zeros([nx,ny])
    pdl_cube_e[np.isnan(pdl_cube_e)]=0
    #print AGE
    #sys.exit()
    for i in range(nt-1, n_ssp-1, -1):
        label=hdr['FILE_'+str(i)]
#        print label,n_age,n_met,n_age-i+n_ssp-1,-i+n_ssp,i
        time=label.replace('_NORM_age.fits.gz','')
        time=float_(time.replace('map.CS.'+name+'_',''))
        ages[n_age-i+n_ssp-1]=np.log10(time)+9
    #print np.log10(time)+9, 38-i+156
    #temp=pdl_cube[i,:,:]
    #temp=temp*10.0**(ML[i-156])*pdl_flux*Factor
    #temp_a=temp+temp_a
#    MassT=MassT+np.sum(temp)
#    mas1=np.sum(temp[n1])+mas1
#    mas2=np.sum(temp[n2][n2a])+mas2
#    mas3=np.sum(temp[n3][n3a])+mas3
#    mas4=np.sum(temp[n4][n4a])+mas4
#    mass[38-i+156,0]=np.log10(mas1)
#    mass[38-i+156,1]=np.log10(mas2)
#    mass[38-i+156,2]=np.log10(mas3)
#    mass[38-i+156,3]=np.log10(mas4)
    f2=open(dir3+"/"+name+"_Ensemble.csv","w")
    f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    fL2=open(dir3+"/"+name+"_Ensemble_L.csv","w")
    fL2.write("#  LOG_AGE  N_LIGHTR_1  N_LIGHTR_2  N_LIGHTR_3  N_LIGHTR_4  LOG_LIGHTR_1  LOG_LIGHTR_2  LOG_LIGHTR_3  LOG_LIGHTR_4 \n")
    mass_age_t=np.zeros([n_age,nx,ny])
    mass_age_t_2=np.zeros([n_age,nx,ny])
    mass_age_t_e=np.zeros([n_age,nx,ny])
    light_age_t=np.zeros([n_age,nx,ny])
    light_age_t_2=np.zeros([n_age,nx,ny])
    light_age_t_e=np.zeros([n_age,nx,ny])
    for i in range(0, n_age):
        age_now=a_age[i]
        pdl_age=np.zeros([nx,ny])
        pdl_age_2=np.zeros([nx,ny])
        pdl_age_e=np.zeros([nx,ny])
        pdl_ageL=np.zeros([nx,ny])
        pdl_age_2L=np.zeros([nx,ny])
        pdl_age_eL=np.zeros([nx,ny])
        for j in range(0, n):
            if age_now == AGE[j]:
                #if AGE[j] <= 2:
                pdl_age=pdl_age+pdl_cube[j,:,:]*10.0**(ML[j])*pdl_flux*Factor*10.0**(0.4*pdl_Av)*pdl_mask#*0.25/np.pi#/1.47
                pdl_age_e=pdl_age_e+((pdl_cube_e[j,:,:]/pdl_cube[j,:,:])**2.0+(pdl_flux_e/pdl_flux)**2.0+(np.log(10.0)*0.4*pdl_Av_e)**2.0)*(pdl_cube[j,:,:]*10.0**(ML[j])*pdl_flux*Factor*pdl_mask*10.0**(0.4*pdl_Av))**2.0
                pdl_age_2=pdl_age_2+pdl_cube[j,:,:]*10.0**(ML[j])*pdl_flux*Factor*pdl_mask*10.0**(0.4*pdl_Av)/COR[j]
                pdl_age_2L=pdl_age_2L+pdl_cube[j,:,:]*pdl_flux*Factor*pdl_mask*10.0**(0.4*pdl_Av)/COR[j]
                pdl_ageL=pdl_ageL+pdl_cube[j,:,:]*pdl_flux*Factor*10.0**(0.4*pdl_Av)*pdl_mask
                pdl_age_eL=pdl_age_eL+((pdl_cube_e[j,:,:]/pdl_cube[j,:,:])**2.0+(pdl_flux_e/pdl_flux)**2.0+(np.log(10.0)*0.4*pdl_Av_e)**2.0)*(pdl_cube[j,:,:]*pdl_flux*Factor*pdl_mask*10.0**(0.4*pdl_Av))**2.0
                #else:
                #    pdl_age=pdl_age+pdl_cube[j,:,:]*10.0**(ML[j])*pdl_flux*Factor*pdl_mask*COR[j]
        pdl_age[np.where(np.isfinite(pdl_age) == False)]=0
        pdl_age_2[np.where(np.isfinite(pdl_age_2) == False)]=0
        pdl_age_e[np.where(np.isfinite(pdl_age_e) == False)]=0
        pdl_ageL[np.where(np.isfinite(pdl_ageL) == False)]=0
        pdl_age_2L[np.where(np.isfinite(pdl_age_2L) == False)]=0
        pdl_age_eL[np.where(np.isfinite(pdl_age_eL) == False)]=0
        #pdl_age_e[np.isnan(pdl_age_e)]=0
        for k in range(0, n_age):
            if np.log10(age_now)+9 == ages[k]:
                mass_age_t[k,:,:]=pdl_age
                mass_age_t_2[k,:,:]=pdl_age_2
                mass_age_t_e[k,:,:]=pdl_age_e
                light_age_t[k,:,:]=pdl_ageL
                light_age_t_2[k,:,:]=pdl_age_2L
                light_age_t_e[k,:,:]=pdl_age_eL
    temp5=np.sum(mass_age_t,axis=0)+0.01
#    temp6=np.log10(np.sum(mass_age_t,axis=0)+1.0)#QUITAR
#    wfits(dir3+"/"+name+"mass_tot.fits",temp6,hdr001)#QUITAR
    upvalue=math.ceil(np.log10(np.amax(temp5))/.05)*.05
    if np.isinf(upvalue):
        upvalue=8
    if upvalue-1 <= 6.5:
        lovalue=math.ceil(np.log10(np.amin(temp5))/.05)*.05
        if upvalue-2 > lovalue:
            lovalue=upvalue-2 
    else:
        lovalue=6.5
    mass_temp_total=0
    light_temp_total=0
    for i in range(0, n_age):
        if i == 0:
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
            age_i=0.0
        elif i == n_age-1:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=2.0*10.0**(ages[i])-age_i
        else:
            age_i=10.0**((ages[i]+ages[i-1])/2.0)
            age_s=10.0**((ages[i]+ages[i+1])/2.0)
        Dt_age=np.abs(age_s-age_i)
        sfdt[i]=Dt_age/1e6
        temp=mass_age_t[i,:,:]
        temp_2=mass_age_t_2[i,:,:]
        temp_e=mass_age_t_e[i,:,:]
        tempL=light_age_t[i,:,:]
        temp_2L=light_age_t_2[i,:,:]
        temp_eL=light_age_t_e[i,:,:]
        #temp[np.where(np.isfinite(temp) == False)]=0
        #temp_e[np.where(np.isfinite(temp_e) == False)]=1
        #temp_2[np.where(np.isfinite(temp_2) == False)]=0
        if i == 0:
            if fits_f == 1:
                [nx,ny]=temp.shape
                MASS_map_cube=np.zeros([n_age,nx,ny])
                MGH_map_cube=np.zeros([n_age,nx,ny])
                SFH_map_cube=np.zeros([n_age,nx,ny])
                LIGHT_map_cube=np.zeros([n_age,nx,ny])
                LGH_map_cube=np.zeros([n_age,nx,ny])
            temp1=temp
            temp1L=tempL
        else:
            temp1=temp1+temp
            temp1L=temp1L+tempL
        if fits_f == 1:
            MASS_map_cube[i,:,:]=temp
            MGH_map_cube[i,:,:]=temp1
            SFH_map_cube[i,:,:]=temp_2/Dt_age
            LIGHT_map_cube[i,:,:]=tempL
            LGH_map_cube[i,:,:]=temp1L
        #if pdf==1:
        #map_plot(temp1,ages[i],pdl_rad,dir=dir_map+"/",pdf=1,title=name,form='pdf',fname=name+'_smap_'+str(i),minval=lovalue,maxval=upvalue)
        MassT=MassT+np.sum(temp)
        LighT=LighT+np.sum(tempL)
        for ii in range(0, nr-1):
           # print ind[ii],inda[ii],ii
#            print temp[ind[ii]]
#            print len(temp[ind[ii]]),np.amax(inda[ii])
            Dt_mass=np.sum(temp[ind[ii]][inda[ii]])
            Dt_mass_e=np.sum(temp_e[ind[ii]][inda[ii]])
            Dt_mass_2=np.sum(temp_2[ind[ii]][inda[ii]])
            Dt_light=np.sum(tempL[ind[ii]][inda[ii]])
            Dt_light_e=np.sum(temp_eL[ind[ii]][inda[ii]])
            Dt_light_2=np.sum(temp_2L[ind[ii]][inda[ii]])
            massN[ii]=Dt_mass+massN[ii]
            massN_e[ii]=Dt_mass_e+massN_e[ii]
            lightN[ii]=Dt_light+lightN[ii]
            lightN_e[ii]=Dt_light_e+lightN_e[ii]
        #mas2=np.sum(temp[n2][n2a])+mas2
        #mas3=np.sum(temp[n3][n3a])+mas3
        #mas4=np.sum(temp[n4][n4a])+mas4
#            print temp[ind[ii]][inda[ii]]
            mass[i,ii]=np.log10(massN[ii])
            mass_e[i,ii]=massN_e[ii]
            light[i,ii]=np.log10(lightN[ii])
            light_e[i,ii]=lightN_e[ii]
            sfrt[i,ii]=Dt_mass_2/Dt_age#/massN[ii]#/(np.pi*(rx[ii+1]**2.0-rx[ii]**2.0)*rad**2.0)#/(float_(len(temp[ind[ii]][inda[ii]]))*(0.5*DA)**2.0)
        #mass[i,1]=np.log10(mas2)
        #mass[i,2]=np.log10(mas3)
        #mass[i,3]=np.log10(mas4)
    mass_temp_total=np.log10(np.sum(10**mass[nt-n_ssp-1,:]))
    light_temp_total=np.log10(np.sum(10**light[nt-n_ssp-1,:]))
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    LighT=np.log10(LighT)
    LighT=light_temp_total
    if fits_f == 1:
        #if ptt.exists() == True:
        h1=pyf.PrimaryHDU(MASS_map_cube)#.header
        h2=pyf.PrimaryHDU(SFH_map_cube)#.header
        h3=pyf.PrimaryHDU(MGH_map_cube)
        h4=pyf.PrimaryHDU(LIGHT_map_cube)
        h5=pyf.PrimaryHDU(LGH_map_cube)
        h=h1.header
        h["NAXIS"]=3
        h["NAXIS3"]=n_age
        h["NAXIS1"]=nx
        h["NAXIS2"]=ny
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"MASS_maps_"+name+".fits",hlist)
        #if ptt.exists() == True:
        h=h2.header
        h["NAXIS"]=3
        h["NAXIS3"]=n_age
        h["NAXIS1"]=nx
        h["NAXIS2"]=ny
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"SFH_maps_"+name+".fits",hlist) 
        #if ptt.exists() == True:   
        h=h3.header
        h["NAXIS"]=3
        h["NAXIS3"]=n_age
        h["NAXIS1"]=nx
        h["NAXIS2"]=ny
        hlist=pyf.HDUList([h3])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"MGH_maps_"+name+".fits",hlist)
        h=h4.header
        h["NAXIS"]=3
        h["NAXIS3"]=n_age
        h["NAXIS1"]=nx
        h["NAXIS2"]=ny
        hlist=pyf.HDUList([h4])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"LIGHT_maps_"+name+".fits",hlist)
        h=h5.header
        h["NAXIS"]=3
        h["NAXIS3"]=n_age
        h["NAXIS1"]=nx
        h["NAXIS2"]=ny
        hlist=pyf.HDUList([h5])
        hlist.update_extend()
        wfits_ext(dir_map+"/"+"LGH_maps_"+name+".fits",hlist)
    #print MassT,name
    if ptt.exists(dir1+"/"+SN_file) == True:
        SN=np.zeros(nr-1)
        sn_l=''
        for ii in range(0, nr-1):
            SN[ii]=np.average(pdl_SN[ind[ii]][inda[ii]])
            sn_l=sn_l+' , '+str(SN[ii])
        if fi != 0:
            fi.write(name+sn_l+' \n')
    if ptt.exists(dir1+"/"+Ha_file) == True:
        Ha=np.zeros(nr-1)
        ha_l=''
        for ii in range(0, nr-1):
            Ha[ii]=np.sum(pdl_ha[ind[ii]][inda[ii]]*10.0**(0.4*pdl_Av[ind[ii]][inda[ii]]))*(L*1e-16)
            ha_l=ha_l+' , '+str(Ha[ii])
        if fii != 0:
            fii.write(name+ha_l+' \n')
    else:
        ha_l=''
        for ii in range(0, nr-1):
            ha_l=ha_l+' , '+str(-100)
        if fii != 0:
            fii.write(name+ha_l+' \n')
    #print Ha,(L*1e-16)
    #sys.exit(0)
    mass_n=10**(10**(mass-mass[nt-n_ssp-1,:]))
    mass_n_e=np.sqrt((10**(mass-mass[nt-n_ssp-1,:]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[nt-n_ssp-1,:]/10**(2.0*mass[nt-n_ssp-1,:]))))
    light_n=10**(10**(light-light[nt-n_ssp-1,:]))
    light_n_e=np.sqrt((10**(light-light[nt-n_ssp-1,:]))**2.0*((light_e/10**(2.0*light))+(light_e[nt-n_ssp-1,:]/10**(2.0*light[nt-n_ssp-1,:]))))
    #mass_n=10**(mass-mass[nt-156-1,:])
    #mass_n=(mass-mass[nt-156-1,:])
    for i in range(0, n_age):
        #print ages[i],a_age[i],"test_ages"
        line=''
        line=line+str(ages[i])
        for ii in range(0, nr-1):
            line=line+';'+str(mass_n[i,ii])
        for ii in range(0, nr-1):
            line=line+';'+str(mass[i,ii])
        for ii in range(0, nr-1):
            line=line+';'+str(sfrt[i,ii])
        for ii in range(0, nr-1):
            line=line+';'+str(mass_n_e[i,ii])
        line=line+';'+str(sfdt[i])
        line=line+' \n'
        f2.write(line)
        lineL=''
        lineL=lineL+str(ages[i])
        for ii in range(0, nr-1):
            lineL=lineL+';'+str(light_n[i,ii])
        for ii in range(0, nr-1):
            lineL=lineL+';'+str(light[i,ii])
        for ii in range(0, nr-1):
            lineL=lineL+';'+str(sfrt[i,ii])
        for ii in range(0, nr-1):
            lineL=lineL+';'+str(light_n_e[i,ii])
        lineL=lineL+';'+str(sfdt[i])
        lineL=lineL+' \n'
        fL2.write(lineL)
    #if not pdf == 0:
        #dev=dir3+"/"+name+"_Relative_Mass.pdf"
        ##if pdf == 1:
        #    #matplotlib.use('Agg')
        #import matplotlib.pyplot as plt
        ##plt.axis([8, 10.5, 0, 10])
        #plt.xlabel("$log_{10}(time/yr)$",fontsize=14)
        #plt.ylabel("$10^{M(t)/M_{0}}$",fontsize=14)
        #plt.title(name+' $\log M_{tot}='+('%7.2f' % MassT)+'$',fontsize=15)
        ##plt.semilogx('log')
        #for ii in range(0, nr-1):
        #    plt.plot(ages,mass_n[:,ii],label='$'+('%6.1f' % rx[ii])+'R_e<R<'+('%6.1f' % rx[ii+1])+'R_e$')
        ##plt.plot(ages,mass_n[:,1],label='$'+('%6.1f' % r1)+'<R<'+('%6.1f' % r2)+'R_e$')
        ##plt.plot(ages,mass_n[:,2],label='$'+('%6.1f' % r2)+'<R<'+('%6.1f' % r3)+'R_e$')
        ##plt.plot(ages,mass_n[:,3],label='$'+('%6.1f' % r3)+'<R<'+('%6.1f' % r4)+'R_e$')
        #plt.legend(loc=3)
        #if pdf == 1:
        #    plt.savefig(dev,dpi = 1000)
        #else:
        #    plt.show()
        #plt.close()
    if not pdf == 0:
        dev=dir3+'/'+name+"_Relative_Mass2.pdf"
        #if pdf == 1:
        #    matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(6,5.5))
        ax.set_xlabel("$log_{10}(time/yr)$",fontsize=14)
        ax.set_ylabel("$M(t)/M_{0}$",fontsize=14)
        #MassT=10.32
        ax.set_title(name+' $\log M_{tot}='+('%7.2f' % MassT)+'$',fontsize=15)
        ax.set_xlim(8.6,10.1)
            #ax.set_ylim(0,12)
        ax.set_ylim(func_plot(np.log10(1.78),ftype=1),func_plot(np.log10(12),ftype=1))
        for ii in range(0, nr-1):
            plt.plot(ages,func_plot(np.log10(mass_n[:,ii]),ftype=1),label='$'+('%6.1f' % rx[ii])+'R_e<R<'+('%6.1f' % rx[ii+1])+'R_e$')
        plt.legend(loc=3)
        plt.plot(np.arange(0,20,.1),np.ones(200)*func_plot(0.95,ftype=1),'--',color='black')
        plt.plot(np.arange(0,20,.1),np.ones(200)*func_plot(0.50,ftype=1),'--',color='green')
        fig.canvas.draw()
        labels = [item.get_text() for item in ax.get_yticklabels()] 
        for i in range(0, len(labels)):
            labels[i]=labels[i].replace(u'\u2212','-')
        for i in range(0, len(labels)):
            if labels[i] != u'':
                if float_(labels[i]) == 0:
                    labels[i]=u'%3.2f' % 10**(0)
                else:
                    labels[i]=u'%3.2f' % 10**(float_(labels[i]))
        ax.set_yticklabels(labels)
        if pdf == 1:
            fig.tight_layout()
            plt.savefig(dev)#,dpi = 1000)
        else:
            plt.show()
        plt.close()
    f2.close()
    fL2.close()
    fo.write(name+" "+str(MassT)+" "+str(redshift)+" ")
    
