import matplotlib
import sys
import numpy as np
import cosmolopy.distance as cd
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import time
import os.path as ptt
<<<<<<< HEAD
import h5py
import illustris_python as il
=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
from scipy.interpolate.interpolate import interp1d
import pyfits as pyf
from pyfits import writeto as wfit
from pyfits import update as ufits 
import scipy.signal as sign
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from pyfits import getdata as gdata
from astropy import wcs
from astropy.io import fits
from pyfits import getheader as ghead
from scipy import stats
from astropy.units.equivalencies import spectral
from astropy.io.fits.hdu.compressed import DITHER_SEED_CHECKSUM
from scipy.ndimage.filters import convolve1d
import warnings
<<<<<<< HEAD
from illustris_python import lhalotree
matplotlib.use('Agg')
=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
warnings.filterwarnings("ignore")

def galfit_param(name,dir='',dir2='',band='',psf=1.2,ra_c=0,dec_c=0,dx=100,dy=100,repro=0):
    dirn=dir.replace(" ","\ ")
    dir2n=dir2.replace(" ","\ ")
    [data, head]=gdata(dir+name+band+".fits", 0, header=True)
    head_n=head
    dpx_ra=head["CD1_1"]*3600.
    dpx_de=head["CD2_1"]*3600.
    dpy_ra=head["CD1_2"]*3600.
    dpy_de=head["CD2_2"]*3600.
    hdulist = fits.open(dir+name+band+".fits")
#    print hdulist[0].header
    w = wcs.WCS(hdulist[0].header)
    pixcrd = np.array([[ra_c, dec_c],[ra_c, dec_c]])
#    print w
#    print pixcrd
    x, y = w.wcs_world2pix(pixcrd, 1)[0]
    x=int(x)
    y=int(y)
    [nxi,nyi]=data.shape
    yo=y-dy
    if yo < 0:
        yo=0
        dy=y
        dx=dy
    yf=y+dy
    if yf > nxi:
        yf=nxi-1
        dy=yf-y
        dx=dy
    xo=x-dx
    if xo < 0:
        xo=0
        dx=x
        dy=dx
    xf=x+dx
    if xf > nyi:
        xf=nyi-1
        dx=xf-x
        dy=dx
    yo=y-dy
    yf=y+dy
    xo=x-dx
    xf=x+dx
    data2=data[yo:yf,xo:xf]/1e-20
    [nx,ny]=data2.shape
    head["NAXIS1"]=nx
    head["NAXIS2"]=ny
    wfits(name+band+"_extr.fits", data2, head)
    dpx=np.sqrt(dpx_ra**2.+dpx_de**2.)
    dpy=np.sqrt(dpy_ra**2.+dpy_de**2.)#0.396
    sig=psf/dpx
    X = np.arange(-nx/2, nx/2, 1)
    Y = np.arange(-nx/2, nx/2, 1)
    X, Y = np.meshgrid(X, Y)
    R = np.exp(-1./2.*(X**2 + Y**2)/sig**2)/(2.*np.pi*sig**2)
    wfits(name+band+"_PSF.fits", R, head)
    #if ptt.exists(dir+name+'_'+band+'_res.fits') == False:
    if ptt.exists(dir+name+band+"_galfit.feedme") == False:# or repro == 1:
        f=open("galfit.feedme","r")
        f2=open(name+band+"_galfit.feedme","w")
        count=0
        for line in f:
            if count == 3:
                line=line.replace('test.fits',name+band+'_extr.fits')
            if count == 4:
                line=line.replace('imgblock.fits',name+band+'_mod.fits')
            if count == 6:
                line=line.replace('test_PSF.fits',name+band+'_PSF.fits')
            if count == 10:
                line=line.replace('1    200   1    200','1    '+str(nx)+'   1    '+str(ny))
            if count == 11:
                line=line.replace('200    200',str(nx)+'    '+str(ny))
            if count == 13:
                line=line.replace('0.038  0.038',str(dpx)+'  '+str(dpy))
            if count == 34:
                line=line.replace('100  100',str(nx/2)+'  '+str(ny/2))
            f2.write(line)
            count=count+1
        f2.close()
        f.close()
    else:
        sycall("cp "+dirn+name+band+"_galfit.feedme "+name+band+"_galfit.feedme")
    #if ptt.exists(dir+name+'_'+band+'_res.fits') == False:
    #if name == "manga-8317-1901":
    #    repro == 1
    not_fit=0
    if ptt.exists(dir+name+band+"_galfit.03") == False or repro == 1:
        sycall("./galfit "+name+band+"_galfit.feedme")
        sycall("mv galfit.01 "+name+band+"_galfit.01")
        if ptt.exists(dir2+name+band+"_galfit.01") == False:
           # not_fit=1
            sycall("cp "+name+band+"_galfit.feedme "+dir2n+name+band+"_galfit.01")
        sycall("./galfit "+name+band+"_galfit.01")
        sycall("mv galfit.01 "+name+band+"_galfit.02")
        if ptt.exists(dir2+name+band+"_galfit.02") == False:
          #  not_fit=1
            sycall("cp "+name+band+"_galfit.01 "+dir2n+name+band+"_galfit.02")         
        sycall("./galfit "+name+band+"_galfit.02")
        sycall("mv galfit.01 "+name+band+"_galfit.03")
        if ptt.exists(dir2+name+band+"_galfit.03") == False:
         #   not_fit=1
            sycall("cp "+name+band+"_galfit.02 "+dir2n+name+band+"_galfit.03")
    else:
        sycall("cp "+dirn+name+band+"_galfit.03 "+dir2n+name+band+"_galfit.03")
        sycall("cp "+dirn+name+band+"_galfit.03 "+name+band+"_galfit.03")
    if ptt.exists(dir2+name+band+"_galfit.01") == False:
        sycall("cp "+name+band+"_galfit.feedme "+dir2n+name+band+"_galfit.03")
        #not_fit=1
    f=open(name+band+"_galfit.03","r")
    f2=open(name+band+"_galfit.04","w")
    count=0
    for line in f:
        if count == 19:
            line=line.replace('P) 0','P) 1')
        if count == 42:
            datl=line.split(" ")
            datl=filter(None,datl)
            e_rad=float_(datl[1])
        if count == 43:
            datl=line.split(" ")
            datl=filter(None,datl)
            n_ser=float_(datl[1])
        if count == 47:
            datl=line.split(" ")
            datl=filter(None,datl)
            ab=float_(datl[1])
        if count == 48:
            datl=line.split(" ")
            datl=filter(None,datl)
            PA=float_(datl[1])
        f2.write(line)
        count=count+1
    f2.close()
    f.close()
    e_rad_f=e_rad*(dpx+dpy)/2.
    if not_fit == 0:
        if ptt.exists(dir+name+band+"_rad.fits") == False:
            sycall("./galfit "+name+band+"_galfit.04")
            sycall("python sub_fits.py "+name+band+'_extr.fits  '+name+band+'_mod.fits  '+name+band+'_res.fits')
            nx_n=head_n["NAXIS1"]
            ny_n=head_n["NAXIS2"]
            radis=np.zeros([nx_n,ny_n])
            for i in range(0, nx_n):
                for j in range(0, ny_n):
                    x_n=-((i+1)-nx_n/2-1)*np.sin(PA*np.pi/180.)+(-(j)+ny_n/2+1)*np.cos(PA*np.pi/180.)
                    y_n=-((i+1)-nx_n/2-1)*np.cos(PA*np.pi/180.)-(-(j)+ny_n/2+1)*np.sin(PA*np.pi/180.)
                    r_n=np.sqrt((y_n)**2.0+(x_n/ab)**2.0)/e_rad
                    radis[i,j]=r_n
            wfits(dir+name+band+"_rad.fits", radis, head_n)
    return [e_rad_f,e_rad,n_ser,ab,PA,not_fit]

def func_plot(x,ftype=0):
    if ftype == 0:
        y=10**(x)
    if ftype == 1:
        y=np.log10(x)
    if ftype == 2:
        y=x
    return y

def imag_conv(outf,x,y,z,input_val,fo,fi,index=[[]],ages_f=[0],dir_o='',photo=1,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0]):
    n_ages=1
    no_nan=0
    if photo != 1:
        n_ages=len(ages_f)
    if n_ages > 1:
        outf=outf+'_mass'
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    if photo == 1:
        imag=np.ones([nl,nl])*(10**(-0.4*22.)*(dap*3600.)**2.)
    else:
        if n_ages == 1:
            imag=np.zeros([nl,nl])
        else:
            imag=np.zeros([n_ages,nl,nl])
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    xi=x
    yi=y
    zi=z
    input_vali=input_val
    for k in range(0, n_ages):
        if n_ages > 1:
            x=xi[index[k]]
            y=yi[index[k]]
            z=zi[index[k]]
            input_val=input_vali[index[k]]
        rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
        reds=reds_cos(rad/1e3)
        radA=rad/(1+reds)
        radL=rad*(1+reds)
        phi=np.arcsin(x/radA)
        the=np.arcsin(y/(radA*np.cos(phi)))
        the=the*180/np.pi
        phi=phi*180/np.pi
        for i in range(0, nl):
            dp0=oap+dap*i
            dp1=oap+dap*(i+1)
            yima[i]=(dp0+dp1)/2.
            xima[i]=(dp0+dp1)/2.
            #print i
            for j in range(0, nl):
                dt0=oap+dap*j
                dt1=oap+dap*(j+1)
                n=np.where((the >= dt0) & (the < dt1) & (phi >= dp0) & (phi < dp1))[0]
                #print i#, dt0, dt1, dp0, dp1, np.amax(the), np.amin(the), np.amax(phi), np.amin(phi)
#                print modes(x),modes(y),modes(z)
#                import matplotlib.pyplot as plt
#                n, bins, patches = plt.hist(x, 50, normed=1, histtype='stepfilled')
#                plt.show()
#                sys.exit(0)
                if len(n) > 0:
                    if photo == 1:
                        flux=10.**(-0.4*input_val[n])
                        dis=radL[n]*1e3
                        nf=flux*(10./dis)**2.
#                        print nf
#                        print dis
#                        sys.exit(0)
                        val_f=np.sum(nf)+10**(-0.4*22.)*(dap*3600.)**2.
                        imag[j,i]=val_f
                        #print val_f
                    else:
                        nf=input_val[n]*1e10
                        val_f=np.sum(nf)
                        if n_ages == 1:
                            imag[j,i]=val_f
                        else:
                            imag[k,j,i]=val_f
#    print np.isinf(imag)
#    print np.isnan(imag)
    [X,Y]=np.meshgrid(xima,yima)
    dv=sig/3600./dap
    di=sig/3600.
    PSF=Gaussian2DKernel(stddev=dv)#,x_size=nl,y_size=nl)
    #print PSF
    if n_ages == 1:
        print imag.shape
        print PSF.shape
        imag_F=convolve(imag, PSF)#, mode='full')#,  boundary='symm')
        #imag_F=imag
    else:
        imag_F=imag
        for k in range(0, n_ages):
            imag_F[k,:,:]=convolve(imag[k,:,:], PSF)
    h=pyf.PrimaryHDU().header
    if n_ages == 1:
        h["NAXIS"]=2
    else:
        h["NAXIS"]=3
        h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="OUTPUT auto_ssp_elines_rnd.pl FITs"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    if n_ages > 1:
        h['CDELT3']=ages_f[1]-ages_f[0]
        h['CRPIX3']=1
        h['CRVAL3']=ages_f[0]+9
        h['CUNIT3']='log(age [yrs])'
#    del h['RADECSYS']
    h['RADECSYSa']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    #print red_0
    out_fit=dir_o+outf+'.fits'
    wfits(out_fit,imag_F,h)
    if n_ages > 1:
        name_pdl=outf.replace('_mass','_rad')
        name=outf.replace('_mass','')
        [pdl_rad, hd]=gdata(dir_o+name_pdl+".fits", 0, header=True)
        ind=[]
        inda=[]
        nr=len(rx)
        for ii in range(0, nr-1):
            nt=np.where(pdl_rad< rx[ii+1])
            nta=np.where(pdl_rad[nt]>= rx[ii])
            ind.extend([nt])
            inda.extend([nta])   
        MassT=0
        massN=np.zeros(nr)
        mass=np.zeros([n_ages,nr])
        ages=np.zeros([n_ages])
        for i in range(n_ages-1, -1, -1):
            temp=imag_F[i,:,:]
            MassT=MassT+np.sum(temp)
            for ii in range(0, nr-1):
                massN[ii]=np.sum(temp[ind[ii]][inda[ii]])+massN[ii]
                mass[n_ages-1-i,ii]=np.log10(massN[ii])
            ages[n_ages-1-i]=ages_f[i]+9
        MassT=np.log10(MassT)
        mass_n=10**(10**(mass-mass[n_ages-1,:]))
        tem=np.sum(mass)
        if np.isnan(tem) or np.isnan(MassT):
            no_nan=1
        f2=open(dir_o+name+"_Ensemble.csv","w")
        f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")        
        for i in range(0, n_ages):
            line=''
            line=line+str(ages[i])
            for ii in range(0, nr-1):
                line=line+';'+str(mass_n[i,ii])
            for ii in range(0, nr-1):
                line=line+';'+str(mass[i,ii])
            line=line+' \n'
            f2.write(line)
        if not pdf == 0:
            dev=dir_o+name+"_Relative_Mass.pdf"
            if pdf == 1:
                matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(6,5.5))
            ax.set_xlabel("$log_{10}(time/yr)$",fontsize=14)
            ax.set_ylabel("$M(t)/M_{0}$",fontsize=14)
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
        if no_nan == 0:
            fo.write(name+" "+str(MassT)+" "+str(red_0)+" ")
        fi.write(name+" "+str(MassT)+" \n")
    return [imag_F,X,Y,oap,no_nan]
   
def sycall(comand):
    import os
    linp=comand
    os.system(comand)

def reds_cos(dis):
    red=np.arange(0,3,.01)
    cosmo = {'omega_M_0' : 0.2726, 'omega_lambda_0' : 0.7274, 'h' : 0.704}
    cosmo = cd.set_omega_k_0(cosmo)
    dist=cd.comoving_distance(red, **cosmo)
    z=interp1d(dist, red,kind='linear',bounds_error=False)(dis)
    return z

def sycallo(comand):
    import os
    out=os.popen(comand, 'r')
    line=out.readline()
    return line

def modes(x,nit=3,n_s=2):
    xt=x
    for i in range(0, nit):
        xo=np.average(xt)
        sx=np.std(xt)
        n=np.where((xt <xo+n_s*sx) & (xt > xo-n_s*sx))[0]
        xt=xt[n]
    xm=np.average(xt)
    return xm

def wfits(name, data, hdr):
    if ptt.exists(name) == False:
        wfit(name,data,hdr)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        wfit(name,data,hdr)

def wfits_ext(name,hlist):
    if ptt.exists(name) == False:
        hlist.writeto(name)
    else:
        name1=name.replace("\ "," ")
        name1=name1.replace(" ","\ ")
        sycall("rm "+name1)
        hlist.writeto(name)

def A_l(Rv,l):
    l=l/10000.; #Amstrongs to Microns
    x=1.0/l
    Arat=np.zeros(len(x))
    for i in range(0, len(x)):
        if x[i] > 1.1 and x[i] <= 3.3:
            y=(x[i]-1.82)
            ax=1+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7
            bx=1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
        if x[i] <= 1.1 and x[i] > 0.3:
            ax=0.574*x[i]**1.61
            bx=-0.527*x[i]**1.61
        if x[i] > 3.3 and x[i] <= 8.0:
            if x[i] > 5.9 and x[i] <= 8.0:
                Fa=-0.04473*(x[i]-5.9)**2.0-0.009779*(x[i]-5.9)**3.0
                Fb=0.2130*(x[i]-5.9)**2.0+0.1207*(x[i]-5.9)**3.0
            else:
                Fa=0.0
                Fb=0.0
            ax=1.752-0.316*x[i]-0.104/((x[i]-4.67)**2.0+0.341)+Fa
            bx=-3.090+1.825*x[i]+1.206/((x[i]-4.62)**2.0+0.263)+Fb
        if x[i] > 8.0:
            ax=-1.073-0.628*(x[i]-8.0)+0.137*(x[i]-8.0)**2.0-0.070*(x[i]-8.0)**3.0
            bx=13.670+4.257*(x[i]-8.0)-0.420*(x[i]-8.0)**2.0+0.374*(x[i]-8.0)**3.0
        val=ax+bx/Rv
        if val < 0:
            val=0
        Arat[i]=val
    return Arat

def ages_redefinition_l(ages,ages_ssp):
    ages_l=ages
    ages_f=np.zeros(len(ages_l))
    n_ages=len(ages_ssp)
    #ind=[]
    for i in range(0, n_ages):
        if i < n_ages-1:
            dage=(ages_ssp[i+1]+ages_ssp[i])/2.
        else:
            dage=ages_ssp[i]+1.0
        if i == 0:
            age1=0
        else:
            age1=age2
        age2=dage
        nt=np.where((ages_l > age1) & (ages_l <= age2))#[0]
        ages_f[nt]=ages_ssp[i]
        #ind.extend([nt[0]])
    return ages_f

def ages_definition_l(ages,ages_ssp):
    ages_l=ages
    n_ages=len(ages_ssp)
    ind=[]
    for i in range(0, n_ages):
        if i < n_ages-1:
            dage=(ages_ssp[i+1]+ages_ssp[i])/2.
        else:
            dage=ages_ssp[i]+1.0
        if i == 0:
            age1=0
        else:
            age1=age2
        age2=dage
        nt=np.where((ages_l > age1) & (ages_l <= age2))#[0]
        ind.extend([nt[0]])
    return ind

def met_definition_l(met,met_ssp):
    met_l=met
    n_met=len(met_ssp)
    ind=[]
    for i in range(0, n_met):
        if i < n_met-1:
            dage=(met_ssp[i+1]+met_ssp[i])/2.
        else:
            dage=met_ssp[i]+1.0
        if i == 0:
            met1=0
        else:
            met1=met2
        met2=dage
        nt=np.where((met_l > met1) & (met_l <= met2))#[0]
        ind.extend([nt[0]])
    return ind

def val_definition_l(val,val_ssp):
    val_l=val
    n_val=len(val_ssp)
    ind=[]
    for i in range(0, n_val):
        if i < n_val-1:
            dval=(val_ssp[i+1]+val_ssp[i])/2.
        else:
            dval=val_ssp[i]+1.0
        if i == 0:
            val1=0
        else:
            val1=val2
        val2=dval
        nt=np.where((val_l > val1) & (val_l <= val2))
        ind.extend([nt[0]])
    return ind

def ages_definition(ages,ages_m=14.2,ages_mi=0.001,n_ages=50):
    ages_l=np.log10(ages)
    age_m=np.log10(ages_m)
    age_mi=np.log10(ages_mi)
    dage=(age_m-age_mi)/float_(n_ages)
    age_f=np.zeros(n_ages)
    ind=[]
    for i in range(0, n_ages):
        age1=age_mi+i*dage
        age2=age_mi+(i+1)*dage
        age_0=(age1+age2)/2.
        #print age_0
        nt=np.where((ages_l > age1) & (ages_l <= age2))#[0]
        age_f[i]=age_0
        ind.extend([nt[0]])
    return [age_f,ind]

<<<<<<< HEAD
def integ_flux_2(dr,x,y,type=0,sig=3,xo=0.0,yo=0.0,dt=0.01):
    nt=int(dr/dt)
    xt=-nt/2*dt+x
    yt=-nt/2*dt+y
    xi=xt
    xf=xt
    Ft=0.0
    for i in range(0, nt):
        xi=xf
        xf=xf+dt
        yi=yt
        yf=yt
        for j in range(0, nt):
            yi=yf
            yf=yf+dt
            if type == 0:
                Rt=np.sqrt(((xi+xf)/2.0-x)**2.0+((yi+yf)/2.0-y)**2.0)
                if Rt <= dr/2.0 :
                    Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0
            else:
                Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0
    print x,y
    return Ft


=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
def integ_flux(dr,x,y,type=0,sig=3,xo=0.0,yo=0.0,dt=0.01):
    nt=int(dr/dt)
    xt=-nt/2*dt+x
    yt=-nt/2*dt+y
    xi=xt
    xf=xt
    Ft=0.0
    for i in range(0, nt):
        xi=xf
        xf=xf+dt
        yi=yt
        yf=yt
        for j in range(0, nt):
            yi=yf
            yf=yf+dt
            if type == 0:
<<<<<<< HEAD
                Rt=np.sqrt(((xi+xf)/2.0-x)**2+((yi+yf)/2.0-y)**2.)
=======
                Rt=np.sqrt( ((xi+xf)/2.0-x)**2+((yi+yf)/2.0-y)**2. )
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
                if Rt <= dr/2.0 :
                    Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0+1*dt**2.0
            else:
                Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0+1*dt**2.0
    print x,y
    return Ft
<<<<<<< HEAD

def don(x,y,r1,r2,type=1):
    from scipy.interpolate.interpolate import interp1d
    tht=np.arange(1000)/999.0*np.pi-np.pi/2.0
    x1=(r2+r1)*np.sin(tht)
    y1=(r2+r1)*np.cos(tht)
    y2=(r2-r1)*np.cos(tht)
    x2=(r2-r1)*np.sin(tht)
    if type == 1:
        y3=interp1d(x2,y2,bounds_error=False,fill_value=0)(x1)
    else:
        y3=np.zeros(len(x1))
    x1=x1+x
    y4=-y3
    y5=-y1
    y1=y1+y
    y3=y3+y
    y5=y5+y
    y4=y4+y
    return x1,y1,y3,y4,y5

def plot_don(x,y,r1=7.4,r2=15.0,color='blue',apl='0.1',type=1):
    import matplotlib.pyplot as plt
    for i in range(0, len(x)):
        x1,y1,y2,y3,y4=don(x[i],y[i],r1,r2,type=type)
        plt.fill_between(x1,y1,y2,alpha=apl,color=color,lw=0)
        plt.fill_between(x1,y3,y4,alpha=apl,color=color,lw=0)

def lst_c(mjd=56000.0,tai=np.inf,lng=-105.820417):
    if np.isfinite(tai) == True:
        jd=2400000.5+tai/(24.0*3600.0)
    else:
        jd=mjd+2400000.5
    c = [280.46061837,360.98564736629,0.000387933,38710000.0]
    jd2000=2451545.0
    t0 = jd - jd2000
    t = t0/36525.0
    theta=c[0]+(c[1]*t0)+t**2.0*(c[2]-t/c[3])
    
    lst=(theta+lng+360.0)/15.0
    if lst < 0.0:
        lst = 24.0+(lst % 24.0)
    lst = lst % 24.0
    return lst

def paralactic_angle(ha,dec=0.0,phi=32.78):
    p=np.arctan2(np.sin(ha/24.0*360.0*np.pi/180.0),(np.tan(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)-np.sin(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)))*180.0/np.pi
    p=p % 360.0
    return p    

def airmas(ha,dec=0.0,phi=32.78):
    cosz=np.sin(phi*np.pi/180.0)*np.sin(dec*np.pi/180.0)+np.cos(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)
    M=(1.002432*cosz**2.0+0.148386*cosz+0.0096467)/(cosz**3.0+0.149864*cosz**2.0+0.0102963*cosz+0.000303978)# Airmas formula from Young 1994
    return M
   
def refrac_dif(wave,ha,dec=0.0,phi=32.78,lo=5070.0,T=7.0,P=600.0,f=8.0,vapor=False):
    cosz=np.sin(phi*np.pi/180.0)*np.sin(dec*np.pi/180.0)+np.cos(phi*np.pi/180.0)*np.cos(dec*np.pi/180.0)*np.cos(ha/24.0*360.0*np.pi/180.0)
    R=64.328+29498.1/(146.0-(1.0/(wave*1e-4))**2.0)+255.4/(41.0-(1.0/(wave*1e-4))**2.0)
    R=R*P*(1.0+P*(1.049-0.0157*T)*1e-6)/(720.883*(1.0+0.003661*T))
    R1=64.328+29498.1/(146.0-(1.0/(lo*1e-4))**2.0)+255.4/(41.0-(1.0/(lo*1e-4))**2.0)
    R1=R1*P*(1.0+P*(1.049-0.0157*T)*1e-6)/(720.883*(1.0+0.003661*T))
    if vapor == True:
        R=R*(0.0624-0.000680/(wave*1e-4)**2.0)/(1.0+0.003661*T)*f
        R1=R1*(0.0624-0.000680/(lo*1e-4)**2.0)/(1.0+0.003661*T)*f
    R=R-R1
    z=np.arccos(cosz)
    R=R*np.tan(z)
    R=R/1e6*3600*180/np.pi
    return R
    
def flat_r(x):
    #x=np.arange(5300.0,10352.0,0.5)
    y1=1.5*np.sqrt(np.exp(-((x-9600.0)/300.0)**2.0/2.0))
    y2=0.01*np.sqrt(np.exp(-((x-9000.0)/100.0)**2.0/2.0))
    y3=1.0/(1.0+np.exp(((x-9000.0)/100.0)))/(1.0+np.exp(((7000.0-x)/100.0)))
    y4=0.3*np.sqrt(np.exp(-((x-7100.0)/200.0)**2.0/2.0))
    y5=0.05*np.sqrt(np.exp(-((x-9200.0)/550.0)**2.0/2.0))
    y6=0.3*np.sqrt(np.exp(-((x-9100.0)/150.0)**2.0/2.0))
    y7=0.2*np.exp(-((x-6350.0)/200.0)**2.0/2.0)
    y8=0.03*np.exp(-((x-6200.0)/50.0)**2.0/2.0)
    y9=0.03*np.exp(-((x-6000.0)/50.0)**2.0/2.0)
    y10=0.03*np.exp(-((x-5850.0)/50.0)**2.0/2.0)
    y11=ran.randn(len(x))*0.005
    y12=1.0/(1.0+np.exp((6600.0-x)/100.0))
    y13=0.25/(1.0+np.exp(((x-10000.0)/100.0)))/(1.0+np.exp(((9000.0-x)/100.0)))
    y14=-0.03*np.exp(-((x-9300.0)/50.0)**2.0/2.0)
    y15=1.0/(1.0+np.exp((x-9800.0)/200.0))
    y=(y1+y2+y3+y4+y5+y6+y7+y8+y9+y10+y11*y12+y13+y14)*y15
    return y
        
def flat_b(x):
    #x=np.arange(5300.0,10352.0,0.5)
    #x=np.arange(3500.0,6510.0,0.5)
    y1=1.4*np.exp(-((x-4700)/520.0)**2.0/2.0)
    y2=0.9*np.exp(-((x-5600)/180.0)**2.0/2.0)
    y3=0.1*np.exp(-((x-5100)/120.0)**2.0/2.0)
    y4=0.1*np.exp(-((x-6080)/100.0)**2.0/2.0)
    y5=0.06*np.exp(-((x-6180)/100.0)**2.0/2.0)
    y6=ran.randn(len(x))*0.005
    y7=-0.1*np.exp(-((x-4000)/30.0)**2.0/2.0)
    y8=-0.05*np.exp(-((x-4500)/30.0)**2.0/2.0)
    y9=-0.05*np.exp(-((x-4700)/10.0)**2.0/2.0)
    y10=1./(1.+np.exp((3830-x)/120.0))
    y11=1./(1.+np.exp((x-6210)/180.0))
    y12=1.45
    y=(y1+y2+y3+y4+y5+y6+y7+y8+y9)*y10*y11*y12        
    return y

def corr_f(wave,dir_tem='libs'):
    f=open(dir_tem+'/corf_0.txt','r')
    corfi=[]
    wavei=[]
    for line in f:
        data=line.replace('\n','').split(' ')
        data=filter(None,data)
        wavei.extend([np.float(data[0])])
        corfi.extend([np.float(data[1])])
    f.close()
    corfi=np.array(corfi)
    wavei=np.array(wavei)
    corr=interp1d(wavei,corfi,bounds_error=False,fill_value=1.)(wave)
    return corr

def flux_distort_corrvec(wavevec, xx, yy, fibid):
    if fibid < 500:
        specid=1
    else:
        specid=2

    coeff0 =  -0.022340535
    coeff1 =  -0.075340805
    coeff2 =  0.00044807866
    coeff3 =  0.0046223343
    coeff4 =  0.0043642268
    coeff5 =  0.022493345
    coeff6 =  0.028808753
    coeff7 =  -0.021108005
    coeff8 =  0.026092067
    coeff9 =  0.016979068
    coeff10 =  -0.013467023
    coeff11 =  0.060558261
    coeff12 =  0.14761290

    lam0 = 5070.0
    lwave2 = 1.0 - (lam0/wavevec)**2.0 # Basically normalized to [0.5,2.0]
    xx = xx / 325.0 # This 325 is mm on the plate
    yy = yy / 325.0 # This 325 is mm on the plate
    if specid == 1:
        cvec=(1.0+coeff0)*np.exp(coeff2*xx+coeff3*yy+coeff4*xx*yy+coeff5*xx**2.0+coeff6*yy**2.0+coeff7*xx*lwave2+coeff8*yy*lwave2+coeff9*lwave2+coeff11*lwave2**2.0)
    else:
        cvec=(1.0+coeff1)*np.exp(coeff2*xx+coeff3*yy+coeff4*xx*yy+coeff5*xx**2.0+coeff6*yy**2.0+coeff7*xx*lwave2+coeff8*yy*lwave2+coeff10*lwave2+coeff12*lwave2**2.0)
#        + coeff[11] * xx*yy * lwave2 $
#        + coeff[12] * xx^2 * lwave2 $
#        + coeff[13] * yy^2 * lwave2)
    cvec=1.0/cvec
    cvec=cvec/corr_f(wavevec)
    return cvec
        
def response_v(x,fac=4.0,dir_tem='libs'):
    #x=np.arange(5000)*1.3+3600.0
    #y1=50.0*(0.5-0.5*np.sin((x-3600.0)/1000.0))
    #y2=(95.0-19)*np.exp(-((x-6500.0)/250.0)**2.0/2.0)
    #y3=0.5*np.sin((x-3600.0)/100.0)
    #y4=23*np.exp(-((x-6000.0)/200.0)**2.0/2.0)
    #y5=8*np.exp(-((x-5250.0)/200.0)**2.0/2.0)
    #y6=-20*np.exp(-((x-9400.0)/150.0)**2.0/2.0)
    #y7=0.0#ran.randn(len(x))*1.0
    #y8=8.0/(1.+np.exp((x-6100)/100))
    #y=(y1+y2+y3+y4+y5+y6+y7+y8)*fac
    
    f=open(dir_tem+'/response_c.txt','r')
    x1=[]
    y1=[]
    for line in f:
        data=line.replace('\n','').split(',')
        data=filter(None,data)
        x1.extend([np.float(data[0])])
        y1.extend([np.float(data[1])])
    y1=np.array(y1)
    x1=np.array(x1)
    y=interp1d(x1,y1,bounds_error=False,fill_value=0.)(x)
    y=y*fac
    return y            

def extintion_c(wave,dir_tem='libs'):
    f=open(dir_tem+'/extintion_curve.txt','r')
    x1=[]
    y1=[]
    for line in f:
        data=line.replace('\n','').split(' ')
        data=filter(None,data)
        x1.extend([np.float(data[0])])
        y1.extend([np.float(data[1])])
    y1=np.array(y1)
    x1=np.array(x1)
    Kl=interp1d(x1,y1,bounds_error=False,fill_value=0.)(wave)
    return Kl
=======
        
        
                
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8

def cube_conv_t(outf,dir_o='',psfi=0,nl=7,fov=30.0,thet=0.0,ifutype="MaNGA"):
    if "MaNGA" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        nl=11
        sigma_inst=25.0
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
    elif "MUSE" in ifutype:
        pix_s=0.2#0.1#0.025#arcsec
        scp_s=300.0#150.0#300.0#1200.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        nl=int(fov*scp_s/fibA/2)+1
        sigma_inst=25.0
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
    else:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi           
    scalep=1.0/scp_s
    dap=scalep
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    ns=3*nl*(nl-1)+1
    Dfib=fibA*scalep
    Rifu=Dfib*((2.0*nl-1.0)/2.0-0.5)
    xfib0=-Rifu
    yfib0=0
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    xifu=np.zeros(ns)
    yifu=np.zeros(ns)
    ini=0
    for i in range(0, nl):
        nt=nl*2-1-i
        yfib=yfib0+i*dyf*Dfib
        for j in range(0, nt):
            xfib=xfib0+(j*dxf+0.5*i)*Dfib
            xifu[ini]=xfib
            yifu[ini]=yfib
            ini=ini+1
        if i > 0:
            for j in range(0, nt):
                xfib=xfib0+(j*dxf+0.5*i)*Dfib
                xifu[ini]=xfib
                yifu[ini]=-yfib
                ini=ini+1
    ndt=35
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.0]
    dit[1,:]=[+0.00,+1.0]
    dit[2,:]=[+0.00,-1.0]
    dit[3,:]=[+0.25,+0.5]
    dit[4,:]=[+0.25,-0.5]
    dit[5,:]=[-0.25,+0.5]
    dit[6,:]=[-0.25,-0.5]
    dit[7,:]=[+0.50,+0.5]
    dit[8,:]=[+0.50,-0.5]
    dit[9,:]=[+0.50,+1.5]
    dit[10,:]=[+0.50,-1.5]
    dit[11,:]=[-0.50,+0.5]
    dit[12,:]=[-0.50,-0.5]
    dit[13,:]=[-0.50,+1.5]
    dit[14,:]=[-0.50,-1.5]
    dit[15,:]=[+1.00,+0.0]
    dit[16,:]=[+1.00,+0.5]
    dit[17,:]=[+1.00,-0.5]
    dit[18,:]=[+1.00,+1.0]
    dit[19,:]=[+1.00,-1.0]    
    dit[20,:]=[-1.00,+1.0]
    dit[21,:]=[-1.00,-1.0]
    dit[22,:]=[+1.50,+0.0]
    dit[23,:]=[+1.50,-0.5]
    dit[24,:]=[+1.50,+0.5]
    dit[25,:]=[+0.50,+0.0]
    dit[26,:]=[-0.50,+0.0]
    dit[27,:]=[+0.35,-0.9]
    dit[28,:]=[-0.35,-0.9]
    dit[29,:]=[+0.85,-1.3]
    dit[30,:]=[-0.85,-1.3]
    dit[31,:]=[+0.78,-0.25]
    dit[32,:]=[+0.78,+0.25]
    dit[33,:]=[-0.35,+0.73]
    dit[34,:]=[+0.50,+0.70]
    dyf=1.0
    dyt=Dfib/2.0/np.cos(30.0*np.pi/180.0)
    dxt=Dfib/2.0
    ndt=3
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.00]+ran.randn(2)*0.025
    dit[1,:]=[+0.00,+dyt/1.0]+ran.randn(2)*0.025
    dit[2,:]=[-dxt, +dyt/2.0]+ran.randn(2)*0.025
    
    spec_ifu=np.zeros(ndt*ns)
    x_ifu=np.zeros(ndt*ns)
    y_ifu=np.zeros(ndt*ns)
    facto=(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)
    con=0
    for i in range(0, ndt):
        for j in range(0, ns):
            xo=xifu[j]+dit[i,0]
            yo=yifu[j]+dyf*dit[i,1]    
            #val_t=1.0
            val_t=integ_flux(fibB*scalep,xo,yo,type=0,sig=2)
            spec_ifu[con]=val_t*facto#+ran.randn(1)*0.001
            x_ifu[con]=xo
            y_ifu[con]=yo
            con=con+1
    nl=int(round((np.amax([np.amax(x_ifu),-np.amin(x_ifu),np.amax(y_ifu),-np.amin(y_ifu)])+1)*2/pix_s))
    print nl
    ifu_f=np.zeros([3,nl,nl])
    ifu=np.zeros([nl,nl])
    ifu2=np.zeros([nl,nl])
    xo=-nl/2*pix_s
    yo=-nl/2*pix_s
    print xo,yo
    xi=xo
    xf=xo
    for i in range(0, nl):
        xi=xf
        xf=xf+pix_s
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+pix_s
            spt_new=0
            Wgt=0
            for k in range(0, len(x_ifu)):
                Rsp=np.sqrt((x_ifu[k]-(xf+xi)/2.0)**2.0+(y_ifu[k]-(yf+yi)/2.0)**2.0)
                if Rsp <= fibA*scalep*1.4/2.0:
#                    Wg=np.exp(-(Rsp/pix_s)**2.0/2.0)
                    Wg=np.exp(-(Rsp/0.6)**2.0/2.0)
                    spt_new=spec_ifu[k]*Wg+spt_new
                    Wgt=Wgt+Wg
#                    print Wgt,spt_new/Wgt,Wg,xi,yi,Rsp,k
            ifu2[j,i]=integ_flux(pix_s,(xf+xi)/2.0,(yf+yi)/2.0,type=1,sig=2)
            if Wgt == 0:
                Wgt=1
            ifu[j,i]=spt_new/Wgt
    #sys.exit()
    ifu_m=np.ones([nl,nl])
    ifu_m[np.where(ifu == 0)]=0
    ifu_r=ifu/ifu2        
    valt0=np.nansum(ifu)
    valt1=np.nansum(ifu2*ifu_m)
    valt2=valt0/valt1
    ifu_f[0,:,:]=ifu
    ifu_f[1,:,:]=ifu2*ifu_m
    ifu_f[2,:,:]=ifu_r
    h1=pyf.PrimaryHDU(ifu_f)
    h=h1.header
    h["NAXIS"]=3 
    h["NAXIS3"]=3
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mock "+ifutype+" IFU"
    h["CRVAL1"]=0
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['FOV']=Rifu*2.0
    h['VAL0']=valt0
    h['VAL1']=valt1
    h['VAL2']=valt2
    h['IFUCON']=(str(np.int(ns))+' ','NFibers')
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
<<<<<<< HEAD

def fib_arc_raw(dirtemp="./",temp_1='arcs.txt'):
    cdelt_w=0.5
    crval_w=3522.0
    crpix_w=1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    arcB=ssp_extract_arc(wave_f_b,dir_tem=dirtemp,col='b')
    flatB=flat_b(wave_f_b)
    spec_ifu_b=arcB#*flatB#*100.0
    arcR=ssp_extract_arc(wave_f_r,dir_tem=dirtemp,col='r')
    flatR=flat_r(wave_f_r)
    spec_ifu_r=arcR#*flatR#*100.0
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r) 
    #nt1=np.where(np.isfinite(spec_b) ==  True)
    #nt2=np.where(np.isfinite(spec_r) ==  True)
    #spec_b=spec_b[nt1]
    #spec_r=spec_r[nt2]
    return spec_b,spec_r

    
def fib_flat_raw(dirtemp="./"):
    cdelt_w=0.5
    crval_w=3522.0
    crpix_w=1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    flatB=flat_b(wave_f_b)
    spec_ifu_b=flatB*20000.0*1.1#*0.9
    flatR=flat_r(wave_f_r)
    spec_ifu_r=flatR*20000.0*1.1#*0.9
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r)
    #nt1=np.where(np.isfinite(spec_b) ==  True)
    #nt2=np.where(np.isfinite(spec_r) ==  True)
    #spec_b=spec_b*wave_b/dwave_b/2.0/1633.0
    #spec_r=spec_r*wave_r/dwave_r/2.0/2216.0
    #spec_b=spec_b[nt1]
    #spec_r=spec_r[nt2]
    return spec_b,spec_r    
    
def fib_sky_raw(outf,xx,yy,fibid,dirtemp="./",template4="sky.txt",dir_o='',Flux_m=20.0,sp_res=0,SNi=15.0,pdf=2,ifutype="SDSS",expt=900.0):
    vel_light=299792.458
    if "SDSS" in ifutype:
=======
    
def fib_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="SDSS"):
    nh=dens
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol
    sup=4.0*np.pi*Rs**2.0
    vel_light=299792.458
    no_nan=0
    if "SDSS" in ifutype:
        pix_s=0.5#arcsec
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        scp_s=60.4#microns per arcsec
        fibA=210.0
        fibB=180.0
        sigma_inst=25.0
<<<<<<< HEAD
        if sp_res <=0:
            sp_res=2000.0
    elif "BOSS" in ifutype:
        scp_s=60.4#microns per arcsec
        fibA=190.0
        fibB=120.0
        sigma_inst=25.0
        if sp_res <=0:
            sp_res=2000.0
    elif "APO" in ifutype:
        scp_s=60.482#microns per arcsec
        fibA=190.0
        fibB=120.0
        sigma_inst=25.0
        if sp_res <=0:
            sp_res=2000.0
    elif "LCO" in ifutype:
        scp_s=91.475#microns per arcsec
        fibA=190.0
        fibB=180.0
        fibB=120.0
        sigma_inst=25.0
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
=======
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        sigma_inst=25.0
<<<<<<< HEAD
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        scp_s=300.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
        if sp_res <=0:
            sp_res=4000.0
    else:
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
        if sp_res <=0:
            sp_res=2000.0           
    scalep=1.0/scp_s
    
    [sky_template,wave_s,crval_s,cdelt_s,crpix_s]=ssp_extract_sky(template4)
    mag_1=band_mag(sky_template*1e-18,crpix_s,cdelt_s,crval_s,dir='legacy/',k=3)
    dmag=-2.5*np.log10(10.0**(-0.4*(Flux_m))*np.pi*(fibB*scalep/2.0)**2.0)-mag_1
    dflux=10.0**(-0.4*dmag)*1e-18
    sky_template=sky_template*dflux/1e-16
    cdelt_w=cdelt_s#1.25
    crval_w=3522.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    #if Flux_m != 20.0:
    #    Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    #else:
    #    Flux_m=20.0
    #s_nr=noise_sig(wave_f,SNi-1.0)
    nw=len(wave_f)
    #t_noise=1.0/Flux_m
    noise=0.0#t_noise*ran.randn(nw)/s_nr*0.0
    fac=1.0*expt/900.0
    spect_sky=interp1d(wave_s,sky_template,bounds_error=False,fill_value=0.)(wave_f)
    mag_f=band_mag(spect_sky*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    spec_ifu=noise+spect_sky
    flux_d=flux_distort_corrvec(wave_f,xx,yy,fibid)
    cal_v=response_v(wave_f,fac=fac)
    spec_ifuc=spec_ifu*cal_v*flux_d/1e-1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    spec_ifu_b=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_b)
    spec_ifu_r=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_r)    
    flatB=flat_b(wave_f_b)
    spec_ifu_b=spec_ifu_b*flatB
    flatR=flat_r(wave_f_r)
    spec_ifu_r=spec_ifu_r*flatR
    
    
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r) 
    #nt1=np.where(np.isfinite(spec_b) ==  True)
    #nt2=np.where(np.isfinite(spec_r) ==  True)
    spec_b=spec_b/(wave_b/dwave_b/2.0/1633.0)
    spec_r=spec_r/(wave_r/dwave_r/2.0/2216.0)
    #spec_b=spec_b[nt1]
    #spec_r=spec_r[nt2]
    
    if pdf ==1 :
        max_val1=np.amax(spec_ifu)*1.25
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        plt.title('mf='+str(np.round(mag_f,3)))
        plt.plot(wave_f,spec_ifu)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'.pdf')
        plt.close()
    
        max_val1=np.amax(spec_ifuc)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.plot(wave_f,spec_ifuc)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'C.pdf')
        plt.close()
    
        max_val1=np.amax(spec_b)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.plot(pix_b,spec_b)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'B.pdf')
        plt.close()
        
        max_val1=np.amax(spec_r)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.plot(pix_r,spec_r)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'R.pdf')
        plt.close()

    return spec_b,spec_r

def plots_hol_type(mjd,name='plot_hole',ra_0=180.0,dec_0=0.0,phi=32.78,scp_s=60.4,fibB=120.0):
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=1,name=name+'_A',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=0,name=name+'_B',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=2,name=name+'_C',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=1,name=name+'_A',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=0,name=name+'_B',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=2,name=name+'_C',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=1,name=name+'_A',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=0,name=name+'_B',lshf=4000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=2,name=name+'_C',lshf=4000.0,scp_s=scp_s,fibB=fibB)

    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=1,name=name+'_A',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=0,name=name+'_B',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=2,name=name+'_C',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=1,name=name+'_A',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=0,name=name+'_B',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=2,name=name+'_C',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=1,name=name+'_A',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=0,name=name+'_B',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=2,name=name+'_C',lshf=4500.0,scp_s=scp_s,fibB=fibB)
    
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=1,name=name+'_A',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=0,name=name+'_B',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=2,name=name+'_C',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=1,name=name+'_A',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=0,name=name+'_B',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=2,name=name+'_C',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=1,name=name+'_A',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=0,name=name+'_B',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=2,name=name+'_C',lshf=5000.0,scp_s=scp_s,fibB=fibB)
    
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=1,name=name+'_A',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=0,name=name+'_B',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-10.0,phi=phi,shf=2,name=name+'_C',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=1,name=name+'_A',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=0,name=name+'_B',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0-00.0,phi=phi,shf=2,name=name+'_C',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=1,name=name+'_A',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=0,name=name+'_B',lshf=5500.0,scp_s=scp_s,fibB=fibB)
    plot_fiberhole(mjd,ra_0=ra_0,dec_0=dec_0+10.0,phi=phi,shf=2,name=name+'_C',lshf=5500.0,scp_s=scp_s,fibB=fibB)

def plot_fiberhole(mjd,name='plot_hole',ra_0=180.0,hat=0.0,expof=0.0,dec_0=0.0,phi=32.78,dr=3.0,dn=3.0,shf=1,lshf=4000.0,scp_s=60.4,fibA=190.0,fibB=120.0):
    expot=dr*2*3600.0
    expt=expot/(dn*2)
    cdelt_w=1.25
    crval_w=3500.0
    crpix_w=1
    wave_f=np.arange(crval_w,6000.0,cdelt_w)
    n_d=5
    col=['blue','deepskyblue','green','lime','gold','red']
    wave_d=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    wave_d=wave_d+(wave_d[1]-wave_d[0])/2.0
    wave_d=wave_d[0:n_d]
    wave_l=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    #print wave_l
    dnt=np.int(dn*2+1)
    xoR=np.zeros([n_d+1,dnt])
    yoR=np.zeros([n_d+1,dnt])
    for i in range(0, dnt):
        expoff=expof+(i-dn)*expt+2000.0
        if i == 0:
            toffs=expoff
        ha=lst_c(tai=((np.float(mjd)+.25)*3600.0*24.0+expoff+expt/2.0))-ra_0/180.0*12.0 
        if shf == 0:
            ha1=ha
        elif shf ==1:
            #ha1=lst_c(tai=((np.float(mjd)+0.25)*3600.0*24.0+toffs+expot/2.0))-ra_0/180.0*12.0     
            ha1=ha
        else:
            ha1=0.0
        R=refrac_dif(wave_l,ha,dec=dec_0,phi=phi,lo=5070.0,T=7.0,P=600.0,f=8.0)
        Rs=refrac_dif(lshf,ha1,dec=dec_0,phi=phi,lo=5070.0,T=7.0,P=600.0,f=8.0)
        pa=paralactic_angle(ha,dec=dec_0,phi=phi)
        #print ha,ha1,expoff#,R[0],Rs,pa,expoff,ha1
        R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
        xo=0#ran.randn(1)*0.025
        yo=0#ran.randn(1)*0.025
        for w in range(0, n_d+1):
            R_adr=np.dot(R2,np.array([0,R[w]-Rs]))
            yoR[w,i]=yo+R_adr[1]#+R_t[1]
            xoR[w,i]=xo+R_adr[0]#+R_t[0]
    dn1=600#300
    expt=expot/(dn1*2)
    expt1=expot/(dn*2)
    dnt=np.int(dn1*2+1)
    xoR1=np.zeros([n_d+1,dnt])
    yoR1=np.zeros([n_d+1,dnt])
    for i in range(0, dnt):
        expoff=expof+(i-dn1)*expt+3600.0
        expoff1=expof+(np.round((i-dn1*0)*expt*1.0/(expt1*1.0))-dn)*expt1+2000.0#2000.0
        if i == 0:
            toffs=expoff
        ha=lst_c(tai=((np.float(mjd)+.25)*3600.0*24.0+expoff+expt/2.0))-ra_0/180.0*12.0 
        if shf == 0:
            ha1=ha
        elif shf == 1:
            #ha1=hat
            #ha1=lst_c(tai=((np.float(mjd)+0.25)*3600.0*24.0+toffs+expot/2.0))-ra_0/180.0*12.0
            ha1=lst_c(tai=((np.float(mjd)+.25)*3600.0*24.0+expoff1+expt1/2.0))-ra_0/180.0*12.0 
        else:
            ha1=0.0     
        R=refrac_dif(wave_l,ha,dec=dec_0,phi=phi,lo=5070.0,T=7.0,P=600.0,f=8.0)
        Rs=refrac_dif(lshf,ha1,dec=dec_0,phi=phi,lo=5070.0,T=7.0,P=600.0,f=8.0)
        pa=paralactic_angle(ha,dec=dec_0,phi=phi)
        #print ha,ha1,expoff,expoff1,(np.round((i-dn1*0)*expt*1.0/(expt1*1.0))-dn)#,R[0],Rs,pa,expoff
        R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
        xo=0#ran.randn(1)*0.025
        yo=0#ran.randn(1)*0.025
        for w in range(0, n_d+1):
            R_adr=np.dot(R2,np.array([0,R[w]-Rs]))
            yoR1[w,i]=yo+R_adr[1]#+R_t[1]
            xoR1[w,i]=xo+R_adr[0]#+R_t[0]
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    #import matplotlib.style
    #matplotlib.style.use('classic')
    fig, ax = plt.subplots(figsize=(6.0,5.65))
    plt.xlim(-fibB/scp_s/2.0*1.05,fibB/scp_s/2.0*1.05)
    plt.ylim(-fibB/scp_s/2.0*1.05,fibB/scp_s/2.0*1.05)
    for w in range(0, n_d+1):
        plt.plot(xoR[w,:],yoR[w,:],'o',color=col[w],lw=0.0,markeredgewidth=0)
        #plt.plot(xoR1[w,:],yoR1[w,:],'.',color=col[w])#lw=0.5)
        
        plt.scatter(xoR1[w,:],yoR1[w,:], c=col[w], alpha=0.8,s=0.1,edgecolors=col[w])#'none')
        
    plt.text(-fibB/scp_s/2.0*0.98,fibB/scp_s/2.0*0.95,r'$\delta='+str(np.round(dec_0,2))+'$',fontsize=15)
    plt.text(-fibB/scp_s/2.0*0.98,fibB/scp_s/2.0*0.82,r'$\lambda_c='+str(np.round(lshf,2))+'$',fontsize=15)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.98,r'$6000\AA$',color=col[5],fontsize=17)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.88,r'$5500\AA$',color=col[4],fontsize=17)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.78,r'$5000\AA$',color=col[3],fontsize=17)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.68,r'$4500\AA$',color=col[2],fontsize=17)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.58,r'$4000\AA$',color=col[1],fontsize=17)
    plt.text(-fibB/scp_s/2.0*0.15,-fibB/scp_s/2.0*0.48,r'$3500\AA$',color=col[0],fontsize=17)
    nr=100
    theta=np.arange(0,nr)/(nr-1.0)*2.0*np.pi
    x=fibB/scp_s/2.0*np.cos(theta)
    y=fibB/scp_s/2.0*np.sin(theta)
    plt.xlabel(r'$x(\lambda,\lambda_c,h,\delta)\ (arcsec)$',fontsize=14)
    plt.ylabel(r'$y(\lambda,\lambda_c,h,\delta)\ (arcsec)$',fontsize=14)
    plt.plot(x,y)#,fillstyle='none')
    fig.tight_layout()
    plt.savefig(name+'_'+str(np.round(lshf,0)).replace('.0','')+'_'+str(np.round(dec_0,2)).replace('.','')+'.jpeg',dpi = 1000)
    plt.close()
    
def conv(xt,ke=2.5):
    from scipy.ndimage.filters import gaussian_filter1d as filt1d
    nsf=len(xt)
    krn=ke
    xf=filt1d(xt,ke)
    return xf
    
def fib_star_raw_plot(outf,Av_g,mag_s,fib_id=1,xx=0,yy=0,dsep=1.0,lshf=5000,star_t='F0II (25291)',alp=0.0,bet=0.0,ha=0.0,ha1=0,dec_0=0.0,dirtemp="libs/",template4="libs/sky_model.txt",template3="libs/spEigenStar-55734.fits",template6="libs/da012500_800.dat",wd=False,dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,vel_0=100.00,pdf=1,ifutype="APO",expt=900.0,outs=0,apot=0):
    vel_light=299792.458
    if "APO" in ifutype:
        scp_s=60.482#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "LCO" in ifutype:
        scp_s=91.475#microns per arcsec
        fibA=190.0
        fibB=180.0
        fibB=120.0
        obs_l=-29.0146
        sigma_inst=25.0
        beta_s=4.6#76
=======
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        pix_s=0.2#0.1#0.025#arcsec
        scp_s=300.0#150.0#300.0#1200.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
<<<<<<< HEAD
            sp_res=2000.0
    else:
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0  
        obs_l=32.78   
        sigma_inst=25.0
        beta_s=4.76
=======
            sp_res=4000.0
    else:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
<<<<<<< HEAD
            sp_res=2000.0
    #seeing=psf_pdf(1,psf_m=seeing)[0]           
    scalep=1.0/scp_s
    #print scalep,scp_s
    dlam=(1.0+(vel_0/vel_light))
    Dfib=fibB*scalep
    n_pt=100000.0
    #n_pt=100000.0
    #n_pt=10000000.0
    
    [sky_template,wave_s,crval_s,cdelt_s,crpix_s]=ssp_extract_sky(template4)
    mag_1=band_mag(sky_template*1e-18,crpix_s,cdelt_s,crval_s,dir='legacy/',k=3)
    dmag=-2.5*np.log10(10.0**(-0.4*(Flux_m))*np.pi*(fibB*scalep/2.0)**2.0)-mag_1
    mag_sk=mag_1+dmag
    dflux=10.0**(-0.4*dmag)*1e-18
    sky_template=sky_template*dflux/1e-16
    
    #[ssp_template,wave,crval,cdelt,crpix]=ssp_extract_star(template3,star_type=star_t)
    if wd == False:
        [ssp_template,wave,crval,cdelt,crpix]=ssp_extract_star(template3,star_type=star_t)
    else:
        [ssp_template,wave,crval,cdelt,crpix]=ssp_extract_wd(template6)
        star_t='DA White D.'
    mag_2=band_mag(ssp_template,crpix,cdelt,crval,dir='legacy/',k=3)
    dmag=mag_s-mag_2
    #print mag_2,dmag
    
    
    dflux=10.0**(-0.4*dmag)
    ssp_template=ssp_template*dflux/n_pt
    dust_rat_ssp=A_l(3.1,wave)
    
    
    cdelt_w=cdelt
    crval_w=3522.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    nw=len(wave_f)
    nw_s=len(wave)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    spect_t=np.zeros(nw_s)
    spect_ot=np.zeros(nw_s)
    spect=np.zeros(nw)
    #noise=0.0#t_noise*ran.randn(nw)/s_nr*0.0
    #Ft=0
    fac=1.0*expt/900.0
    
    n_d=100
    wave_d=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    wave_d=wave_d+(wave_d[1]-wave_d[0])/2.0
    wave_d=wave_d[0:n_d]
    wave_l=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    R=refrac_dif(wave_l,ha,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)#,vapor=True)
    Rs=refrac_dif(lshf,ha1,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)#,vapor=True)
    pa=paralactic_angle(ha,dec=dec_0,phi=obs_l)
    R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
    R1=np.array([[np.cos((alp+bet)*np.pi/180.0),-np.sin((alp+bet)*np.pi/180.0)],[np.sin((alp+bet)*np.pi/180.0),np.cos((alp+bet)*np.pi/180.0)]])
    Rt=0.0
    if apot == 1:
        Rt=Rt+dsep*1e3*scalep
    R_t=np.dot(R1,np.array([0,Rt]))
    print ha,ha1,Rs,pa#,lshf
    phie=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s) 
    thee=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s)
    
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    flux_lost=np.zeros(n_d+1)
    for w in range(0, n_d+1):
        R_adr=np.dot(R2,np.array([0,R[w]-Rs]))
        yoR=yo+R_adr[1]+R_t[1]
        xoR=xo+R_adr[0]+R_t[0]
        r=np.sqrt((xoR-phie)**2.0+(yoR-thee)**2.0)
        nt=np.where(r <= fibB*scalep/2.0)[0]
        flux_lost[w]=np.float(len(nt))
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    print len(r),len(nt),fibB*scalep/2.0,fibB,1/scalep
    print np.float(len(r))/np.float(len(nt))," Fiber-PSF factor",ifutype,beta_s,psfi
    flux_lost=flux_lost/np.float(len(nt))
    Av=Av_g
    dust=10**(-0.4*Av_g*dust_rat_ssp)
    if len(nt) > 0:
        for k in range(0, len(nt)):
            spect_s=ssp_template/1e-16*dust#/1e-16
            spect_o=ssp_template*dust
            spect_ot=spect_o+spect_ot 
            spect_t=spect_s+spect_t
    spect_t[np.isnan(spect_t)]=0
    spect_ot[np.isnan(spect_ot)]=0        
    spect_tf=shifts(spect_t,wave,dlam)
    spect_otf=shifts(spect_ot,wave,dlam)        
    spect_i=inst_disp(wave,spect_tf,sigma_inst)
    spect_oi=inst_disp(wave,spect_otf,sigma_inst)        
    spect_ii=spec_resol(wave,spect_i,sp_res)
    spect_oii=spec_resol(wave,spect_oi,sp_res)
    
    flux_lost_f=interp1d(wave_l,flux_lost,bounds_error=False,fill_value=0.)(wave_f)
    flux_lost_f=conv(flux_lost_f,ke=60)
    Kvl=extintion_c(wave_f)
    Xa=airmas(ha, dec=dec_0,phi=obs_l)
    print Xa
    at_ext=10.0**(-0.4*Xa*Kvl)
            
    spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
    specto=interp1d(wave,spect_oii,bounds_error=False,fill_value=0.)(wave_f)
    spect[np.isnan(spect)]=0
    specto[np.isnan(specto)]=0
    specto=specto/1e-16
    spect_sky=interp1d(wave_s,sky_template,bounds_error=False,fill_value=0.)(wave_f)
    spec_ifu=spect
    
    mag_2=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    spec_ifu=spec_ifu*at_ext*flux_lost_f+spect_sky
    mag_f=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    #print 10**(-0.4*(mag_s-mag_2)),mag_s,mag_2,seeing
    flux_d=flux_distort_corrvec(wave_f,xx,yy,fib_id)
    cal_v=response_v(wave_f,fac=fac)
    spec_ifuc=spec_ifu*cal_v*flux_d/1e-1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    spec_ifu_b=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_b)
    spec_ifu_r=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_r)    
    flatB=flat_b(wave_f_b)
    spec_ifu_b=spec_ifu_b*flatB
    flatR=flat_r(wave_f_r)
    spec_ifu_r=spec_ifu_r*flatR
    cal_v_b=interp1d(wave_f,cal_v*flux_d/1e-1,bounds_error=False,fill_value=0.)(wave_f_b)
    cal_v_r=interp1d(wave_f,cal_v*flux_d/1e-1,bounds_error=False,fill_value=0.)(wave_f_r)    
    
    #cal_v_b1=cal_v_b/(wave_f_b/cdelt_w/2.0/1633.0)
    #cal_v_r1=cal_v_r/(wave_f_r/cdelt_w/2.0/2216.0)
    
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r)
    cal_b=interp1d(wave_f_b,cal_v_b,bounds_error=False,fill_value=0.0)(wave_b)
    cal_r=interp1d(wave_f_r,cal_v_r,bounds_error=False,fill_value=0.0)(wave_r)
    spec_b=spec_b/(wave_b/dwave_b/2.0/1633.0)
    spec_r=spec_r/(wave_r/dwave_r/2.0/2216.0)
    
    flatB_b=flat_b(wave_b)
    flatR_r=flat_r(wave_r)
    cal_b=cal_b/(wave_b/dwave_b/2.0/1633.0)
    cal_r=cal_r/(wave_r/dwave_r/2.0/2216.0)
    
    
    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    ifu_o=np.zeros([nw,1])
    
    
    if pdf ==1 :
        max_val1=np.amax(spec_ifu)*1.25
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(star_t+' mf='+str(np.round(mag_f,3))+' mr='+str(np.round(mag_2,3))+' mrs='+str(np.round(mag_s,3))+' ms='+str(np.round(mag_sk,3)))
        plt.plot(wave_f,spec_ifu)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'.jpg')
        plt.close()
    
        max_val1=np.amax(specto)*1.1
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[10^{-16} erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_f,spect_sky,alpha=0.4,color='grey')
        plt.plot(wave_f,spec_ifu,alpha=1.0,color='black')#spectro final
        plt.plot(wave_f,specto,color='blue',alpha=0.8)#spectro original
        plt.plot(wave_f,spect*at_ext*flux_lost_f,alpha=0.9,color='red')
  #      plt.plot(wave_f,at_ext)
  #      plt.plot(wave_f,flux_lost_f)
  #      plt.plot(wave_f,flux_lost_f*at_ext)#*flux_d)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_O.jpg',dpi=1000)
        plt.close()

        max_val1=np.amax(spec_ifuc)*0.07
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_b,spec_b,alpha=1.0,color='blue',lw=0.8)
        plt.plot(wave_r,spec_r,alpha=1.0,color='red',lw=0.80)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_BR.jpg',dpi=1000)
        plt.close()

        max_val1=1000#np.amax(spec_ifuc)*0.07
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("$Counts/10^{-16} erg s^{-1}cm^{-2}\AA^{-1}$",fontsize=14)
        plt.title("Response")
        #plt.plot(wave_r,0*cal_r,alpha=1.0,linestyle='-'  ,color='black',label='Raw')
        plt.plot(wave_r,0*cal_r,alpha=1.0,linestyle='-.',color='black',label='Flatted')
        plt.plot(wave_r,cal_r,alpha=1.0,color='red', linestyle='--',lw=1.7)
        plt.plot(wave_b,cal_b,alpha=1.0,color='blue',linestyle='--',lw=1.7)
        plt.plot(wave_r,cal_r*flatR_r,alpha=1.0,color='red',lw=1.7)
        plt.plot(wave_b,cal_b*flatB_b,alpha=1.0,color='blue',lw=1.7)    
    #    plt.plot(wave_b,spec_b/flatB_b/cal_b,alpha=1.0,color='blue',lw=0.8)
    #    plt.plot(wave_r,spec_r/flatR_r/cal_r,alpha=1.0,color='red',lw=0.8)
        plt.legend(loc=1,frameon=False)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_Res.jpg',dpi=1000)
        plt.close()


        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,0.04)#max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("$10^{-16} erg s^{-1}cm^{-2}\AA^{-1}/Counts$",fontsize=14)
        plt.title("Response")        
        plt.plot(wave_r,0/cal_r,alpha=1.0,linestyle='-'  ,color='black',label='Raw')
        plt.plot(wave_r,0/cal_r,alpha=1.0,linestyle='-.',color='black',label='Flatted')
        plt.plot(wave_r,1/cal_r,alpha=1.0,color='red')
        plt.plot(wave_b,1/cal_b,alpha=1.0,color='blue')
        plt.plot(wave_b,1/flatB_b/cal_b,linestyle='-.',alpha=1.0,color='blue',lw=1.0)
        plt.plot(wave_r,1/flatR_r/cal_r,linestyle='-.',alpha=1.0,color='red',lw=1.0)
        plt.legend(loc=1,frameon=False) 
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_Res_1.jpg',dpi=1000)
        plt.close()
    
        max_val1=np.amax(spec_ifuc)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_f,spec_ifuc)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'C.jpg')
        plt.close()
    
        max_val1=np.amax(spec_b)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t+' Blue')
        plt.plot(pix_b,spec_b)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'B.jpg')
        plt.close()
        
        max_val1=np.amax(spec_r)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t+' Red')
        plt.plot(pix_r,spec_r)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'R.jpg')
        plt.close()
        
        
    if outs == 1:    
        ifu[:,0]=spec_ifu
        ifu_o[:,0]=specto
        h1=pyf.PrimaryHDU(ifu)
        h2=pyf.ImageHDU(ifu_e)
        h3=pyf.ImageHDU(ifu_1)
        h4=pyf.ImageHDU(ifu_m)
        h5=pyf.ImageHDU(ifu_o)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=nw
        h["NAXIS2"]=1
        h["COMMENT"]="Mock "+ifutype+" single spectra"
        h['CDELT1']=cdelt_w
        h['CRPIX1']=crpix_w
        h['CRVAL1']=crval_w
        h['CUNIT1']='Wavelength [A]'
        h['PSF']=seeing
        h['DFIB']=Dfib
        h['TYPE']=star_t
        h['REDSHIFT']=float(vel_0/vel_light)
        h['UNITS']='1E-16 erg/s/cm^2/A'
        hlist=pyf.HDUList([h1,h2,h3,h4,h5])
        hlist.update_extend()
        out_fit=dir_o+outf+'.fits'
        wfits_ext(out_fit,hlist)
        dir_o1=dir_o.replace(" ","\ ")
        out_fit1=dir_o1+outf+'.fits'
        sycall('gzip  '+out_fit1)
        
    
def boss_apogee_cont(lshf=5000,ha=0.0,dspmin=0.0,dspmax=5.0,dec_0=32.78,dir_o='./',psfi=[0],ifutype="APO",adrt=False):
    if "APO" in ifutype:
        scp_s=60.482#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        beta_s=2.0
        if psfi[0] <= 0:
            seeing=1.43
        else:
            seeing=psfi
    if "LCO" in ifutype:
        scp_s=91.475#microns per arcsec
        fibA=190.0
        fibB=180.0
        fibB=120.0
        obs_l=-29.0146
        beta_s=4.76
        if psfi[0] <= 0:
            seeing=0.6
        else:
            seeing=psfi
    scalep=1.0/scp_s
    #if adrt == True:
    n_pt=10000000.0
    #else:
    #n_pt=100000000.0
    ndsp=100
    dsep=np.arange(ndsp)/np.float(ndsp-1)*(dspmax-dspmin)+dspmin 
    n_psf=len(seeing)
    if adrt == True:
        magv1=[]
        magv2=[]
        magv3=[]
        magv4=[]
        magv5=[]
    else:
        magvf=[]
    for kt in range(0, n_psf):
        if adrt == True:
            mag1=np.zeros(ndsp)
            mag2=np.zeros(ndsp)
            mag3=np.zeros(ndsp)
            mag4=np.zeros(ndsp)
            mag5=np.zeros(ndsp)
            cdelt_w=1.2
            crval_w=1522.0
            crpix_w=1
            wave_f=np.arange(crval_w,11352.0,cdelt_w)
            n_d=100
            wave_d=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
            wave_d=wave_d+(wave_d[1]-wave_d[0])/2.0
            wave_d=wave_d[0:n_d]
            wave_l=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
            R=refrac_dif(wave_l,ha,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
            Rs=refrac_dif(lshf,ha,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
            phie=psf_moffat(np.int(n_pt),psf=seeing[kt],beta=beta_s) #ran.randn(np.int(n_pt))*seeing/2.0
            thee=psf_moffat(np.int(n_pt),psf=seeing[kt],beta=beta_s) #ran.randn(np.int(n_pt))*seeing/2.0    
            xo=ran.randn(1)*0.025
            yo=ran.randn(1)*0.025
            r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
            nt0=np.where(r <= fibB*scalep/2.0)[0]
        else:
            n_d=0
            magf=np.zeros(ndsp)
            #phie=ran.randn(np.int(n_pt))*seeing/2.0
            #thee=ran.randn(np.int(n_pt))*seeing/2.0
        for i in range(0, len(dsep)):
            Rt=dsep[i]*1e3*scalep
            if adrt == True:        
                flux_lost=np.zeros(n_d+1)
                for w in range(0, n_d+1):
                    yoR=yo-(R[w]-Rs)+Rt
                    xoR=xo
                    r=np.sqrt((xoR-phie)**2.0+(yoR-thee)**2.0)
                    nt=np.where(r <= fibB*scalep/2.0)[0]
                    flux_lost[w]=np.float(len(nt))
                #flux_lost=flux_lost/np.float(n_pt)*100.0
                flux_lost=flux_lost/np.float(len(nt0))*100.0
                flux_lost_f=interp1d(wave_l,flux_lost,bounds_error=False,fill_value=0.)(wave_f)
                #print band_val(flux_lost_f,wave_f,dir='legacy/',k=1)
                mag1[i]=band_val(flux_lost_f,wave_f,dir='legacy/',k=1)
                mag2[i]=band_val(flux_lost_f,wave_f,dir='legacy/',k=2)
                mag3[i]=band_val(flux_lost_f,wave_f,dir='legacy/',k=3)
                mag4[i]=band_val(flux_lost_f,wave_f,dir='legacy/',k=4)
                mag5[i]=band_val(flux_lost_f,wave_f,dir='legacy/',k=5)
                print Rt,dsep[i],mag3[i]
            else:
                magf[i]=psf_moffat_flux(psf=seeing[kt],beta=beta_s,rfib=fibB*scalep/2.0,dt=Rt)#*100.0
                #magf[i]=magf[i]/magf[0]
                #magf[i]=flux_lost[0]
                print Rt,dsep[i],magf[i]/magf[0]*100.0
                #magf[i]=magt
            #sys.exit()
        #print dsep
        if adrt == False:
            magf=magf*100.0/magf[0]
            magvf.extend([magf])
        else:
            magv1.extend([mag1])
            magv2.extend([mag2])
            magv3.extend([mag3])
            magv4.extend([mag4])
            magv5.extend([mag5])
    matplotlib.use('agg')
    color_v=['blue','deepskyblue','green','gold','red']
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,5.5))
    ax.set_xlim(0,dspmax*1e3)
    ax.set_ylim(1e-5,1e2)
    ax.set_xlabel("$Fiber\ Separation\ [\mu m]$",fontsize=14)
    #ax.set_ylabel("Flux % $[F_{inFib}/F_{psf}]$",fontsize=14)
    ax.set_ylabel("Flux % $[F/F_{0}]$",fontsize=14)
    if adrt == True:
        plt.title('ha='+str(np.round(ha,3))+' dec='+str(np.round(dec_0,3))+' PSF_fwhm='+str(np.round(seeing[0],3))+' arcsec S.plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        for i in range(0, n_psf):
            plt.semilogy(dsep*1e3,magv1[i],color='blue',label=r'$\%\ in\ u\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3,magv2[i],color='deepskyblue',label=r'$\%\ in\ g\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3,magv3[i],color='green',label=r'$\%\ in\ r\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3,magv4[i],color='gold',label=r'$\%\ in\ i\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3,magv5[i],color='red',label=r'$\%\ in\ z\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
    else:
        if n_psf == 1:
            plt.title('PSF_fwhm='+str(np.round(seeing[0],3))+' arcsec, scale plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        else:
            plt.title('scale plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        for i in range(0, n_psf):
            plt.semilogy(dsep*1e3,magvf[i],color=color_v[i % 5],label=r'$PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
    plt.semilogy([fibB/1.0,fibB/1.0],[1e-5,1e3],color='black',linestyle='-',label='$Fiber\ D.size$')
    plt.semilogy([2.6*1e3,2.6*1e3],[1e-5,1e3],color='black',linestyle='--',label='$Fiber\ separation$')
    plt.semilogy([4*1e3,4*1e3],[1e-5,1e3],color='black',linestyle='-.',label='$Fiber\ collision$')
    plt.legend(loc=1,fontsize=11)
    fig.tight_layout()
    plt.savefig(dir_o+'flux_ratio_microns_'+str(np.int(np.abs(ha)))+'_'+str(np.round(seeing[0],3)).replace('.','')+'_'+ifutype+'.jpg')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(6,5.5))
    ax.set_xlim(0,dspmax*1e3*scalep)
    ax.set_ylim(1e-5,1e2)
    ax.set_xlabel("$Fiber\ Separation\ [arcsec]$",fontsize=14)
    #ax.set_ylabel("Flux % $[F_{inFib}/F_{psf}]$",fontsize=14)
    ax.set_ylabel("Flux % $[F/F_{0}]$",fontsize=14)
    if adrt == True:
        plt.title('ha='+str(np.round(ha,3))+' dec='+str(np.round(dec_0,3))+' PSF_fwhm='+str(np.round(seeing[0],3))+' arcsec S.plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        for i in range(0, n_psf):
            plt.semilogy(dsep*1e3*scalep,magv1[i],color='blue',label=r'$\%\ in\ u\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3*scalep,magv2[i],color='deepskyblue',label=r'$\%\ in\ g\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3*scalep,magv3[i],color='green',label=r'$\%\ in\ r\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3*scalep,magv4[i],color='gold',label=r'$\%\ in\ i\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
            plt.semilogy(dsep*1e3*scalep,magv5[i],color='red',label=r'$\%\ in\ z\ band\ PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
    else:
        if n_psf == 1:
            plt.title('PSF_fwhm='+str(np.round(seeing[0],3))+', arcsec S.plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        else:
            plt.title('scale plate='+str(np.round(scp_s,3))+' microns/arcsec, '+ifutype,fontsize=11)
        for i in range(0, n_psf):
            plt.semilogy(dsep*1e3*scalep,magvf[i],color=color_v[i % 5],label=r'$PSF='+str(np.round(seeing[i],3))+'\ arcsec$')
    plt.semilogy([fibB/1.0*scalep,fibB/1.0*scalep],[1e-5,1e3],color='black',linestyle='-',label='$Fiber\ D.size$')
    plt.semilogy([2.6*1e3*scalep,2.6*1e3*scalep],[1e-5,1e3],color='black',linestyle='--',label='$Fiber\ separation$')
    plt.semilogy([4.0*1e3*scalep,4.0*1e3*scalep],[1e-5,1e3],color='black',linestyle='-.',label='$Fiber\ collision$')
    plt.legend(loc=1,fontsize=11)
    fig.tight_layout()
    plt.savefig(dir_o+'flux_ratio_arcsec_'+str(np.int(np.abs(ha)))+'_'+str(np.round(seeing[0],3)).replace('.','')+'_'+ifutype+'.jpg')
    plt.close()
    

def RM_spectra(wave,flux,bhmass,bhacre,mjd=56008.0,mjd0=56008.0,tp=75.0,dla=50.0,dl=1.0,z_0=0.0,alpha=0.5,Ewhb=0.75,Ewha=2.0,pw_a=0.68,beta=0.85):#beta=0.5
    dwave=np.zeros(len(wave))
    nt=np.argmax(bhmass)
    bh_mass=np.amax(bhmass)
    if np.nansum(flux) > 0:
        if bh_mass > 0:
            bh_acret=bhacre[nt]        
            for i in range(0, len(wave)):
                if i < len(wave)-1:
                    dwave[i]=wave[i+1]-wave[i]
                if i == len(wave)-1:
                    dwave[i]=dwave[i-10]
            L500=band_val(flux,wave,dir='legacy/',k=2)
            L500=L500*dla
            #print 'L500 [log10]='+str(np.log10(L500*1e44))
            #L600=band_val(flux,wave,dir='legacy/',k=3)
            logL5100=np.log10(L500)+4.0#4.2
            vel_light=299792.458
            G=6.67408e-11
            fblr=5.0
            logRblr=1.44+0.49*logL5100
            Tau=10.0**logRblr
            #print Tau*24*3600.0
            Rblr=Tau*24.0*3600.0*vel_light/3.08567758e13*10#parcecs
            Vel=np.sqrt(bh_mass/fblr*G/Rblr/(3.08567758e16*(1e3)**2.0/1.9891e30))
            #Lha=L500**(1.010)*Ewha
            #Lhb=L500**(1.251)*Ewhb
            Lha=L500**(1.00)*Ewha*0.5#1.1#*0.5#/5100.0
            Lhb=L500**(1.00)*Ewhb*0.5#*1.1#*0.5#/5100.0
            #print np.log10(Lha*1e44/1e-16/dla)
            L500t=L500*alpha*(1.0+np.sin((mjd-mjd0-Tau)*np.pi/tp))*beta/dla
            Lhat=Lha*(1.0+alpha*(1.0+np.sin((mjd-mjd0)*np.pi/tp)))
            Lhbt=Lhb*(1.0+alpha*(1.0+np.sin((mjd-mjd0)*np.pi/tp)))
            dlam_hb=(Vel/vel_light)*4861.0
            dlam_ha=(Vel/vel_light)*6563.0
            Lhbtf=(Lhbt/np.sqrt(2.0*np.pi*dlam_hb**2.0))*np.exp(-0.5*((wave-4861.0)/dlam_hb)**2.0)#*dwave
            Lhatf=(Lhat/np.sqrt(2.0*np.pi*dlam_ha**2.0))*np.exp(-0.5*((wave-6563.0)/dlam_ha)**2.0)#*dwave
            Spec=Lhbtf+Lhatf
            dlam=(1+z_0)
            
            fl=(wave/(5100.0*dlam))**(pw_a-2)*L500t
            
            Spec_s=shifts(Spec,wave,dlam)+fl#5100.0
            #print L500t/(4.0*np.pi*dl**2.0)*1e44/1e-16
            Spec_s=Spec_s/(4.0*np.pi*dl**2.0)*1e44/1e-16
            print "Delta L,vel",dlam_hb,Vel*2
            #print (L500-L500t)/L500
            #print 'L500 ='+str(L500/(4.0*np.pi*dl**2.0)*1e44/1e-16/dla)
            #print 'L500t ='+str(L500t/(4.0*np.pi*dl**2.0)*1e44/1e-16/dla)
            #print 'Lhat ='+str(Lhat/(4.0*np.pi*dl**2.0)*1e44/1e-16/dlam_ha/4.0)
            #print 'Lhbt ='+str(Lhbt/(4.0*np.pi*dl**2.0)*1e44/1e-16/dlam_hb/4.0)
            #print Spec_s[2000:3000]  
            L50044=logL5100+44
        else:
            Tau=0.0 
            bh_acret=0.0
            Spec_s=np.zeros(len(wave)) 
            Vel=0.0
            L50044=0.0
    else:
        Tau=0.0 
        bh_acret=0.0
        Spec_s=np.zeros(len(wave)) 
        Vel=0.0
        L50044=0.0
            
    #flux_f=flux+Spec_s
    #return flux_f
    return Spec_s,Tau,bh_mass,bh_acret,Vel,L50044
    
def fib_star_raw(outf,Av_g,mag_s,xx,yy,fib_id,dsep=1.0,lshf=5000,star_t='F0II (25291)',alp=0.0,bet=0.0,ha=0.0,ha1=0,dec_0=0.0,dirtemp="./",template4="sky.txt",template3="spEigenStar-55734.fits",template6="da012500_800.dat",wd=False,dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,vel_0=100.00,pdf=2,ifutype="SDSS",expt=900.0,outs=1,apot=0):
    #dsep=2.6
    vel_light=299792.458
    if "SDSS" in ifutype:
        scp_s=60.4#microns per arcsec
        fibA=210.0
        fibB=180.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "BOSS" in ifutype:
        scp_s=60.4#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "APO" in ifutype:
        scp_s=60.482#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "LCO" in ifutype:
        scp_s=91.475#microns per arcsec
        fibA=190.0
        fibB=180.0
        fibB=120.0
        obs_l=-29.0146
        sigma_inst=25.0
        beta_s=4.6#4.76
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=4.76
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        scp_s=300.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=4.76
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=4000.0
    else:
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0  
        obs_l=32.78   
        sigma_inst=25.0
        beta_s=4.76
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    #seeing=psf_pdf(1,psf_m=seeing)[0]           
    scalep=1.0/scp_s
    #print scalep,scp_s,beta_s
    dlam=(1.0+(vel_0/vel_light))
    Dfib=fibB*scalep
    n_pt=100000.0
    
    [sky_template,wave_s,crval_s,cdelt_s,crpix_s]=ssp_extract_sky(template4)
    #print band_mag(sky_template,crpix_s,cdelt_s,crval_s,dir='legacy/',k=3)+40
    mag_1=band_mag(sky_template*1e-18,crpix_s,cdelt_s,crval_s,dir='legacy/',k=3)
    #print mag_1+40
    dmag=-2.5*np.log10(10.0**(-0.4*(Flux_m))*np.pi*(fibB*scalep/2.0)**2.0)-mag_1
    mag_sk=mag_1+dmag
    dflux=10.0**(-0.4*dmag)*1e-18
    #print dflux/1e-16
    sky_template=sky_template*dflux/1e-16
    
    if wd == False:
        [ssp_template,wave,crval,cdelt,crpix]=ssp_extract_star(template3,star_type=star_t)
        #crval,cdelt,crpix
    else:
        [ssp_template,wave,crval,cdelt,crpix]=ssp_extract_wd(template6)
        #crval,cdelt,crpix
        star_t='DA White D.'
    mag_2=band_mag(ssp_template,crpix,cdelt,crval,dir='legacy/',k=3)
    dmag=mag_s-mag_2
    #print mag_2,dmag
    #mag_0=mag_0+dmag
    #mag_1=mag_1+dmag
    #mag_3=mag_3+dmag
    #mag_2=mag_2+dmag
    #mag_4=mag_4+dmag
    
    dflux=10.0**(-0.4*dmag)
    ssp_template=ssp_template*dflux/n_pt
    dust_rat_ssp=A_l(3.1,wave)
    
    
    cdelt_w=cdelt#1.25
    crval_w=3522.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    #if Flux_m != 20.0:
    #    Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    #else:
    #    Flux_m=20.0
    #s_nr=noise_sig(wave_f,SNi-1.0)
    nw=len(wave_f)
    nw_s=len(wave)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    #t_noise=1.0/Flux_m
    spect_t=np.zeros(nw_s)
    spect_ot=np.zeros(nw_s)
    spect=np.zeros(nw)
    noise=0.0#t_noise*ran.randn(nw)/s_nr*0.0
    Ft=0
    fac=1.0*expt/900.0
    
    n_d=100
    wave_d=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    wave_d=wave_d+(wave_d[1]-wave_d[0])/2.0
    wave_d=wave_d[0:n_d]
    wave_l=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    #R=refrac_dif(wave_d,ha,dec=0.0,phi=32.78,lo=5070.0,T=7.0,P=600.0,f=8.0)
    R=refrac_dif(wave_l,ha,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
    Rs=0.0
    Rs=refrac_dif(lshf,ha1,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
    pa=paralactic_angle(ha,dec=dec_0,phi=obs_l)
    R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
    R1=np.array([[np.cos((alp+bet)*np.pi/180.0),-np.sin((alp+bet)*np.pi/180.0)],[np.sin((alp+bet)*np.pi/180.0),np.cos((alp+bet)*np.pi/180.0)]])
    Rt=0.0
    if apot == 1:
        Rt=Rt+dsep*1e3*scalep
    R_t=np.dot(R1,np.array([0,Rt]))
    Xa=airmas(ha, dec=dec_0,phi=obs_l)
    print ha,ha1,Rs,pa,beta_s,seeing,Xa,fib_id
    #Rs=interp1d(wave_l,R,bounds_error=False,fill_value=0)(lshf)
    #print R[30]
    #phie=np.zeros(np.int(n_pt))#ran.randn(np.int(n_pt))*seeing/2.0
    #thee=np.zeros(np.int(n_pt))#ran.randn(np.int(n_pt))*seeing/2.0
    #phie=ran.randn(np.int(n_pt))*seeing/2.0
    #thee=ran.randn(np.int(n_pt))*seeing/2.0
    phie=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s) #ran.randn(np.int(n_pt))*seeing/2.0
    thee=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s) #ran.randn(np.int(n_pt))*seeing/2.0 
    
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    flux_lost=np.zeros(n_d+1)
    for w in range(0, n_d+1):
        R_adr=np.dot(R2,np.array([0,R[w]-Rs]))
        yoR=yo+R_adr[1]+R_t[1]
        xoR=xo+R_adr[0]+R_t[0]
        r=np.sqrt((xoR-phie)**2.0+(yoR-thee)**2.0)
        nt=np.where(r <= fibB*scalep/2.0)[0]
        flux_lost[w]=np.float(len(nt))
    #    nt_s=np.where((wave <= wave_l[w+1]) & (wave >= wave_l[w]))[0]
    #    if (w % 10) == 0:
    #        sycall('echo '+str(len(nt))+'  STARS'+' wl='+str(wave_d[w])+', '+str(Av_g))
    #    if len(nt) > 0:
    #        for k in range(0, len(nt)):
    #            Av=Av_g
    #            dust=10**(-0.4*Av*dust_rat_ssp)
    #            spect_s=ssp_template*dust/1e-16
    #            spect_t[nt_s]=spect_s[nt_s]+spect_t[nt_s]
    #yoR=yo+0
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    flux_lost=flux_lost/np.float(len(nt))
    #sycall('echo '+str(len(nt))+'  STARS comp')
    #print len(r),len(nt),fibB*scalep/2.0,fibB,scalep
    #print np.float(len(r))/np.float(len(nt))," Fiber-PSF factor",ifutype
    Av=Av_g
    dust=10**(-0.4*Av_g*dust_rat_ssp)
    if len(nt) > 0:
        for k in range(0, len(nt)):
            spect_s=ssp_template/1e-16*dust#/1e-16
            spect_o=ssp_template*dust
            spect_ot=spect_o+spect_ot 
            spect_t=spect_s+spect_t
    spect_t[np.isnan(spect_t)]=0
    spect_ot[np.isnan(spect_ot)]=0        
    spect_tf=shifts(spect_t,wave,dlam)
    spect_otf=shifts(spect_ot,wave,dlam)        
    spect_i=inst_disp(wave,spect_tf,sigma_inst)
    spect_oi=inst_disp(wave,spect_otf,sigma_inst)        
    spect_ii=spec_resol(wave,spect_i,sp_res)
    spect_oii=spec_resol(wave,spect_oi,sp_res)
    
    #dust_rat_sspf=A_l(3.1,wave_f)
    #dustf=10**(-0.4*Av_g*dust_rat_sspf)
    flux_lost_f=interp1d(wave_l,flux_lost,bounds_error=False,fill_value=0.)(wave_f)
    flux_lost_f=conv(flux_lost_f,ke=60)
    Kvl=extintion_c(wave_f)
    #Xa=airmas(ha, dec=dec_0,phi=obs_l)
    at_ext=10.0**(-0.4*Xa*Kvl)
            
    spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
    #spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
    specto=interp1d(wave,spect_oii,bounds_error=False,fill_value=0.)(wave_f)
    #specto=interp1d(wave,spect_ot,bounds_error=False,fill_value=0.)(wave_f)
    spect[np.isnan(spect)]=0
    specto[np.isnan(specto)]=0
    spect_ii[np.isnan(spect_ii)]=0
    #spect=spect
    specto=specto/1e-16#*dustf
    spect_sky=interp1d(wave_s,sky_template,bounds_error=False,fill_value=0.)(wave_f)
    spec_ifu=spect#*dustf
    
    
    #dust=10**(-0.4*Av_g*dust_rat_ssp)
    #print Av_g
    
    #mag_0=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=1)
    mag_0=band_mag(spect_ii*1e-16,crpix,cdelt,crval,dir='legacy/',k=1)
    mag_1=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=2)
    mag_2=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    mag_3=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=4)
    mag_4=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=5)
    mag_vec=np.array([mag_0,mag_1,mag_2,mag_3,mag_4])
    spec_ifu=spec_ifu*at_ext*flux_lost_f+spect_sky
    mag_f=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    #print mag_2
    #print 10**(-0.4*(mag_s-mag_2)),mag_s,mag_2,seeing
    #sys.exit()
    flux_d=flux_distort_corrvec(wave_f,xx,yy,fib_id)
    cal_v=response_v(wave_f,fac=fac)
    spec_ifuc=spec_ifu*cal_v*flux_d/1e-1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    spec_ifu_b=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_b)
    spec_ifu_r=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_r)    
    flatB=flat_b(wave_f_b)
    spec_ifu_b=spec_ifu_b*flatB
    flatR=flat_r(wave_f_r)
    spec_ifu_r=spec_ifu_r*flatR
    
    
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r)
    #nt1=np.where(np.isfinite(spec_b) ==  True)
    #nt2=np.where(np.isfinite(spec_r) ==  True)
    spec_b=spec_b/(wave_b/dwave_b/2.0/1633.0)
    spec_r=spec_r/(wave_r/dwave_r/2.0/2216.0)
    #spec_b=spec_b[nt1]
    #spec_r=spec_r[nt2] 
    
    #spec_ifu_e=noise
    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    ifu_o=np.zeros([nw,1])
    
    
    if pdf ==1 :
        max_val1=np.amax(spec_ifu)*1.25
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(star_t+' mf='+str(np.round(mag_f,3))+' mr='+str(np.round(mag_2,3))+' mrs='+str(np.round(mag_s,3))+' ms='+str(np.round(mag_sk,3)))
        plt.plot(wave_f,spec_ifu)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'.pdf')
        plt.close()
    
        max_val1=np.amax(specto)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_f,specto)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_O.pdf')
        plt.close()
        
        max_val1=np.amax(specto)*1.1
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[10^{-16} erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_f,spect_sky,alpha=0.4,color='grey')
        plt.plot(wave_f,spec_ifu,alpha=1.0,color='black')#spectro final
        plt.plot(wave_f,specto,color='blue',alpha=0.8)#spectro original
        plt.plot(wave_f,spect*at_ext*flux_lost_f,alpha=0.9,color='red')
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_On.jpg')#,dpi=1000)
        plt.close()
    
        max_val1=np.amax(spec_ifuc)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t)
        plt.plot(wave_f,spec_ifuc)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'C.pdf')
        plt.close()
    
        max_val1=np.amax(spec_b)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t+' Blue')
        plt.plot(pix_b,spec_b)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'B.pdf')
        plt.close()
        
        max_val1=np.amax(spec_r)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title(star_t+' Red')
        plt.plot(pix_r,spec_r)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'R.pdf')
        plt.close()
        
    if outs == 1:    
        ifu[:,0]=spec_ifu
        ifu_o[:,0]=specto
        h1=pyf.PrimaryHDU(ifu)
        h2=pyf.ImageHDU(ifu_e)
        h3=pyf.ImageHDU(ifu_1)
        h4=pyf.ImageHDU(ifu_m)
        h5=pyf.ImageHDU(ifu_o)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=nw
        h["NAXIS2"]=1
        h["COMMENT"]="Mock "+ifutype+" single spectra"
        h['CDELT1']=cdelt_w
        h['CRPIX1']=crpix_w
        h['CRVAL1']=crval_w
        h['CUNIT1']='Wavelength [A]'
        h['PSF']=seeing
        h['DFIB']=Dfib
        h['TYPE']=star_t
        h['REDSHIFT']=float(vel_0/vel_light)
        h['UNITS']='1E-16 erg/s/cm^2/A'
        hlist=pyf.HDUList([h1,h2,h3,h4,h5])
        hlist.update_extend()
        out_fit=dir_o+outf+'.fits'
        wfits_ext(out_fit,hlist)
        dir_o1=dir_o.replace(" ","\ ")
        out_fit1=dir_o1+outf+'.fits'
        sycall('gzip  '+out_fit1)
    
    return spec_b,spec_r,mag_vec


def fib_conv_raw(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,xx,yy,fib_id,bhmass,bhacre,mjd=56008,lshf=5000,Av_gal=0.0,ha=0,ha1=0,dec_0=0.0,ra_0=180.0,dirtemp="./",template4="sky.txt",outs=1,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="SDSS",expt=900.0):
    nh=dens
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol
    sup=4.0*np.pi*Rs**2.0
    vel_light=299792.458
    no_nan=0
    if "SDSS" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=210.0
        fibB=180.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "BOSS" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "APO" in ifutype:
        scp_s=60.482#microns per arcsec
        fibA=190.0
        fibB=120.0
        obs_l=32.78
        sigma_inst=25.0
        beta_s=2.0
=======
            sp_res=2000.0           
    scalep=1.0/scp_s
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=scalep
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    d_r=0.10/dkpcs
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600
    phi_g=phi_g*180/np.pi*3600
    Dfib=fibB*scalep
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template5)
    [ssp_template3,wave3,age_ssp3,met_ssp3,ml_ssp3,crval_w3,cdelt_w3,crpix_w3]=ssp_extract(template3)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)    
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g=asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    #pht_g=asosiate_pho2(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_s,met_s,age_s,x_g,y_g,z_g,x,y,z,Rs)
    
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=1.25
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    if Flux_m != 20.0:
        Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    else:
        Flux_m=20.0
    s_nr=noise_sig(wave_f,SNi-1.0)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    spec_val=np.zeros([18])
    n_ages=num_ages(age_ssp3)
    ages_r=arg_ages(age_ssp3)
    sim_imag=np.zeros([n_ages])
    t_noise=1.0/Flux_m
    con=0
    phie=phi+ran.randn(len(rad))*seeing
    thee=the+ran.randn(len(rad))*seeing
    phieg=phi_g+ran.randn(len(rad_g))*seeing 
    theeg=the_g+ran.randn(len(rad_g))*seeing
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    r_g=np.sqrt((xo-phieg)**2.0+(yo-theeg)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
    spect_t=np.zeros(nw_s)
    spect=np.zeros(nw)
    spect_g=np.zeros(nw_g)
    spect_gf=np.zeros(nw)
    noise=t_noise*ran.randn(nw)/s_nr
    mass_t=0
    ml_t=0
    vel_t=0
    sfr_t=0
    Av_s=0
    Av_sg=0
    sve_t=0
    Lt=0
    ml_ts=0
    met_ligt=0
    met_mas=0
    age_ligt=0
    age_mas=0
    Av_ligt=0
    Av_flux=0
    Ft=0
    sycall('echo '+str(len(nt))+'  STARS')
    if len(nt) > 0:
        mass_t=np.sum(mass_s[nt])
        vel_t=np.average(v_rad[nt])
        sve_t=np.std(v_rad[nt])
        mass_t_t=asosiate_ages(age_ssp3,age_s[nt],mass_s[nt])
        sim_imag=mass_t_t
        for k in range(0, len(nt)):
            nt_e=np.where((abs(phi[nt[k]]-phi_g) <= d_r) & (abs(the[nt[k]]-the_g) <= d_r) & (rad_g <= rad[nt[k]]))[0]#DECOMENTAR
            if len(nt_e) > 0:
                Av=np.sum(Av_g[nt_e])
            else:
                Av=0
            Av_s=Av+Av_s
            if np.isnan(in_ssp[nt[k]]):
                spect=spect
            else:
                if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                    dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                    spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust/1e-16#/20.0
                    spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                    spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                    spect_t=spect_sf+spect_t
                    ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]
                    Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]
                    Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                    met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                    met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                    age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                    age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas
                    Av_ligt=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+Av_ligt
                    Av_flux=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)+Av_flux
        ml_t=ml_t/len(nt)
        ml_ts=mass_t/Lt
        if Lt > 0:
            met_ligt=10.0**(met_ligt/Lt)
            age_ligt=10.0**(age_ligt/Lt)
            Av_ligt=(Av_ligt/Lt)
            Av_flux=(Av_flux/Ft)
            met_mas=10.0**(met_mas/mass_t)
            age_mas=10.0**(age_mas/mass_t)
        spect_t[np.isnan(spect_t)]=0
        spect_i=inst_disp(wave,spect_t,sigma_inst)
        spect_ii=spec_resol(wave,spect_i,sp_res)
        spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect[np.isnan(spect)]=0
        spec_val[0]=Av_s/len(nt)
    sycall('echo '+str(len(nt_g))+'  GAS')
    if len(nt_g) > 0:
        sfr_t=np.sum(sfri[nt_g])
        for k in range(0, len(nt_g)):
            if band_g[nt_g[k]] > 0:
                nt_e=np.where((abs(phi_g[nt_g[k]]-phi_g) <= d_r) & (abs(the_g[nt_g[k]]-the_g) <= d_r) & (rad_g <= rad_g[nt_g[k]]))[0]#DECOMENTAR
                if len(nt_e) > 0:
                    Av=np.sum(Av_g[nt_e])
                else:
                    Av=0
                Av_sg=Av+Av_sg
                if np.isnan(in_gas[nt_g[k]]):
                    spect_gf=spect_gf
                else:
                    if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:
                        dust=10**(-0.4*Av*dust_rat_gas)
                        spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#*0.01#*mass_g[nt_g[k]]
                        spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                        spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                        spect_g=spect_sfg+spect_g
        spec_val[6]=np.sum(Av_g[nt_g])         
        spec_val[4]=Av_sg/len(nt_g)
        spect_g[np.isnan(spect_g)]=0
        spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
        spect_g_ii=spec_resol(wave_g,spect_g_i,sp_res)
        spect_gf=interp1d(wave_g,spect_g_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect_gf[np.isnan(spect_gf)]=0
    spec_ifu=(spect_gf+spect+noise)
    spec_ifu_e=noise
    spec_val[1]=mass_t
    spec_val[2]=vel_t
    spec_val[3]=sfr_t
    spec_val[5]=spec_val[0]+spec_val[4]
    spec_val[7]=sve_t
    spec_val[8]=ml_t
    spec_val[10]=Lt
    spec_val[9]=ml_ts
    spec_val[11]=met_ligt
    spec_val[12]=met_mas
    spec_val[13]=age_ligt
    spec_val[14]=age_mas
    spec_val[15]=Av_ligt
    spec_val[16]=Ft
    spec_val[17]=Av_flux

    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    ifu_v=np.zeros([18,1])
    ifu_a=np.zeros([n_ages,1])

    max_val1=np.amax(spec_ifu)*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#0.75e-14)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[1E-16 erg/s/cm^2]$",fontsize=14)
    plt.plot(wave_f,spec_ifu)
    fig.tight_layout()
    plt.savefig(dir_o+outf+'.pdf')
    plt.close()
    

    ifu[:,0]=spec_ifu
    ifu_v[:,0]=spec_val
    ifu_a[:,0]=sim_imag
    ifu_e[:,0]=spec_ifu_e
                
           
             
    h1=pyf.PrimaryHDU(ifu)#.header
    h2=pyf.ImageHDU(ifu_e)
    h3=pyf.ImageHDU(ifu_1)
    h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=nw
    h["NAXIS2"]=1
    h["COMMENT"]="Mock "+ifutype+" single spectra"
    h['CDELT1']=cdelt_w
    h['CRPIX1']=crpix_w
    h['CRVAL1']=crval_w
    h['CUNIT1']='Wavelength [A]'
    h['PSF']=seeing
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='1E-16 erg/s/cm^2'
    hlist=pyf.HDUList([h1,h2,h3,h4])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
    ifu_v[1,0]=np.log10(ifu_v[1,0]+1.0)
    ifu_v[10,0]=np.log10(ifu_v[10,0]+1.0)
    ifu_v[15,0]=-2.5*np.log10(ifu_v[15,0]+0.0001)
    ifu_v[17,0]=-2.5*np.log10(ifu_v[17,0]+0.0001)
    h1t=pyf.PrimaryHDU(ifu_v)
    h=h1t.header
    h["NAXIS"]=2
    h["NAXIS1"]=18
    h["COMMENT"]="Real Values "+ifutype+" single spectra"
    h['Type0']=('DUST_S  ','Av')
    h['Type1']=('MASS    ','log10(Msun)')
    h['Type2']=('VEL     ','km/s')
    h['Type3']=('SFR     ','Msun/yr')
    h['Type4']=('DUST_G  ','Av')
    h['Type5']=('DUST_T  ','Av')
    h['Type6']=('DUST_Av ','Av')
    h['Type7']=('DISP    ','km/s')
    h['Type8']=('aML     ','Msun/Lsun')
    h['Type9']= ('tML     ','Msun/Lsun')
    h['Type10']=('LUM     ','log10(Lsun)')
    h['Type11']=('Z_lw    ','Z/H')
    h['Type12']=('Z_mw    ','Z/H')
    h['Type13']=('AGE_lw  ','Gyr')
    h['Type14']=('AGE_mw  ','Gyr')
    h['Type15']=('Av_lw   ','Mag')
    h['Type16']=('Flux    ','1e-16 ergs/s/cm2')
    h['Type17']=('Av_fw   ','Mag')
    h['PSF']=sig
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist1=pyf.HDUList([h1t])
    hlist1.update_extend()
    out_fit=dir_o+outf+'_val.fits'
    wfits_ext(out_fit,hlist1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val.fits'
    sycall('gzip -f '+out_fit1)
    h1tt=pyf.PrimaryHDU(ifu_a)
    h=h1tt.header
    h["NAXIS"]=2
    h["NAXIS1"]=n_ages 
    h["COMMENT"]="Real Values "+ifutype+" single spectra"
    h['PSF']=sig
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist1t=pyf.HDUList([h1tt])
    hlist1t.update_extend()
    out_fit=dir_o+outf+'_val_mass_t.fits'
    wfits_ext(out_fit,hlist1t)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val_mass_t.fits'
    sycall('gzip -f '+out_fit1)
    
   


def cube_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,sp_samp=1.25,sp_res=0.0,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=7,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="MaNGA"):
    #if not "MaNGA" in ifutype or not "CALIFA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    #print mass_gssp
    #print met_g
    #sys.exit()
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    if "MaNGA" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
<<<<<<< HEAD
    elif "LCO" in ifutype:
        scp_s=91.475#microns per arcsec
        fibA=190.0
        fibB=180.0
        fibB=120.0
        obs_l=-29.0146
        sigma_inst=25.0
        beta_s=4.6#76
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
<<<<<<< HEAD
        obs_l=32.78
        sigma_inst=25.0
        beta_s=4.76
=======
        nl=11
        sigma_inst=25.0
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        pix_s=0.2#0.1#0.025#arcsec
        scp_s=300.0#150.0#300.0#1200.0#microns per arcsec
        fibA=150.0
        fibB=120.0
<<<<<<< HEAD
        obs_l=32.78
        sigma_inst=25.0
        beta_s=4.76
=======
        nl=int(fov*scp_s/fibA/2)+1
        sigma_inst=25.0
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=4000.0
    else:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
<<<<<<< HEAD
        fibB=120.0  
        obs_l=32.78   
        sigma_inst=25.0
        beta_s=4.76
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0   
    #seeing=psf_pdf(1,psf_m=seeing)[0]        
    scalep=1.0/scp_s
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    cam_l=cam*(1+red_0)*(3.08567758e19*100)
    dap=scalep
=======
        fibB=120.0     
        sigma_inst=25.0
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi   
        if sp_res <=0:
            sp_res=2000.0       
    if sp_samp <= 0:
        sp_samp=5000./sp_res/2.0
    scalep=1.0/scp_s#2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=scalep
    xima=np.zeros(nl)
    yima=np.zeros(nl)
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    d_r=0.10/dkpcs
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
<<<<<<< HEAD
    the=the*180/np.pi*3600#/10.0
    phi=phi*180/np.pi*3600#/10.0
    #print np.amax(the),np.amin(the),cam
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600#/10.0
    phi_g=phi_g*180/np.pi*3600#/10.0
    Dfib=fibB*scalep
    fac=1.0*expt/900.0
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template5)
    [ssp_template3,wave3,age_ssp3,met_ssp3,ml_ssp3,crval_w3,cdelt_w3,crpix_w3]=ssp_extract(template3)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)    
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g=asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=1.25
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    #if Flux_m != 20.0:
    #    Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    #else:
    #    Flux_m=20.0
    #s_nr=noise_sig(wave_f,SNi-1.0)
    
    
    [sky_template,wave_s,crval_s,cdelt_s,crpix_s]=ssp_extract_sky(template4)
    mag_1=band_mag(sky_template*1e-18,crpix_s,cdelt_s,crval_s,dir='legacy/',k=3)
    dmag=-2.5*np.log10(10.0**(-0.4*(Flux_m))*np.pi*(fibB*scalep/2.0)**2.0)-mag_1
    mag_sk=mag_1+dmag
    dflux=10.0**(-0.4*dmag)*1e-18
    sky_spec=sky_template*dflux/1e-16
    
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    spec_val=np.zeros([18])
    n_ages=num_ages(age_ssp3)
    ages_r=arg_ages(age_ssp3)
    sim_imag=np.zeros([n_ages])
    #t_noise=1.0/Flux_m
    con=0
    #phie=phi+ran.randn(len(rad))*seeing/2.0
    #thee=the+ran.randn(len(rad))*seeing/2.0
    phie=phi+psf_moffat(len(rad),psf=seeing,beta=beta_s)
    thee=the+psf_moffat(len(rad),psf=seeing,beta=beta_s)
    #phieg=phi_g+ran.randn(len(rad_g))*seeing/2.0 
    #theeg=the_g+ran.randn(len(rad_g))*seeing/2.0
    phieg=phi_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
    theeg=the_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
    
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    
    
    #n_pt=1000
    n_pt=100000.0
    n_d=100
    #wave_d=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    #wave_d=wave_d+(wave_d[1]-wave_d[0])/2.0
    #wave_d=wave_d[0:n_d]
    wave_l=np.arange(n_d+1)/np.float(n_d)*(np.amax(wave_f)-np.amin(wave_f))+np.amin(wave_f)
    #R=refrac_dif(wave_d,ha,dec=0.0,phi=32.78,lo=5070.0,T=7.0,P=600.0,f=8.0)
    R=refrac_dif(wave_l,ha,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
    Rs=refrac_dif(lshf,ha1,dec=dec_0,phi=obs_l,lo=5070.0,T=7.0,P=600.0,f=8.0)
    pa=paralactic_angle(ha,dec=dec_0,phi=obs_l)
    R2=np.array([[np.cos(pa*np.pi/180.0),-np.sin(pa*np.pi/180.0)],[np.sin(pa*np.pi/180.0),np.cos(pa*np.pi/180.0)]])
    print ha,ha1,fib_id
    #Rs=interp1d(wave_l,R,bounds_error=False,fill_value=0)(lshf)
    #print R[30]
    #phie_f=ran.randn(np.int(n_pt))*seeing/2.0
    #thee_f=ran.randn(np.int(n_pt))*seeing/2.0    
    phie_f=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s)
    thee_f=psf_moffat(np.int(n_pt),psf=seeing,beta=beta_s)
    
    #xo=ran.randn(1)*0.025
    #yo=ran.randn(1)*0.025
    flux_lost=np.zeros(n_d+1)
    for w in range(0, n_d+1):        
        R_adr=np.dot(R2,np.array([0,R[w]-Rs]))
        yoR=yo+R_adr[1]
        xoR=xo+R_adr[0]#+Rt
        #yoR=yo+R[w]-Rs
        #xoR=xo
        r=np.sqrt((xoR-phie_f)**2.0+(yoR-thee_f)**2.0)
        nt1=np.where(r <= fibB*scalep/2.0)[0]
        flux_lost[w]=np.float(len(nt1))
    #    nt_s=np.where((wave <= wave_l[w+1]) & (wave >= wave_l[w]))[0]
    #    if (w % 10) == 0:
    #        sycall('echo '+str(len(nt))+'  STARS'+' wl='+str(wave_d[w])+', '+str(Av_g))
    #    if len(nt) > 0:
    #        for k in range(0, len(nt)):
    #            Av=Av_g
    #            dust=10**(-0.4*Av*dust_rat_ssp)
    #            spect_s=ssp_template*dust/1e-16
    #            spect_t[nt_s]=spect_s[nt_s]+spect_t[nt_s]
    yoR=yo
    xoR=xo
    r=np.sqrt((xoR-phie_f)**2.0+(yoR-thee_f)**2.0)
    nt1=np.where(r <= fibB*scalep/2.0)[0]
    flux_lost=flux_lost/np.float(len(nt1))
    flux_lost_f=interp1d(wave_l,flux_lost,bounds_error=False,fill_value=0.)(wave_f)
    
    
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    r_g=np.sqrt((xo-phieg)**2.0+(yo-theeg)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
    spect_t=np.zeros(nw_s)
    spect_to=np.zeros(nw_s)
    spect=np.zeros(nw)
    spect_g=np.zeros(nw_g)
    spect_go=np.zeros(nw_g)
    spect_gf=np.zeros(nw)
    #noise=t_noise*ran.randn(nw)/s_nr
    mass_t=0
    ml_t=0
    vel_t=0
    sfr_t=0
    Av_s=0
    Av_sg=0
    sve_t=0
    Lt=0
    ml_ts=0
    met_ligt=0
    met_mas=0
    age_ligt=0
    age_mas=0
    Av_ligt=0
    Av_flux=0
    Ft=0
    #print phi
    #print phie
    #print len(r)
    #print np.amin(r)
    #print np.amax(r)
    #print r
    #print fibB*scalep/2.0
    sycall('echo '+str(len(nt))+'  STARS')
    #if len(nt) > 500:#COMENTAR
    #    nt=nt[0:500]#COMENTAR
    if len(nt) > 0:
        mass_t=np.sum(mass_s[nt])
        #print np.log10(mass_t)
        vel_t=np.average(v_rad[nt])
        sve_t=np.std(v_rad[nt])
        mass_t_t=asosiate_ages(age_ssp3,age_s[nt],mass_s[nt])
        sim_imag=mass_t_t
        for k in range(0, len(nt)):
            nt_e=np.where((abs(phi[nt[k]]-phi_g) <= d_r) & (abs(the[nt[k]]-the_g) <= d_r) & (rad_g <= rad[nt[k]]))[0]#DECOMENTAR
            if len(nt_e) > 0:
                Av=np.sum(Av_g[nt_e])
            else:
                Av=0
            Av_s=Av+Av_s
            if np.isnan(in_ssp[nt[k]]):
                spect=spect
            else:
                if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                    #print k
                    dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                    spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust/1e-16#/20.0
                    spect_so=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33*dust/1e44
                    spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                    #spect_sfo=shifts(spect_so,wave,dlam[nt[k]])
                    spect_t=spect_sf+spect_t
                    #spect_to=spect_sfo+spect_to
                    spect_to=spect_so+spect_to
                    ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]
                    Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]
                    Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                    met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                    met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                    age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                    age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas
                    Av_ligt=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+Av_ligt
                    Av_flux=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)+Av_flux
        ml_t=ml_t/len(nt)
        ml_ts=mass_t/Lt
        if Lt > 0:
            met_ligt=10.0**(met_ligt/Lt)
            age_ligt=10.0**(age_ligt/Lt)
            Av_ligt=(Av_ligt/Lt)
            Av_flux=(Av_flux/Ft)
            met_mas=10.0**(met_mas/mass_t)
            age_mas=10.0**(age_mas/mass_t)
        spect_t[np.isnan(spect_t)]=0
        spect_to[np.isnan(spect_to)]=0
        spect_i=inst_disp(wave,spect_t,sigma_inst)
        #spect_io=inst_disp(wave,spect_to,sigma_inst)
        spect_ii=spec_resol(wave,spect_i,sp_res)
        #spect_iio=spec_resol(wave,spect_io,sp_res)
        spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
        #specto=interp1d(wave,spect_iio,bounds_error=False,fill_value=0.)(wave_f)
        specto=interp1d(wave,spect_to,bounds_error=False,fill_value=0.)(wave_f)
        spect[np.isnan(spect)]=0
        specto[np.isnan(specto)]=0
        spec_val[0]=Av_s/len(nt)
    else:
        ml_t=0
        ml_ts=0
        mass_t=0
        Av_flux=0
        Av_ligt=0
        Ft=0
        Lt=0
        met_mas=0
        age_mas=0
        met_ligt=0
        age_ligt=0
        sve_t=0
        specto=np.copy(wave_f)*0
        spect=np.copy(wave_f)*0
    sycall('echo '+str(len(nt_g))+'  GAS')
    if len(nt_g) > 0:
        sfr_t=np.sum(sfri[nt_g])
        for k in range(0, len(nt_g)):
            if band_g[nt_g[k]] > 0:
                nt_e=np.where((abs(phi_g[nt_g[k]]-phi_g) <= d_r) & (abs(the_g[nt_g[k]]-the_g) <= d_r) & (rad_g <= rad_g[nt_g[k]]))[0]#DECOMENTAR
                if len(nt_e) > 0:
                    Av=np.sum(Av_g[nt_e])
                else:
                    Av=0
                Av_sg=Av+Av_sg
                if np.isnan(in_gas[nt_g[k]]):
                    spect_gf=spect_gf
                else:
                    if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:
                        dust=10**(-0.4*Av*dust_rat_gas)
                        spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16*10.0**(-2.18+2.18-3.18)#+3.18)#+1.18)#4.18)#+0.3+0.6)#*0.01#*mass_g[nt_g[k]]
                        #spect_sgo=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)/1e-16*10.0**(-2.18+2.18-3.18)
                        spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                        #spect_sfgo=shifts(spect_sgo,wave_g,dlam_g[nt_g[k]])
                        #spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                        spect_g=spect_sfg+spect_g
                        #spect_go=spect_sfgo+spect_go
        spec_val[6]=np.sum(Av_g[nt_g])         
        spec_val[4]=Av_sg/len(nt_g)
        spect_g[np.isnan(spect_g)]=0
        #spect_go[np.isnan(spect_go)]=0
        spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
        #spect_go_i=inst_disp(wave_g,spect_go,sigma_inst)
        spect_g_ii=spec_resol(wave_g,spect_g_i,sp_res)
        #spect_go_ii=spec_resol(wave_g,spect_go_i,sp_res)
        spect_gf=interp1d(wave_g,spect_g_ii,bounds_error=False,fill_value=0.)(wave_f)
        #spect_gfo=interp1d(wave_g,spect_go_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect_gf[np.isnan(spect_gf)]=0
        #spect_gfo[np.isnan(spect_gfo)]=0
    else:
        sfr_t=0
        spect_gf=np.copy(wave_f)*0
    dust_rat_fin=A_l(3.1,wave_f)
    dust=10**(-0.4*Av_gal*dust_rat_fin)
    #bhmass=0.0 #this one desavilites the RM_spectra
    rm_spec,Tau,bh_mass,bh_acret,fhvel,l5100=RM_spectra(wave_f,specto,bhmass,bhacre,mjd=mjd,dl=cam_l,z_0=red_0)
    
    spec_ifu=(spect_gf+spect+rm_spec)*dust
    spec_ifu_o=spect_gf+spect#+rm_spec#spect_gfo+specto
    
    mag_0=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=1)
    mag_1=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=2)
    mag_2=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    mag_3=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=4)
    mag_4=band_mag(spec_ifu*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=5)
    mag_vec=np.array([mag_0,mag_1,mag_2,mag_3,mag_4])
    mag_f=band_mag(spec_ifu*flux_lost_f*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
    
    Kvl=extintion_c(wave_f)
    Xa=airmas(ha, dec=dec_0,phi=obs_l)
    at_ext=10.0**(-0.4*Xa*Kvl)
    spec_ifu=spec_ifu*at_ext*flux_lost_f
    #specto=specto/1e-16
    spect_sky=interp1d(wave_s,sky_spec,bounds_error=False,fill_value=0.)(wave_f)
    spec_ifu=spec_ifu+spect_sky
    flux_d=flux_distort_corrvec(wave_f,xx,yy,fib_id)
    cal_v=response_v(wave_f,fac=fac)
    spec_ifuc=spec_ifu*cal_v*flux_d/1e-1
    wave_f_b=np.arange(crval_w,6510.0,cdelt_w)
    wave_f_r=np.arange(5300.0,10352.0,cdelt_w)
    spec_ifu_b=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_b)
    spec_ifu_r=interp1d(wave_f,spec_ifuc,bounds_error=False,fill_value=0.)(wave_f_r)    
    flatB=flat_b(wave_f_b)
    spec_ifu_b=spec_ifu_b*flatB
    flatR=flat_r(wave_f_r)
    spec_ifu_r=spec_ifu_r*flatR
    
    
    [pix_b,wave_b,dwave_b]=ssp_extract_lambdpix(dir_tem=dirtemp,col="b")
    [pix_r,wave_r,dwave_r]=ssp_extract_lambdpix(dir_tem=dirtemp,col="r")
    spec_b=interp1d(wave_f_b,spec_ifu_b,bounds_error=False,fill_value=0.0)(wave_b)
    spec_r=interp1d(wave_f_r,spec_ifu_r,bounds_error=False,fill_value=0.0)(wave_r) 
#    nt1=np.where(np.isfinite(spec_b) ==  True)
#    nt2=np.where(np.isfinite(spec_r) ==  True)
    spec_b=spec_b/(wave_b/dwave_b/2.0/1633.0)
    spec_r=spec_r/(wave_r/dwave_r/2.0/2216.0)
    #spec_b=spec_b[nt1]
    #spec_r=spec_r[nt2]
    
    
    spec_val[1]=mass_t
    spec_val[2]=vel_t
    spec_val[3]=sfr_t
    spec_val[5]=spec_val[0]+spec_val[4]
    spec_val[7]=sve_t
    spec_val[8]=ml_t
    spec_val[10]=Lt
    spec_val[9]=ml_ts
    spec_val[11]=met_ligt
    spec_val[12]=met_mas
    spec_val[13]=age_ligt
    spec_val[14]=age_mas
    spec_val[15]=Av_ligt
    spec_val[16]=Ft
    spec_val[17]=Av_flux

    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    ifu_v=np.zeros([18,1])
    ifu_a=np.zeros([n_ages,1])
    
    if pdf == 1:
        max_val1=np.amax(spec_ifu)*1.25
        matplotlib.use('agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(' mf='+str(np.round(mag_f,3))+' mr='+str(np.round(mag_2,3))+' ms='+str(np.round(mag_sk,3)))
        plt.plot(wave_f,spec_ifu)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'.pdf')
        plt.close()
    
        max_val1=np.amax(spec_ifu_o)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
        #plt.title(star_t)
        plt.plot(wave_f,spec_ifu_o)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_O.pdf')
        plt.close()
        
        id_obj=target_idform(ra_0,dec_0,typ=1)
        max_val1=np.amax(spec_ifu_o)*1.25
        #=-np.amax(spec_ifu_o)*0.1
        fig, ax = plt.subplots(figsize=(6.5,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlim(3500,10100)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Flux $[10^{-16} erg/s/cm^2/\AA]$",fontsize=14)
        plt.title(id_obj)
        plt.plot(wave_f,spect_sky,alpha=0.4,color='grey')
        plt.plot(wave_f,spec_ifu,alpha=1.0,color='black')#spectro final
        plt.plot(wave_f,spect,color='blue',alpha=0.8)#spectro original
        plt.plot(wave_f,spect_gf,alpha=0.7,color='violet')
        plt.plot(wave_f,rm_spec,alpha=0.5,color='green')
        plt.plot(wave_f,(spect_gf+spect+rm_spec)*dust*at_ext*flux_lost_f,alpha=0.9,color='red')
        fig.tight_layout()
        plt.savefig(dir_o+outf+'_On.jpg')#,dpi=1000)
        plt.close()
    
        max_val1=np.amax(spec_ifuc)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("$Wavelength [A]$",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        #plt.title(star_t)
        plt.plot(wave_f,spec_ifuc)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'C.pdf')
        plt.close()
    
        max_val1=np.amax(spec_b)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title('Blue')
        plt.plot(pix_b,spec_b)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'B.pdf')
        plt.close()
        
        max_val1=np.amax(spec_r)*1.25
        fig, ax = plt.subplots(figsize=(8,5.5))
        ax.set_ylim(0,max_val1)
        ax.set_xlabel("Pixels",fontsize=14)
        ax.set_ylabel("Counts",fontsize=14)
        plt.title('Red')
        plt.plot(pix_r,spec_r)
        fig.tight_layout()
        plt.savefig(dir_o+outf+'R.pdf')
        plt.close() 
     
     
    if outs == 1: 
        ifu[:,0]=spec_ifu
        ifu_v[:,0]=spec_val
        ifu_a[:,0]=sim_imag
        ifu_e[:,0]=spec_ifu_o
        ifu_1[:,0]=rm_spec
        ifu_m[:,0]=spect         
        h1=pyf.PrimaryHDU(ifu)#.header
        h2=pyf.ImageHDU(ifu_e)
        h3=pyf.ImageHDU(ifu_1)
        h4=pyf.ImageHDU(ifu_m)
        h5=pyf.ImageHDU(ifu_e)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=nw
        h["NAXIS2"]=1
        h["COMMENT"]="Mock "+ifutype+" single spectra"
        h['CDELT1']=cdelt_w
        h['CRPIX1']=crpix_w
        h['CRVAL1']=crval_w
        h['CUNIT1']='Wavelength [A]'
        h['TAU']=(Tau,'RM Tau [days]')
        h['FVEL']=(fhvel,'Sigma_line [km/s]')
        h['BHMASS']=(np.log10(bh_mass+1.0),'BH mass [log10(Msun)]')
        h['BHACCR']=(bh_acret,'BH accret [Msun/yr]')
        h['L5100']=(l5100, 'logL5100 [erg/s]')
        h['PSF']=seeing
        h['DFIB']=Dfib
        h['CAMX']=0
        h['CAMY']=0
        h['CAMZ']=cam
        h['REDSHIFT']=float(red_0)
        h['H0']=ho
        h['Lambda_0']=Lam
        h['Omega_m']=Om
        h['UNITS']='1E-16 erg/s/cm^2'
        hlist=pyf.HDUList([h1,h2,h3,h4,h5])
        hlist.update_extend()
        out_fit=dir_o+outf+'.fits'
        wfits_ext(out_fit,hlist)
        dir_o1=dir_o.replace(" ","\ ")
        out_fit1=dir_o1+outf+'.fits'
        sycall('gzip  '+out_fit1)
        ifu_v[1,0]=np.log10(ifu_v[1,0]+1.0)
        ifu_v[10,0]=np.log10(ifu_v[10,0]+1.0)
        ifu_v[15,0]=-2.5*np.log10(ifu_v[15,0]+0.0001)
        ifu_v[17,0]=-2.5*np.log10(ifu_v[17,0]+0.0001)
        h1t=pyf.PrimaryHDU(ifu_v)
        h=h1t.header
        h["NAXIS"]=2
        h["NAXIS1"]=18
        h["COMMENT"]="Real Values "+ifutype+" single spectra"
        h['Type0']=('DUST_S  ','Av')
        h['Type1']=('MASS    ','log10(Msun)')
        h['Type2']=('VEL     ','km/s')
        h['Type3']=('SFR     ','Msun/yr')
        h['Type4']=('DUST_G  ','Av')
        h['Type5']=('DUST_T  ','Av')
        h['Type6']=('DUST_Av ','Av')
        h['Type7']=('DISP    ','km/s')
        h['Type8']=('aML     ','Msun/Lsun')
        h['Type9']= ('tML     ','Msun/Lsun')
        h['Type10']=('LUM     ','log10(Lsun)')
        h['Type11']=('Z_lw    ','Z/H')
        h['Type12']=('Z_mw    ','Z/H')
        h['Type13']=('AGE_lw  ','Gyr')
        h['Type14']=('AGE_mw  ','Gyr')
        h['Type15']=('Av_lw   ','Mag')
        h['Type16']=('Flux    ','1e-16 ergs/s/cm2')
        h['Type17']=('Av_fw   ','Mag')
        h['TAU']=(Tau,'RM Tau [days]')
        h['FVEL']=(fhvel,'Sigma_line [km/s]')
        h['BHMASS']=(np.log10(bh_mass+1.0),'BH mass [log10(Msun)]')
        h['BHACCR']=(bh_acret,'BH accret [Msun/yr]')
        h['L5100']=(l5100, 'logL5100 [erg/s]')
        h['PSF']=seeing
        h['DFIB']=Dfib
        h['CAMX']=0
        h['CAMY']=0
        h['CAMZ']=cam
        h['REDSHIFT']=float(red_0)
        h['H0']=ho
        h['Lambda_0']=Lam
        h['Omega_m']=Om
        hlist1=pyf.HDUList([h1t])
        hlist1.update_extend()
        out_fit=dir_o+outf+'_val.fits'
        wfits_ext(out_fit,hlist1)
        dir_o1=dir_o.replace(" ","\ ")
        out_fit1=dir_o1+outf+'_val.fits'
        sycall('gzip -f '+out_fit1)
        h1tt=pyf.PrimaryHDU(ifu_a)
        h=h1tt.header
        h["NAXIS"]=2
        h["NAXIS1"]=n_ages 
        h["COMMENT"]="Real Values "+ifutype+" single spectra"
        h['PSF']=seeing
        h['DFIB']=Dfib
        h['CAMX']=0
        h['CAMY']=0
        h['CAMZ']=cam
        h['REDSHIFT']=float(red_0)
        h['H0']=ho
        h['Lambda_0']=Lam
        h['Omega_m']=Om
        h['UNITS']='Msun'
        for kk in range(0, n_ages):
            h['AGE'+str(kk)]=ages_r[kk]
        hlist1t=pyf.HDUList([h1tt])
        hlist1t.update_extend()
        out_fit=dir_o+outf+'_val_mass_t.fits'
        wfits_ext(out_fit,hlist1t)
        dir_o1=dir_o.replace(" ","\ ")
        out_fit1=dir_o1+outf+'_val_mass_t.fits'
        sycall('gzip -f '+out_fit1)
    #if np.log10(bh_mass+1.0) > 6:
        print "BHMAS "+str(np.log10(bh_mass+1.0))+" Tau="+str(Tau)
        #if len(nt) > 3000:
        #    sys.exit()
    #else:
        #print "BHMAS "+str(np.log10(bh_mass+1.0))
    return spec_b,spec_r,mag_vec
    
def fib_star(outf,Av_g,mag_s,star_t='F0II (25291)',template3="spEigenStar-55734.fits",dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,vel_0=100.00,pdf=2,ifutype="SDSS"):
    vel_light=299792.458
    if "SDSS" in ifutype:
        scp_s=60.4#microns per arcsec
        fibA=210.0
        fibB=180.0
        sigma_inst=25.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    if "BOSS" in ifutype:
        scp_s=60.4#microns per arcsec
        fibA=190.0
        fibB=120.0
        sigma_inst=25.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        sigma_inst=25.0
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        scp_s=300.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=4000.0
    else:
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0           
    scalep=1.0/scp_s
    dlam=(1.0+(vel_0/vel_light))
    Dfib=fibB*scalep
    n_pt=1000.0
    [ssp_template,wave,crval,cdelt,crpix]=ssp_extract_star(template3,star_type=star_t)
    mag_1=band_mag(ssp_template,crpix,cdelt,crval,dir='legacy/',k=3)
    dmag=mag_s-mag_1
    dflux=10.0**(-0.4*dmag)
    ssp_template=ssp_template*dflux/n_pt
    dust_rat_ssp=A_l(3.1,wave)
    cdelt_w=cdelt#1.25
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    if Flux_m != 20.0:
        Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    else:
        Flux_m=20.0
    s_nr=noise_sig(wave_f,SNi-1.0)
    nw=len(wave_f)
    nw_s=len(wave)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    t_noise=1.0/Flux_m
    phie=ran.randn(np.int(n_pt))*seeing
    thee=ran.randn(np.int(n_pt))*seeing
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    spect_t=np.zeros(nw_s)
    spect=np.zeros(nw)
    noise=t_noise*ran.randn(nw)/s_nr
    Ft=0
    sycall('echo '+str(len(nt))+'  STARS')
    if len(nt) > 0:
        for k in range(0, len(nt)):
            Av=Av_g
            dust=10**(-0.4*Av*dust_rat_ssp)
            spect_s=ssp_template*dust/1e-16
            spect_s=spect_s+ran.randn(nw_s)*np.median(spect_s)*0.01
            spect_t=spect_s+spect_t
        spect_t[np.isnan(spect_t)]=0
        spect_tf=shifts(spect_t,wave,dlam)
        spect_i=inst_disp(wave,spect_tf,sigma_inst)
        spect_ii=spec_resol(wave,spect_i,sp_res)
        spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect[np.isnan(spect)]=0
    spec_ifu=spect+noise
    spec_ifu_e=noise
    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    max_val1=np.amax(spec_ifu)*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[1E-16 erg/s/cm^2/\AA]$",fontsize=14)
    plt.plot(wave_f,spec_ifu)
    fig.tight_layout()
    plt.savefig(dir_o+outf+'.pdf')
    plt.close()
    ifu[:,0]=spec_ifu
    h1=pyf.PrimaryHDU(ifu)
    h2=pyf.ImageHDU(ifu_e)
    h3=pyf.ImageHDU(ifu_1)
    h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=nw
    h["NAXIS2"]=1
    h["COMMENT"]="Mock "+ifutype+" single spectra"
    h['CDELT1']=cdelt_w
    h['CRPIX1']=crpix_w
    h['CRVAL1']=crval_w
    h['CUNIT1']='Wavelength [A]'
    h['PSF']=seeing
    h['DFIB']=Dfib
    h['REDSHIFT']=float(vel_0/vel_light)
    h['UNITS']='1E-16 erg/s/cm^2/A'
    hlist=pyf.HDUList([h1,h2,h3,h4])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
    
    
def fib_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,sp_res=0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="SDSS"):
    nh=dens
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol
    sup=4.0*np.pi*Rs**2.0
    vel_light=299792.458
    no_nan=0
    if "SDSS" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=210.0
        fibB=180.0
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "BOSS" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=190.0
        fibB=120.0
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        sigma_inst=25.0
        beta_s=4.7
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        pix_s=0.2#0.1#0.025#arcsec
        scp_s=300.0#150.0#300.0#1200.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
        beta_s=4.7
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=4000.0
    else:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
        beta_s=2.0
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0           
    scalep=1.0/scp_s
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=scalep
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    d_r=0.10/dkpcs
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600
    phi_g=phi_g*180/np.pi*3600
    Dfib=fibB*scalep
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template5)
    [ssp_template3,wave3,age_ssp3,met_ssp3,ml_ssp3,crval_w3,cdelt_w3,crpix_w3]=ssp_extract(template3)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)    
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g=asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    #pht_g=asosiate_pho2(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_s,met_s,age_s,x_g,y_g,z_g,x,y,z,Rs)
    
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=1.25
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    if Flux_m != 20.0:
        Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    else:
        Flux_m=20.0
    s_nr=noise_sig(wave_f,SNi-1.0)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    spec_ifu=np.zeros([nw])
    spec_ifu_e=np.zeros([nw])
    spec_val=np.zeros([18])
    n_ages=num_ages(age_ssp3)
    ages_r=arg_ages(age_ssp3)
    sim_imag=np.zeros([n_ages])
    t_noise=1.0/Flux_m
    con=0
    #phie=phi+ran.randn(len(rad))*seeing
    #thee=the+ran.randn(len(rad))*seeing
    #phieg=phi_g+ran.randn(len(rad_g))*seeing 
    #theeg=the_g+ran.randn(len(rad_g))*seeing
    
    phie=phi+psf_moffat(len(rad),psf=seeing,beta=beta_s)
    thee=the+psf_moffat(len(rad),psf=seeing,beta=beta_s)
    phieg=phi_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
    theeg=the_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
    
    xo=ran.randn(1)*0.025
    yo=ran.randn(1)*0.025
    r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
    r_g=np.sqrt((xo-phieg)**2.0+(yo-theeg)**2.0)
    nt=np.where(r <= fibB*scalep/2.0)[0]
    nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
    spect_t=np.zeros(nw_s)
    spect=np.zeros(nw)
    spect_g=np.zeros(nw_g)
    spect_gf=np.zeros(nw)
    noise=t_noise*ran.randn(nw)/s_nr
    mass_t=0
    ml_t=0
    vel_t=0
    sfr_t=0
    Av_s=0
    Av_sg=0
    sve_t=0
    Lt=0
    ml_ts=0
    met_ligt=0
    met_mas=0
    age_ligt=0
    age_mas=0
    Av_ligt=0
    Av_flux=0
    Ft=0
    sycall('echo '+str(len(nt))+'  STARS')
    if len(nt) > 0:
        mass_t=np.sum(mass_s[nt])
        vel_t=np.average(v_rad[nt])
        sve_t=np.std(v_rad[nt])
        mass_t_t=asosiate_ages(age_ssp3,age_s[nt],mass_s[nt])
        sim_imag=mass_t_t
        for k in range(0, len(nt)):
            nt_e=np.where((abs(phi[nt[k]]-phi_g) <= d_r) & (abs(the[nt[k]]-the_g) <= d_r) & (rad_g <= rad[nt[k]]))[0]#DECOMENTAR
            if len(nt_e) > 0:
                Av=np.sum(Av_g[nt_e])
            else:
                Av=0
            Av_s=Av+Av_s
            if np.isnan(in_ssp[nt[k]]):
                spect=spect
            else:
                if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                    dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                    spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust/1e-16#/20.0
                    spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                    spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                    spect_t=spect_sf+spect_t
                    ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]
                    Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]
                    Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                    met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                    met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                    age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                    age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas
                    Av_ligt=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+Av_ligt
                    Av_flux=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)+Av_flux
        ml_t=ml_t/len(nt)
        ml_ts=mass_t/Lt
        if Lt > 0:
            met_ligt=10.0**(met_ligt/Lt)
            age_ligt=10.0**(age_ligt/Lt)
            Av_ligt=(Av_ligt/Lt)
            Av_flux=(Av_flux/Ft)
            met_mas=10.0**(met_mas/mass_t)
            age_mas=10.0**(age_mas/mass_t)
        spect_t[np.isnan(spect_t)]=0
        spect_i=inst_disp(wave,spect_t,sigma_inst)
        spect_ii=spec_resol(wave,spect_i,sp_res)
        spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect[np.isnan(spect)]=0
        spec_val[0]=Av_s/len(nt)
    sycall('echo '+str(len(nt_g))+'  GAS')
    if len(nt_g) > 0:
        sfr_t=np.sum(sfri[nt_g])
        for k in range(0, len(nt_g)):
            if band_g[nt_g[k]] > 0:
                nt_e=np.where((abs(phi_g[nt_g[k]]-phi_g) <= d_r) & (abs(the_g[nt_g[k]]-the_g) <= d_r) & (rad_g <= rad_g[nt_g[k]]))[0]#DECOMENTAR
                if len(nt_e) > 0:
                    Av=np.sum(Av_g[nt_e])
                else:
                    Av=0
                Av_sg=Av+Av_sg
                if np.isnan(in_gas[nt_g[k]]):
                    spect_gf=spect_gf
                else:
                    if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:
                        dust=10**(-0.4*Av*dust_rat_gas)
                        spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#*0.01#*mass_g[nt_g[k]]
                        spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                        spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                        spect_g=spect_sfg+spect_g
        spec_val[6]=np.sum(Av_g[nt_g])         
        spec_val[4]=Av_sg/len(nt_g)
        spect_g[np.isnan(spect_g)]=0
        spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
        spect_g_ii=spec_resol(wave_g,spect_g_i,sp_res)
        spect_gf=interp1d(wave_g,spect_g_ii,bounds_error=False,fill_value=0.)(wave_f)
        spect_gf[np.isnan(spect_gf)]=0
    spec_ifu=(spect_gf+spect+noise)
    spec_ifu_e=noise
    spec_val[1]=mass_t
    spec_val[2]=vel_t
    spec_val[3]=sfr_t
    spec_val[5]=spec_val[0]+spec_val[4]
    spec_val[7]=sve_t
    spec_val[8]=ml_t
    spec_val[10]=Lt
    spec_val[9]=ml_ts
    spec_val[11]=met_ligt
    spec_val[12]=met_mas
    spec_val[13]=age_ligt
    spec_val[14]=age_mas
    spec_val[15]=Av_ligt
    spec_val[16]=Ft
    spec_val[17]=Av_flux

    ifu=np.zeros([nw,1])
    ifu_e=np.ones([nw,1])
    ifu_1=np.ones([nw,1])
    ifu_m=np.zeros([nw,1])
    ifu_v=np.zeros([18,1])
    ifu_a=np.zeros([n_ages,1])

    max_val1=np.amax(spec_ifu)*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#0.75e-14)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[1E-16 erg/s/cm^2]$",fontsize=14)
    plt.plot(wave_f,spec_ifu)
    fig.tight_layout()
    plt.savefig(dir_o+outf+'.pdf')
    plt.close()
    

    ifu[:,0]=spec_ifu
    ifu_v[:,0]=spec_val
    ifu_a[:,0]=sim_imag
    ifu_e[:,0]=spec_ifu_e
                
           
             
    h1=pyf.PrimaryHDU(ifu)#.header
    h2=pyf.ImageHDU(ifu_e)
    h3=pyf.ImageHDU(ifu_1)
    h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=nw
    h["NAXIS2"]=1
    h["COMMENT"]="Mock "+ifutype+" single spectra"
    h['CDELT1']=cdelt_w
    h['CRPIX1']=crpix_w
    h['CRVAL1']=crval_w
    h['CUNIT1']='Wavelength [A]'
    h['PSF']=seeing
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='1E-16 erg/s/cm^2'
    hlist=pyf.HDUList([h1,h2,h3,h4])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
    ifu_v[1,0]=np.log10(ifu_v[1,0]+1.0)
    ifu_v[10,0]=np.log10(ifu_v[10,0]+1.0)
    ifu_v[15,0]=-2.5*np.log10(ifu_v[15,0]+0.0001)
    ifu_v[17,0]=-2.5*np.log10(ifu_v[17,0]+0.0001)
    h1t=pyf.PrimaryHDU(ifu_v)
    h=h1t.header
    h["NAXIS"]=2
    h["NAXIS1"]=18
    h["COMMENT"]="Real Values "+ifutype+" single spectra"
    h['Type0']=('DUST_S  ','Av')
    h['Type1']=('MASS    ','log10(Msun)')
    h['Type2']=('VEL     ','km/s')
    h['Type3']=('SFR     ','Msun/yr')
    h['Type4']=('DUST_G  ','Av')
    h['Type5']=('DUST_T  ','Av')
    h['Type6']=('DUST_Av ','Av')
    h['Type7']=('DISP    ','km/s')
    h['Type8']=('aML     ','Msun/Lsun')
    h['Type9']= ('tML     ','Msun/Lsun')
    h['Type10']=('LUM     ','log10(Lsun)')
    h['Type11']=('Z_lw    ','Z/H')
    h['Type12']=('Z_mw    ','Z/H')
    h['Type13']=('AGE_lw  ','Gyr')
    h['Type14']=('AGE_mw  ','Gyr')
    h['Type15']=('Av_lw   ','Mag')
    h['Type16']=('Flux    ','1e-16 ergs/s/cm2')
    h['Type17']=('Av_fw   ','Mag')
    h['PSF']=sig
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist1=pyf.HDUList([h1t])
    hlist1.update_extend()
    out_fit=dir_o+outf+'_val.fits'
    wfits_ext(out_fit,hlist1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val.fits'
    sycall('gzip -f '+out_fit1)
    h1tt=pyf.PrimaryHDU(ifu_a)
    h=h1tt.header
    h["NAXIS"]=2
    h["NAXIS1"]=n_ages 
    h["COMMENT"]="Real Values "+ifutype+" single spectra"
    h['PSF']=sig
    h['DFIB']=Dfib
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist1t=pyf.HDUList([h1tt])
    hlist1t.update_extend()
    out_fit=dir_o+outf+'_val_mass_t.fits'
    wfits_ext(out_fit,hlist1t)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val_mass_t.fits'
    sycall('gzip -f '+out_fit1)
    
   


def cube_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,sp_samp=1.25,sp_res=0.0,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=7,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="MaNGA"):
    #if not "MaNGA" in ifutype or not "CALIFA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    #print mass_gssp
    #print met_g
    #sys.exit()
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    if "MaNGA" in ifutype:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0
        sigma_inst=25.0
        beta_s=2.0
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=2000.0
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
        nl=11
        sigma_inst=25.0
        beta_s=4.7
        if psfi <= 0:
            seeing=0.7
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=1700.0
    elif "MUSE" in ifutype:
        pix_s=0.2#0.1#0.025#arcsec
        scp_s=300.0#150.0#300.0#1200.0#microns per arcsec
        fibA=150.0
        fibB=120.0
        nl=int(fov*scp_s/fibA/2)+1
        sigma_inst=25.0
        beta_s=4.7
        if psfi <= 0:
            seeing=0.6
        else:
            seeing=psfi
        if sp_res <=0:
            sp_res=4000.0
    else:
        pix_s=0.5#arcsec
        scp_s=60.4#microns per arcsec
        fibA=150.0
        fibB=120.0     
        sigma_inst=25.0
        beta_s=2.0
        if psfi == 0:
            seeing=1.43
        else:
            seeing=psfi   
        if sp_res <=0:
            sp_res=2000.0       
    if sp_samp <= 0:
        sp_samp=5000./sp_res/2.0
    scalep=1.0/scp_s#2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=scalep
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    d_r=0.10/dkpcs
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600#+ran.randn(len(rad))*2.0
    phi=phi*180/np.pi*3600#+ran.randn(len(rad))*2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600#+ran.randn(len(rad_g))*2.0
    phi_g=phi_g*180/np.pi*3600#+ran.randn(len(rad_g))*2.0
    ns=3*nl*(nl-1)+1
    Dfib=fibA*scalep
    Rifu=Dfib*((2.0*nl-1.0)/2.0-0.5)
    xfib0=-Rifu
    yfib0=0
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    xifu=np.zeros(ns)
    yifu=np.zeros(ns)
    ini=0
    for i in range(0, nl):
        nt=nl*2-1-i
        yfib=yfib0+i*dyf*Dfib
        for j in range(0, nt):
            xfib=xfib0+(j*dxf+0.5*i)*Dfib
            xifu[ini]=xfib
            yifu[ini]=yfib
            ini=ini+1
        if i > 0:
            for j in range(0, nt):
                xfib=xfib0+(j*dxf+0.5*i)*Dfib
                xifu[ini]=xfib
                yifu[ini]=-yfib
                ini=ini+1
    ndt=35
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.0]
    dit[1,:]=[+0.00,+1.0]
    dit[2,:]=[+0.00,-1.0]
    dit[3,:]=[+0.25,+0.5]
    dit[4,:]=[+0.25,-0.5]
    dit[5,:]=[-0.25,+0.5]
    dit[6,:]=[-0.25,-0.5]
    dit[7,:]=[+0.50,+0.5]
    dit[8,:]=[+0.50,-0.5]
    dit[9,:]=[+0.50,+1.5]
    dit[10,:]=[+0.50,-1.5]
    dit[11,:]=[-0.50,+0.5]
    dit[12,:]=[-0.50,-0.5]
    dit[13,:]=[-0.50,+1.5]
    dit[14,:]=[-0.50,-1.5]
    dit[15,:]=[+1.00,+0.0]
    dit[16,:]=[+1.00,+0.5]
    dit[17,:]=[+1.00,-0.5]
    dit[18,:]=[+1.00,+1.0]
    dit[19,:]=[+1.00,-1.0]    
    dit[20,:]=[-1.00,+1.0]
    dit[21,:]=[-1.00,-1.0]
    dit[22,:]=[+1.50,+0.0]
    dit[23,:]=[+1.50,-0.5]
    dit[24,:]=[+1.50,+0.5]
    dit[25,:]=[+0.50,+0.0]
    dit[26,:]=[-0.50,+0.0]
    dit[27,:]=[+0.35,-0.9]
    dit[28,:]=[-0.35,-0.9]
    dit[29,:]=[+0.85,-1.3]
    dit[30,:]=[-0.85,-1.3]
    dit[31,:]=[+0.78,-0.25]
    dit[32,:]=[+0.78,+0.25]
    dit[33,:]=[-0.35,+0.73]
    dit[34,:]=[+0.50,+0.70]
    dyf=1.0
    dyt=Dfib/2.0/np.cos(30.0*np.pi/180.0)
    dxt=Dfib/2.0
    ndt=3
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.00]+ran.randn(2)*0.025
    dit[1,:]=[+0.00,+dyt/1.0]+ran.randn(2)*0.025
    dit[2,:]=[-dxt, +dyt/2.0]+ran.randn(2)*0.025
#    ndt=35
#    print Dfib/2.0/np.cos(30.0*np.pi/180.0)
#    import matplotlib.pyplot as plt
#    matplotlib.use('Agg')
    #plt.plot(phi,the,'o')
#    for i in range(0, ndt):
#        plt.plot(xifu+dit[i,0],yifu+dyf*dit[i,1],'o')
#    plt.show()
#    sys.exit()
#    Av=1.1
    
    #template="../../Base_bc03/templete_bc03.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template5)
    [ssp_template3,wave3,age_ssp3,met_ssp3,ml_ssp3,crval_w3,cdelt_w3,crpix_w3]=ssp_extract(template3)
#    [ssp_template5,wave5,age_ssp5,met_ssp5,ml_ssp5,crval_w5,cdelt_w5,crpix_w5]=ssp_extract(template5)
    ml_ssp=1.0/ml_ssp
    #ml_ssp5=1.0/ml_ssp5
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    #pht_g=asosiate_pho2(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_s,met_s,age_s,x_g,y_g,z_g,x,y,z,Rs)
    #print np.amax(pht_g[np.where(np.log10(pht_g) != 0)[0]]),"Q(H) max"
    #print np.amin(pht_g[np.where(np.log10(pht_g) != 0)[0]]),"Q(H) min"
    #print np.amax(met_g)/0.02,"Z_g max"
    #print np.amin(met_g)/0.02,"Z_g min"
    #print np.amax(temp_g),"T max"
    #print np.amin(temp_g),"T min"
    #print np.amax(nh),"nh max"
    #print np.amin(nh),"nh min"
    #print np.amax(met_s),"Z_s max"
    #print np.amin(met_s),"Z_s min"
    #print np.amax(met_ssp),"Z_ssp max"
    #print np.amin(met_ssp),"Z_ssp min"
    #print tem_gas[0]
    #print den_gas[0]
    #print met_gas[0]
    #print pht_gas[0]
    #sys.exit()
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=sp_samp
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    if Flux_m != 20.0:
        Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    else:
        Flux_m=20.0
    #print Flux_m
    #sys.exit()
    s_nr=noise_sig(wave_f,SNi-1.0)#15.0#35#55....70,105#140#80#30#15#60.0
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    spec_ifu=np.zeros([nw,ndt*ns])
    spec_ifu_e=np.zeros([nw,ndt*ns])
    spec_val=np.zeros([35,ndt*ns])
    n_ages=num_ages(age_ssp3)
    ages_r=arg_ages(age_ssp3)
    sim_imag=np.zeros([n_ages,ndt*ns])
    x_ifu=np.zeros(ndt*ns)
    y_ifu=np.zeros(ndt*ns)
    facto=(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)#*np.pi
    t_noise=1.0/Flux_m#0.007#2.0
    con=0
    for i in range(0, ndt):
        #print i
        #phie=phi+ran.randn(len(rad))*seeing# 1.43#/2.0
        #thee=the+ran.randn(len(rad))*seeing# 1.43#/2.0
        #phieg=phi_g+ran.randn(len(rad_g))*seeing #1.43#/2.0
        #theeg=the_g+ran.randn(len(rad_g))*seeing #1.43#/2.0
        
        phie=phi+psf_moffat(len(rad),psf=seeing,beta=beta_s)
        thee=the+psf_moffat(len(rad),psf=seeing,beta=beta_s)
        phieg=phi_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
        theeg=the_g+psf_moffat(len(rad_g),psf=seeing,beta=beta_s)
    
        for j in range(0, ns):
            #print i, j
            xo=xifu[j]+dit[i,0]
            yo=yifu[j]+dyf*dit[i,1]    
            r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
            r_g=np.sqrt((xo-phieg)**2.0+(yo-theeg)**2.0)
            #print np.amin(phi),np.amax(phi),np.amin(the),np.amax(the),Rifu
#            print np.amin(r),fibB*scalep/2.0,fibB,scalep
            nt=np.where(r <= fibB*scalep/2.0)[0]
            nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            noise=t_noise*ran.randn(nw)/s_nr#*0.0
            mass_t=0
            ml_t=0
            vel_t=0
            sfr_t=0
            Av_s=0
            Av_sg=0
            sve_t=0
            Lt=0
            Ltg=0
            ml_ts=0
            met_ligt=0
            met_mas=0
            age_ligt=0
            age_mas=0
            Av_ligt=0
            Av_flux=0
            Ve_ligt=0
            Ve_flux=0
            Avg_ligt=0
            Avg_flux=0
            Veg_ligt=0
            Veg_flux=0
            Ft=0
            Ftg=0
            Sig_flux=0
            Sig_ligt=0
            Sig_flux_g=0
            Sig_ligt_g=0
            wl_t=[]
            wf_t=[]
            wl_tg=[]
            wf_tg=[]
            va_1=[]
            va_1g=[]
            Mft=0
            age_flux=0
            age_Mflux=0
            met_ligt_g=0
            met_flux_g=0
#            noise2=t_noise/s_nr
#            plt.plot(wave_f,np.abs(noise*16.0))
#            plt.plot(wave_f,np.abs(noise2*16.0))
#            plt.show()
#            sys.exit()
            #print i,j,len(nt),"STARS"
            sycall('echo '+str(i)+'  '+str(j)+'  '+str(len(nt))+'  STARS')
            if len(nt) > 0:
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                sve_t=np.std(v_rad[nt])
                mass_t_t=asosiate_ages(age_ssp3,age_s[nt],mass_s[nt])
                sim_imag[:,con]=mass_t_t*facto
                for k in range(0, len(nt)):
                    #r_ge=np.sqrt((phi[nt[k]]-phi_g)**2.0+(the[nt[k]]-the_g)**2.0)
                    #nt_e=np.where((r_ge <= 0.05) & (rad_g <= rad[nt[k]]))[0]
                    #nt_e=np.where((r_g <= fibB*scalep/2.0) & (rad_g <= rad[nt[k]]))[0]
                    nt_e=np.where((abs(phi[nt[k]]-phi_g) <= d_r) & (abs(the[nt[k]]-the_g) <= d_r) & (rad_g <= rad[nt[k]]))[0]#DECOMENTAR
                    if len(nt_e) > 0:#DECOMENTAR
                        Av=np.sum(Av_g[nt_e])#DECOMENTAR
                    else:#DECOMENTAR
                        Av=0#DECOMENTAR
                    #Av=0#COMENTAR
                    Av_s=10**(-0.4*Av)+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        #print in_ssp[nt[k]]
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust/1e-16#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]#*20.0
                            Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]#/20.0
                            Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                            met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                            met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                            age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                            age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas                            
                            #if Av > 0:
                            ft_w=mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)
                            lt_w=mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]
                            fm_w=mass_s[nt[k]]*10**(-0.4*Av*0.44)
                            Mft=Mft+fm_w
                            age_flux=np.log10(age_s[nt[k]])*ft_w+age_flux
                            age_Mflux=np.log10(age_s[nt[k]])*fm_w+age_Mflux
                            Ve_ligt=v_rad[nt[k]]*lt_w+Ve_ligt
                            Ve_flux=v_rad[nt[k]]*ft_w+Ve_flux
                            Av_ligt=10**(-0.4*Av)*lt_w+Av_ligt
                            Av_flux=10**(-0.4*Av)*ft_w+Av_flux
                            va_1.extend([v_rad[nt[k]]])
                            wf_t.extend([ft_w])
                            wl_t.extend([lt_w])
                            #plt.plot(wave,ssp_template[in_ssp[nt[k]],:])
                            #plt.show()
                           # plt.plot(wave,spect_t)
                            #plt.show()
                            #print radL[nt[k]],mass_s[nt[k]],age_s[nt[k]],met_s[nt[k]]
                            #sys.exit()
                ml_t=ml_t/len(nt)
                ml_ts=mass_t/Lt
                if Lt > 0:
                    #met_ligt=10.0**(met_ligt/Lt)
                    met_ligt=met_ligt/Lt
                    age_ligt=10.0**(age_ligt/Lt)
                    age_flux=10.0**(age_flux/Ft)
                    age_Mflux=10.0**(age_Mflux/Mft)
                    Av_ligt=(Av_ligt/Lt)
                    Av_flux=(Av_flux/Ft)
                    Ve_ligt=(Ve_ligt/Lt)
                    Ve_flux=(Ve_flux/Ft)
                    #met_mas=10.0**(met_mas/mass_t)
                    met_mas=met_mas/mass_t
                    age_mas=10.0**(age_mas/mass_t)
                    va_1=np.array(va_1)
                    wf_t=np.array(wf_t)
                    wl_t=np.array(wl_t)
                    Sig_flux=np.sqrt(np.nansum(np.abs(wf_t)*(Ve_flux-va_1)**2.0)/(np.nansum(np.abs(wf_t))-np.nansum(wf_t**2.0)/np.nansum(np.abs(wf_t))))
                    Sig_ligt=np.sqrt(np.nansum(np.abs(wl_t)*(Ve_ligt-va_1)**2.0)/(np.nansum(np.abs(wl_t))-np.nansum(wl_t**2.0)/np.nansum(np.abs(wl_t))))
                spect_t[np.isnan(spect_t)]=0
                spect_i=inst_disp(wave,spect_t,sigma_inst)
                spect_ii=spec_resol(wave,spect_i,sp_res)
                spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
                spec_val[0,con]=Av_s/len(nt)#*facto
#                if j == 5:
#                    plt.plot(wave_f,spect*facto)
#                    plt.plot(wave_f,noise*facto)
#                    plt.savefig(dir_o+outf+'test5.pdf')
#                    plt.show()
#                    sys.exit()
            #print i,j,len(nt_g),"GAS"
            sycall('echo '+str(i)+'  '+str(j)+'  '+str(len(nt_g))+'  GAS')
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                #Ha_f=0
                #SFH_a=np.sum(sfri[nt_g])
                #Tem_a=np.average(temp_g[nt_g])
                #phi_a=np.average(pht_g[nt_g])
                #print SFH_a/(7.9e-42)
                for k in range(0, len(nt_g)):
                    if band_g[nt_g[k]] > 0:
                        #r_ge=np.sqrt((phi_g[nt_g[k]]-phi_g)**2.0+(the_g[nt_g[k]]-the_g)**2.0)
                        #nt_e=np.where((r_ge <= 0.05) & (rad_g <= rad_g[nt_g[k]]))[0]
                        #nt_e=np.where((r_g <= fibB*scalep/2.0) & (rad_g <= rad_g[nt_g[k]]))[0]
                        nt_e=np.where((abs(phi_g[nt_g[k]]-phi_g) <= d_r) & (abs(the_g[nt_g[k]]-the_g) <= d_r) & (rad_g <= rad_g[nt_g[k]]))[0]#DECOMENTAR
                        if len(nt_e) > 0:#DECOMENTAR
                            Av=np.sum(Av_g[nt_e])#DECOMENTAR
                        else:#DECOMENTAR
                            Av=0#DECOMENTAR
                    #    Av=0#COMENTAR
                        Av_sg=Av+Av_sg
                        if np.isnan(in_gas[nt_g[k]]):
                            spect_gf=spect_gf
                        else:
                            if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:
                                dust=10**(-0.4*Av*dust_rat_gas)
                                spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#*0.01#*mass_g[nt_g[k]]
                                spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                                lt_wg=np.nansum(gas_template[in_gas[nt_g[k]],:])*10.0**(-3.18)
                                ft_wg=np.nansum(spect_sfg)
                                Ltg=Ltg+lt_wg
                                Ftg=Ftg+ft_wg
                                Veg_ligt=v_rad_g[nt_g[k]]*lt_wg+Veg_ligt
                                Veg_flux=v_rad_g[nt_g[k]]*ft_wg+Veg_flux
                                Avg_ligt=10**(-0.4*Av)*lt_wg+Avg_ligt
                                Avg_flux=10**(-0.4*Av)*ft_wg+Avg_flux
                                met_ligt_g=np.log10(met_g[nt_g[k]])*lt_wg+met_ligt_g   
                                met_flux_g=np.log10(met_g[nt_g[k]])*ft_wg+met_flux_g                    
                                spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                                spect_g=spect_sfg+spect_g    
                                va_1g.extend([v_rad_g[nt_g[k]]])
                                wf_tg.extend([ft_wg])
                                wl_tg.extend([lt_wg])                            
                                #plt.xlim(3400,4000)
                            #plt.plot(wave_g,dust)
                            #plt.show()
                            #plt.xlim(3400,4000)                            
                            #plt.show()
                            #plt.plot(wave_g,gas_template[in_gas[nt_g[k]],:])
                            #plt.xlim(3400,4000)
                            #plt.plot(wave_g,spect_sg)
                            #plt.show()
                            #plt.xlim(3400,4000)
                            #plt.plot(wave_g,spect_g)
                            #plt.show()
                            #sys.exit()
                            #temp=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*mass_g[nt_g[k]]*3.846e33
                            #nt=np.where((wave_g < (6562+170)) & (wave_g > (6562-170)))[0]
                            #nt1=np.where((wave_g < (6262+170)) & (wave_g > (6262-170)))[0]
                            #Ha_f=np.sum(temp[nt])-np.sum(temp[nt1])+Ha_f
                #print Ha_f,(SFH_a/(7.9e-42)/Ha_f)**(1/3.)
                #print Tem_a,phi_a
#                plt.xlim(3400,4000)
#                plt.plot(wave_g,spect_g)   
                if Ltg > 0:
                    Avg_ligt=(Avg_ligt/Ltg)
                    Avg_flux=(Avg_flux/Ftg)
                    Veg_ligt=(Veg_ligt/Ltg)
                    Veg_flux=(Veg_flux/Ftg)
                    #met_ligt_g=10.0**(met_ligt_g/Ltg)
                    #met_flux_g=10.0**(met_flux_g/Ftg)
                    met_ligt_g=met_ligt_g/Ltg
                    met_flux_g=met_flux_g/Ftg
                    va_1g=np.array(va_1g)
                    wf_tg=np.array(wf_tg)
                    wl_tg=np.array(wl_tg)
                    Sig_flux_g=np.sqrt(np.nansum(np.abs(wf_tg)*(Veg_flux-va_1g)**2.0)/(np.nansum(np.abs(wf_tg))-np.nansum(wf_tg**2.0)/np.nansum(np.abs(wf_tg))))
                    Sig_ligt_g=np.sqrt(np.nansum(np.abs(wl_tg)*(Veg_ligt-va_1g)**2.0)/(np.nansum(np.abs(wl_tg))-np.nansum(wl_tg**2.0)/np.nansum(np.abs(wl_tg))))    
                spec_val[6,con]=np.sum(Av_g[nt_g])         
                spec_val[4,con]=Av_sg/len(nt_g)
                spect_g[np.isnan(spect_g)]=0
                spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
                spect_g_ii=spec_resol(wave_g,spect_g_i,sp_res)
                spect_gf=interp1d(wave_g,spect_g_ii,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
                #spect_gf=conv_sp(wave_f,spect_gf,ke=3)#coment this line
                #plt.plot(wave_f,spect_gf+spect)
                #plt.show()
                #sys.exit()
#            spec_ifu[:,con]=(spect_gf+spect)+noise#np.median(spect_gf+spect)*noise
#            spec_ifu[:,con]=np.sqrt((spect_gf+spect)**2.0+noise**2.0)            
            #spec_ifu[:,con]=np.sqrt((facto**2.0)*(spect_gf+spect)**2.0+noise**2.0)
            spec_ifu[:,con]=facto*(spect_gf+spect+noise)
            spec_ifu_e[:,con]=noise*facto
            spec_val[1,con]=mass_t*facto
            spec_val[2,con]=vel_t#*facto
            spec_val[3,con]=sfr_t*facto
            spec_val[5,con]=spec_val[0,con]*0+spec_val[4,con]
            spec_val[7,con]=sve_t
            spec_val[8,con]=ml_t#*facto
            spec_val[10,con]=Lt*facto
            spec_val[9,con]=ml_ts
            spec_val[11,con]=met_ligt
            spec_val[12,con]=met_mas
            spec_val[13,con]=age_ligt
            spec_val[14,con]=age_mas
            spec_val[15,con]=Av_ligt
            spec_val[16,con]=Ft*facto
            spec_val[17,con]=Av_flux
            spec_val[18,con]=Ve_ligt
            spec_val[19,con]=Ve_flux
            spec_val[20,con]=Sig_ligt
            spec_val[21,con]=Sig_flux
            spec_val[22,con]=Avg_ligt
            spec_val[23,con]=Avg_flux
            spec_val[24,con]=Veg_ligt
            spec_val[25,con]=Veg_flux
            spec_val[26,con]=Sig_ligt_g
            spec_val[27,con]=Sig_flux_g
            spec_val[28,con]=Ftg*facto
            spec_val[29,con]=Ltg*facto
            spec_val[30,con]=Mft*facto
            spec_val[31,con]=age_flux
            spec_val[32,con]=age_Mflux
            spec_val[33,con]=met_ligt_g
            spec_val[34,con]=met_flux_g
            #+ran.randn(nw)*np.median(spect)*0.05
            x_ifu[con]=xo
            y_ifu[con]=yo
            con=con+1
    nl=int(round((np.amax([np.amax(x_ifu),-np.amin(x_ifu),np.amax(y_ifu),-np.amin(y_ifu)])+1)*2/pix_s))
    print nl
    ifu=np.zeros([nw,nl,nl])
    ifu_e=np.ones([nw,nl,nl])
    ifu_1=np.ones([nw,nl,nl])
    ifu_m=np.zeros([nw,nl,nl])
    ifu_v=np.zeros([35,nl,nl])
    ifu_a=np.zeros([n_ages,nl,nl])
    xo=-nl/2*pix_s
    yo=-nl/2*pix_s
    print xo,yo
    xi=xo
    xf=xo
#    import matplotlib.pyplot as plt
                #plt.show()
    int_spect=np.zeros(nw)
    for i in range(0, nl):
        xi=xf
        xf=xf+pix_s
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+pix_s
            spt_new=np.zeros(nw)
            spt_err=np.zeros(nw)
            spt_val=np.zeros(35)
            spt_mas=np.zeros(n_ages)
            Wgt=0
            for k in range(0, len(x_ifu)):
                V1=np.sqrt((x_ifu[k]-xi)**2.0+(y_ifu[k]-yf)**2.0)
                V2=np.sqrt((x_ifu[k]-xf)**2.0+(y_ifu[k]-yf)**2.0)
                V3=np.sqrt((x_ifu[k]-xi)**2.0+(y_ifu[k]-yi)**2.0)
                V4=np.sqrt((x_ifu[k]-xf)**2.0+(y_ifu[k]-yi)**2.0)
                Vt=np.array([V1,V2,V3,V4])
                Rsp=np.sqrt((x_ifu[k]-(xf+xi)/2.0)**2.0+(y_ifu[k]-(yf+yi)/2.0)**2.0)
                Vmin=np.amin(Vt)
                Vmax=np.amax(Vt)
                #if Vmin <= fibB*scalep/2.0:
                #    if Vmax <= fibB*scalep/2.0:
                #        Wg=(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)
                #    else:
                #        Wg=(1.0-(Vmax-fibB*scalep/2.0)/(np.sqrt(2.0)*pix_s))*(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)
                #    spt_new=spec_ifu[:,k]*Wg+spt_new
                #    spt_err=(spec_ifu_e[:,k]*Wg)**2.0+spt_err**2.0
                #    Wgt=Wgt+Wg
                if Rsp <= fibA*scalep*1.4/2.0:
                    Wg=np.exp(-(Rsp/pix_s)**2.0/2.0)
                    spt_new=spec_ifu[:,k]*Wg+spt_new
                    spt_err=(spec_ifu_e[:,k]*Wg)**2.0+spt_err**2.0
                    spt_val=spec_val[:,k]*Wg+spt_val
                    spt_mas=sim_imag[:,k]*Wg+spt_mas
                    Wgt=Wgt+Wg
            if Wgt == 0:
                Wgt=1
            ifu[:,j,i]=spt_new/Wgt
            ifu_v[:,j,i]=spt_val/Wgt
            ifu_a[:,j,i]=spt_mas/Wgt
            #ifu_imag[:,j,i]=spt_imag/Wgt
            if np.sum(np.sqrt(spt_err/Wgt**2.0)) == 0:
                ifu_e[:,j,i]=1.0
            else:
                ifu_e[:,j,i]=np.sqrt(spt_err/Wgt**2.0)
               # ifu_m[:,j,i]=1.0
            int_spect=int_spect+spt_new/Wgt
                #plt.plot(wave,int_spect)
    #plt.plot(wave_f,int_spect)
    ##plt.show() 
    ##fig.tight_layout()
    #plt.savefig(dir_o+outf+'_int.pdf')
    #plt.close()
             
    h1=pyf.PrimaryHDU(ifu)#.header
    h2=pyf.ImageHDU(ifu_e)
    h3=pyf.ImageHDU(ifu_1)
    h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mock "+ifutype+" IFU"
    h["CRVAL1"]=0
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['R']=(sp_res,'Spectral Resolution')
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['IFUCON']=(str(np.int(ns))+' ','NFibers')
    h['UNITS']='1E-16 erg/s/cm^2'
    hlist=pyf.HDUList([h1,h2,h3,h4])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
    ifu_v[0,:,:]=-2.5*np.log10(ifu_v[0,:,:]+0.0001)
    ifu_v[1,:,:]=np.log10(ifu_v[1,:,:]+1.0)
    ifu_v[30,:,:]=np.log10(ifu_v[30,:,:]+1.0)
    ifu_v[10,:,:]=np.log10(ifu_v[10,:,:]+1.0)
    ifu_v[29,:,:]=np.log10(ifu_v[29,:,:]+1.0)
    ifu_v[15,:,:]=-2.5*np.log10(ifu_v[15,:,:]+0.0001)
    ifu_v[17,:,:]=-2.5*np.log10(ifu_v[17,:,:]+0.0001)
    ifu_v[22,:,:]=-2.5*np.log10(ifu_v[22,:,:]+0.0001)
    ifu_v[23,:,:]=-2.5*np.log10(ifu_v[23,:,:]+0.0001)
    h1t=pyf.PrimaryHDU(ifu_v)
    h=h1t.header
    h["NAXIS"]=3
    h["NAXIS3"]=35
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Real Values "+ifutype+" IFU"
    h["CRVAL1"]=0
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['Type0']=('Av_T    ','Mag')
    h['Type1']=('MASS    ','log10(Msun)')
    h['Type2']=('VEL     ','km/s')
    h['Type3']=('SFR     ','Msun/yr')
    h['Type4']=('DUST_G  ','Av BETA')
    h['Type5']=('DUST_T  ','Av BETA')
    h['Type6']=('DUST_Av ','Av BETA')
    h['Type7']=('DISP    ','km/s')
    h['Type8']=('aML     ','Msun/Lsun BETA')
    h['Type9']= ('tML     ','Msun/Lsun BETA')
    h['Type10']=('LUM     ','log10(Lsun)')
    h['Type11']=('Z_lw    ','log10(Z/H) add 1.77 to convert to log10(Z/Z_sun)')
    h['Type12']=('Z_mw    ','log10(Z/H) add 1.77 to convert to log10(Z/Z_sun)')
    h['Type13']=('AGE_lw  ','Gyr')
    h['Type14']=('AGE_mw  ','Gyr')
    h['Type15']=('Av_lw   ','Mag')
    h['Type16']=('FLUX    ','1e-16 ergs/s/cm2')
    h['Type17']=('Av_fw   ','Mag')
    h['Type18']=('VEL_lw  ','km/s')
    h['Type19']=('VEL_fw  ','km/s')
    h['Type20']=('DIS_lw   ','km/s BETA')
    h['Type21']=('DIS_fw   ','km/s BETA')
    h['Type22']=('Av_lw_g  ','Mag')
    h['Type23']=('Av_fw_g  ','Mag')
    h['Type24']=('VEL_lw_g ','km/s')
    h['Type25']=('VEL_fw_g ','km/s')
    h['Type26']=('DIS_l_gas','km/s BETA')
    h['Type27']=('DIS_f_gas','km/s BETA')
    h['Type28']=('FLUX_gas','1e-16 ergs/s/cm2 Bolometric')
    h['Type29']=('LUM_gas ','log10(Lsun) Bolometric')
    h['Type30']=('MASS_fw ','log10(Msun) BETA')
    h['Type31']=('AGE_fw  ','Gyr')
    h['Type32']=('AGE_mfw ','Gyr BETA')
    h['Type33']=('Z_lw_gas ','log10(Z/H) add 10.46 to convert to 12+log10(O/H) or add 1.77 to convert to log10(Z/Z_sun)')
    h['Type34']=('Z_fw_gas ','log10(Z/H) add 10.46 to convert to 12+log10(O/H) or add 1.77 to convert to log10(Z/Z_sun)')#8.69 
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['R']=(sp_res,'Spectral Resolution')
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist1=pyf.HDUList([h1t])
    hlist1.update_extend()
    out_fit=dir_o+outf+'_val.fits'
    wfits_ext(out_fit,hlist1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val.fits'
    sycall('gzip -f '+out_fit1)
    h1tt=pyf.PrimaryHDU(ifu_a)
    h=h1tt.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Real Values "+ifutype+" IFU"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist1t=pyf.HDUList([h1tt])
    hlist1t.update_extend()
    out_fit=dir_o+outf+'_val_mass_t.fits'
    wfits_ext(out_fit,hlist1t)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val_mass_t.fits'
    sycall('gzip -f '+out_fit1)
    
def photo_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    seeing=sig    
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    #radL=10.0*3.08567758e16*100.0
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600+ran.randn(len(rad))*seeing#/2.0
    phi=phi*180/np.pi*3600+ran.randn(len(rad))*seeing#/2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600+ran.randn(len(rad_g))*seeing#/2.0
    phi_g=phi_g*180/np.pi*3600+ran.randn(len(rad_g))*seeing#/2.0
#    print np.amax(phi),np.amin(np.abs(phi))
#    print np.amax(x),np.amin(np.abs(x))
#    sys.exit()
    Av=1.1
    #template="../home/sanchez/ppak/legacy/gsd61_156.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=30.0#5.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)#25000
    s_nr=noise_sig(wave_f,25)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    #spec_ifu=np.zeros([nw,ndt*ns])
    #spec_ifu_e=np.zeros([nw,ndt*ns])
    #x_ifu=np.zeros(ndt*ns)
    #y_ifu=np.zeros(ndt*ns)
    #t_noise=2.0
    #con=0
    #for i in range(0, ndt):
    #    for j in range(0, ns):
    #        xo=xifu[j]+dit[i,0]
    #        yo=yifu[j]+dyf*dit[i,1]    
    #        r=np.sqrt((xo-phi)**2.0+(yo-the)**2.0)
    #        r_g=np.sqrt((xo-phi_g)**2.0+(yo-the_g)**2.0)
    #        nt=np.where(r <= fibB*scalep/2.0)[0]
    #        nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
    #        spect_t=np.zeros(nw_s)
    #        spect=np.zeros(nw)
    #        spect_g=np.zeros(nw_g)
    #        spect_gf=np.zeros(nw)
    #        noise=t_noise*ran.randn(nw)/s_nr
            #if len(nt) > 0:
                #for k in range(0, len(nt)):
                #    if np.isnan(in_ssp[nt[k]]):
                #        spect=spect
                #    else:
                #        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                #            #dust=10**(-0.4*Av*dust_rat_ssp)
                #            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*1e10*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16#*dust
                #            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                #            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                #            spect_t=spect_sf+spect_t
                #spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                #spect[np.isnan(spect)]=0
            #if len(nt_g) > 0:
            #    for k in range(0, len(nt_g)):
            #        if np.isnan(in_gas[nt_g[k]]):
            #            spect_gf=spect_gf
            #        else:
            #            if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 157:
            #                dust=10**(-0.4*Av*dust_rat_gas)
            #                spect_sg=gas_template[in_gas[nt_g[k]],:]*10**(ha_gas[in_gas[nt_g[k]]])*sup[nt_g[k]]*fact[nt_g[k]]*band_g[nt_g[k]]/1e1/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16
            #                spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
            #                spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
            #                spect_g=spect_sfg+spect_g
            #    spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
            #    spect_gf[np.isnan(spect_gf)]=0
            #spec_ifu[:,con]=(spect_gf+spect)+noise
            #spec_ifu_e[:,con]=noise
            #x_ifu[con]=xo
            #y_ifu[con]=yo
            #con=con+1
    #import matplotlib.pyplot as plt
    #plt.plot(phi,the,'o')
    #plt.ylabel('some numbers')
    #plt.show()
    photo_imag=np.zeros([nw,nl,nl])
    ext_temp=np.zeros([3,nl,nl])
    #mas_temp=np.zeros([nl,nl])
    #vel_temp=np.zeros([nl,nl])
    #sfr_temp=np.zeros([nl,nl])
    #ra1=ran.randn(len(phi))*leng_s/dkpcs
    #de1=ran.randn(len(the))*leng_s/dkpcs
    #ra2=ran.randn(len(phi_g))*leng_s/dkpcs
    #de2=ran.randn(len(the_g))*leng_s/dkpcs
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            nt_g=np.where((phi_gt < xf) & (phi_gt >= xi) & (the_gt < yf) & (the_gt >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            #noise=t_noise*ran.randn(nw)/s_nr
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                for k in range(0, len(nt)):
                    nt_e=np.where(rad_g[nt_g] <= rad[nt[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
#                    Av=0
                    Av_s=Av+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            #print spect_sf
                            #print dlam[nt[k]]
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            #plt.plot(wave,ssp_template[in_ssp[nt[k]],:])
                            #plt.show()
                            #plt.plot(wave,spect_sf)
                            #plt.show()
                            #sys.exit()
                #ext_temp[0,j,i]=Av_s/len(nt)
                ext_temp[0,j,i]=mass_t
                ext_temp[1,j,i]=vel_t
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                for k in range(0, len(nt_g)):
#                    nt_e=np.where((phi_g < xf) & (phi_g >= xi) & (the_g < yf) & (the_g >= yi) & (rad_g <= rad_g[nt_g[k]]))[0]
                    nt_e=np.where(rad_g[nt_g] <= rad_g[nt_g[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
#                    Av=0
                    if np.isnan(in_gas[nt_g[k]]):
                        spect_gf=spect_gf
                    else:             
                        if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:   
                            dust=10**(-0.4*Av*dust_rat_gas)
                            spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#0.01#*mass_g[nt_g[k]]
                            spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                            spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                            spect_g=spect_sfg+spect_g  
                ext_temp[2,j,i]=sfr_t
                spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
            photo_imag[:,j,i]=spect+spect_gf
             
    h1=pyf.PrimaryHDU(photo_imag)#.header
    #h2=pyf.ImageHDU(ifu_e)
    #h3=pyf.ImageHDU(ifu_1)
    #h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dus_m=1
    if dus_m==1:
        h2=pyf.PrimaryHDU(ext_temp)
        hf=h2.header
        hf["NAXIS"]=3
        hf["NAXIS1"]=nl
        hf["NAXIS2"]=nl
        hf["NAXIS3"]=3
        hf["CRVAL1"]=0#oap
        hf["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
        hf["CRPIX1"]=nl/2
        hf["CTYPE1"]='RA---TAN'
        hf["CRVAL2"]=0#oap
        hf["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
        hf["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CRPIX2"]=nl/2
        hf["CTYPE2"]='DEC--TAN'
        hf['CUNIT1']='deg     '                                           
        hf['CUNIT2']='deg     '
        hf['RADECSYS']='ICRS    '
        hf['SYSTEM']='FK5     '
        hf['EQUINOX']=2000.00
       # hf['Type1']=('DUST    ','Av')
        hf['Type1']=('MASS    ','Msun')
        hf['Type2']=('VEL     ','km/s')
        hf['Type3']=('SFR     ','Msun/yr')
        hf['PSF']=sig
        hf['FOV']=fov
        hf['KPCSEC']=dkpcs
        hf['CAMX']=observer[0]#0
        hf['CAMY']=observer[1]#0
        hf['CAMZ']=observer[2]#cam
        hf['REDSHIFT']=float(red_0)
        hf['H0']=ho
        hf['Lambda_0']=Lam
        hf['Omega_m']=Om
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        out_fit_d=dir_o+outf+'_val.fits'
        wfits_ext(out_fit_d,hlist)
        dir_od1=dir_o.replace(" ","\ ")
        out_fitd1=dir_od1+outf+'_val.fits'
        sycall('gzip -f '+out_fitd1)
    #band_photo_r(photo_imag, h, name=outf+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir_o)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)

def sim_conv_vel(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,x2,y2,z2,vx2,vy2,vz2,x_g2,y_g2,z_g2,vx_g2,vy_g2,vz_g2,dir_o='',red_0=0.01,red_02=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0],observer2=[0,0,0],ang=0.0):
    vel_light=299792.458
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    cam2=cd.comoving_distance(red_02, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    r_disk=np.sqrt(x**2.+y**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    v_r_disk=(vx*x+vy*y)/r_disk#radial disk velocity   
    v_t_disk=(-vx*y+vy*x)/r_disk#tangential disk velocity
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    r_g_disk=np.sqrt(x_g**2.+y_g**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    v_r_g_disk=(vx_g*x_g+vy_g*y_g)/r_g_disk#radial disk velocity
    v_t_g_disk=(-vx_g*y_g+vy_g*x_g)/r_g_disk#tangential disk velocity
    
    v_r_disk_x=v_r_disk*x/r_disk
    v_r_disk_y=v_r_disk*y/r_disk
    v_r_disk_z=np.copy(v_r_disk)*0.0
    v_t_disk_x=-v_t_disk*y/r_disk
    v_t_disk_y=v_t_disk*x/r_disk
    v_t_disk_z=np.copy(v_t_disk)*0.0
    v_z_disk_x=np.copy(vz)*0.0
    v_z_disk_y=np.copy(vz)*0.0
    v_z_disk_z=np.copy(vz)
    v_r_g_disk_x=v_r_g_disk*x_g/r_g_disk
    v_r_g_disk_y=v_r_g_disk*y_g/r_g_disk
    v_r_g_disk_z=np.copy(v_r_g_disk)*0.0
    v_t_g_disk_x=-v_t_g_disk*y_g/r_g_disk
    v_t_g_disk_y=v_t_g_disk*x_g/r_g_disk
    v_t_g_disk_z=np.copy(v_t_g_disk)*0.0
    v_z_g_disk_x=np.copy(vz_g)*0.0
    v_z_g_disk_y=np.copy(vz_g)*0.0
    v_z_g_disk_z=np.copy(vz_g)
    
    
    #R12=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]]) #z_rot
    #R12=np.array([[1,0,0],[0,np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0)],[0,np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0)]]) #x_rot
    R12=np.array([[np.cos(ang*np.pi/180.0),0,np.sin(ang*np.pi/180.0)],[0,1,0],[-np.sin(ang*np.pi/180.0),0,np.cos(ang*np.pi/180.0)]]) #y_rot
    R22=np.array([[1,0,0],[0,1,0],[0,0,1]]) 
    Ve2=np.array([v_r_disk_x,v_r_disk_y,v_r_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_r_disk_x=Vf2[0]
    v_r_disk_y=Vf2[1]
    v_r_disk_z=Vf2[2]
    Ve2=np.array([v_t_disk_x,v_t_disk_y,v_t_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_t_disk_x=Vf2[0]
    v_t_disk_y=Vf2[1]
    v_t_disk_z=Vf2[2]
    Ve2=np.array([v_z_disk_x,v_z_disk_y,v_z_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_z_disk_x=Vf2[0]
    v_z_disk_y=Vf2[1]
    v_z_disk_z=Vf2[2]
    Ve2=np.array([v_r_g_disk_x,v_r_g_disk_y,v_r_g_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_r_g_disk_x=Vf2[0]
    v_r_g_disk_y=Vf2[1]
    v_r_g_disk_z=Vf2[2]
    Ve2=np.array([v_t_g_disk_x,v_t_g_disk_y,v_t_g_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_t_g_disk_x=Vf2[0]
    v_t_g_disk_y=Vf2[1]
    v_t_g_disk_z=Vf2[2]
    Ve2=np.array([v_z_g_disk_x,v_z_g_disk_y,v_z_g_disk_z])
    Vf2=np.dot(np.dot(R22,R12),Ve2)
    v_z_g_disk_x=Vf2[0]
    v_z_g_disk_y=Vf2[1]
    v_z_g_disk_z=Vf2[2]
    
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    reds_g=reds_cos(rad_g/1e3)
    radA_g=rad_g/(1+reds_g)
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600
    phi_g=phi_g*180/np.pi*3600
    
    vxo=modes(vx2,nit=7,n_s=0.8)
    vyo=modes(vy2,nit=7,n_s=0.8)
    vzo=modes(vz2,nit=7,n_s=0.8)
    
    rad2=np.sqrt(x2**2.+y2**2.+(cam2-z2)**2.)
    r_disk2_t=np.sqrt(x2**2.+y2**2.)
    r_disk2_r=np.sqrt((x2*(cam2-z2))**2.+(y2*(cam2-z2))**2.+(x2**2.+y2**2.)**2.0)    
    dkpcs2=cam2*(1./3600.)*(np.pi/180.)
    v_rad2=(vx2*x2+vy2*y2+vz2*(z2-cam2))/rad2   
    v_r_disk2=(vx2*(-x2*(cam2-z2))+vy2*(-y2*(cam2-z2))+vz2*(x2**2.+y2**2.))/r_disk2_r#radial disk velocity projected
    v_t_disk2=(-vx2*y2+vy2*x2)/r_disk2_t#tangential disk velocity projected
    
    v_rad2_t=(v_t_disk_x*x2+v_t_disk_y*y2+v_t_disk_z*(z2-cam2))/rad2#tangential velociy projected on the line of sight
    v_rad2_r=(v_r_disk_x*x2+v_r_disk_y*y2+v_r_disk_z*(z2-cam2))/rad2#radial velociy projected on the line of sight
    v_rad2_z=(v_z_disk_x*x2+v_z_disk_y*y2+v_z_disk_z*(z2-cam2))/rad2#z velociy projected on the line of sight
    v_rad2_nc=v_rad2_z+v_rad2_r#non-circular velociy projected on the line of sight
    v_rad2_tt=v_rad2_z+v_rad2_r+v_rad2_t#total velociy projected on the line of sight
    
    vx2c=vx2-vxo
    vy2c=vy2-vyo
    vz2c=vz2-vzo
    #Velocity without proper motions
    v_rad2c=(vx2c*x2+vy2c*y2+vz2c*(z2-cam2))/rad2   
    v_r_disk2c=(vx2c*(-x2*(cam2-z2))+vy2c*(-y2*(cam2-z2))+vz2c*(x2**2.+y2**2.))/r_disk2_r#radial disk velocity projected
    v_t_disk2c=(-vx2c*y2+vy2c*x2)/r_disk2_t#tangential disk velocity projected
    
    rad_g2=np.sqrt(x_g2**2.+y_g2**2.+(cam2-z_g2)**2.)
    r_g_disk2_t=np.sqrt(x_g2**2.+y_g2**2.)
    r_g_disk2_r=np.sqrt((x_g2*(cam2-z_g2))**2.+(y_g2*(cam2-z_g2))**2.+(x_g2**2.+y_g2**2.)**2.0)    
    v_rad_g2=(vx_g2*x_g2+vy_g2*y_g2+vz_g2*(z_g2-cam2))/rad_g2 
    v_r_g_disk2=(vx_g2*(-x_g2*(cam2-z_g2))+vy_g2*(-y_g2*(cam2-z_g2))+vz_g2*(x_g2**2.+y_g2**2.))/r_g_disk2_r#radial disk velocity projected
    v_t_g_disk2=(-vx_g2*y_g2+vy_g2*x_g2)/r_g_disk2_t#tangential disk velocity projected

    v_rad_g2_t=(v_t_g_disk_x*x_g2+v_t_g_disk_y*y_g2+v_t_g_disk_z*(z_g2-cam2))/rad_g2#tangential velociy projected on the line of sight
    v_rad_g2_r=(v_r_g_disk_x*x_g2+v_r_g_disk_y*y_g2+v_r_g_disk_z*(z_g2-cam2))/rad_g2#radial velociy projected on the line of sight
    v_rad_g2_z=(v_z_g_disk_x*x_g2+v_z_g_disk_y*y_g2+v_z_g_disk_z*(z_g2-cam2))/rad_g2#z velociy projected on the line of sight
    v_rad_g2_nc=v_rad_g2_z+v_rad_g2_r#non-circular velociy projected on the line of sight
    v_rad_g2_tt=v_rad_g2_z+v_rad_g2_r+v_rad_g2_t#total velociy projected on the line of sight
    
    #Velocity without proper motions
    vx_g2c=vx_g2-vxo
    vy_g2c=vy_g2-vyo
    vz_g2c=vz_g2-vzo
    v_rad_g2c=(vx_g2c*x_g2+vy_g2c*y_g2+vz_g2c*(z_g2-cam2))/rad_g2 
    v_r_g_disk2c=(vx_g2c*(-x_g2*(cam2-z_g2))+vy_g2c*(-y_g2*(cam2-z_g2))+vz_g2c*(x_g2**2.+y_g2**2.))/r_g_disk2_r#radial disk velocity projected
    v_t_g_disk2c=(-vx_g2c*y_g2+vy_g2c*x_g2)/r_g_disk2_t#tangential disk velocity projected
    
    reds2=reds_cos(rad2/1e3)
    radA2=rad2/(1+reds2)
    reds_g2=reds_cos(rad_g2/1e3)
    radA_g2=rad_g2/(1+reds_g2)
    phi2=np.arcsin(x2/radA2)
    the2=np.arcsin(y2/(radA2*np.cos(phi2)))
    the2=the2*180/np.pi*3600
    phi2=phi2*180/np.pi*3600
    phi_g2=np.arcsin(x_g2/radA_g2)
    the_g2=np.arcsin(y_g2/(radA_g2*np.cos(phi_g2)))
    the_g2=the_g2*180/np.pi*3600
    phi_g2=phi_g2*180/np.pi*3600
    
    n_vel=24
    sim_vel=np.zeros([n_vel,nl,nl])
    phit=phi
    the1=the
    phi_gt=phi_g
    the_gt=the_g
    
    phit2=phi2
    the12=the2
    phi_gt2=phi_g2
    the_gt2=the_g2
    
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            if len(nt) > 0:
                vel_tan=np.nanmean(v_t_disk[nt])
                vel_rad=np.nanmean(v_r_disk[nt])
                vel_alt=np.nanmean(vz[nt])
                sim_vel[0,j,i]=vel_tan
                sim_vel[1,j,i]=vel_rad
                sim_vel[2,j,i]=vel_alt   
            nt_g=np.where((phi_gt < xf) & (phi_gt >= xi) & (the_gt < yf) & (the_gt >= yi))[0]
            if len(nt_g) > 0:
                vel_g_tan=np.nanmean(v_t_g_disk[nt_g])
                vel_g_rad=np.nanmean(v_r_g_disk[nt_g])
                vel_g_alt=np.nanmean(vz_g[nt_g])
                sim_vel[3,j,i]=vel_g_tan
                sim_vel[4,j,i]=vel_g_rad
                sim_vel[5,j,i]=vel_g_alt
            nt2=np.where((phit2 < xf) & (phit2 >= xi) & (the12 < yf) & (the12 >= yi))[0]
            if len(nt2) > 0:
                vel_tan2=np.nanmean(v_t_disk2[nt2])
                vel_rad2=np.nanmean(v_r_disk2[nt2])
                vel_los2=np.nanmean(v_rad2[nt2])
                sim_vel[6,j,i]=vel_tan2
                sim_vel[7,j,i]=vel_rad2
                sim_vel[8,j,i]=vel_los2   
                
                vel_tan2=np.nanmean(v_t_disk2c[nt2])
                vel_rad2=np.nanmean(v_r_disk2c[nt2])
                vel_los2=np.nanmean(v_rad2c[nt2])
                sim_vel[12,j,i]=vel_tan2
                sim_vel[13,j,i]=vel_rad2
                sim_vel[14,j,i]=vel_los2

                vel_tan2=np.nanmean(v_rad2_t[nt2])
                vel_rad2=np.nanmean(v_rad2_nc[nt2])
                vel_los2=np.nanmean(v_rad2_tt[nt2])-np.nanmean(v_rad2c[nt2])
                sim_vel[18,j,i]=vel_tan2
                sim_vel[19,j,i]=vel_rad2
                sim_vel[20,j,i]=vel_los2
                                
            nt_g2=np.where((phi_gt2 < xf) & (phi_gt2 >= xi) & (the_gt2 < yf) & (the_gt2 >= yi))[0]
            if len(nt_g2) > 0:
                vel_g_tan2=np.nanmean(v_t_g_disk2[nt_g2])
                vel_g_rad2=np.nanmean(v_r_g_disk2[nt_g2])
                vel_g_los2=np.nanmean(v_rad_g2[nt_g2])
                sim_vel[9,j,i]=vel_g_tan2
                sim_vel[10,j,i]=vel_g_rad2
                sim_vel[11,j,i]=vel_g_los2

                vel_g_tan2=np.nanmean(v_t_g_disk2c[nt_g2])
                vel_g_rad2=np.nanmean(v_r_g_disk2c[nt_g2])
                vel_g_los2=np.nanmean(v_rad_g2c[nt_g2])
                sim_vel[15,j,i]=vel_g_tan2
                sim_vel[16,j,i]=vel_g_rad2
                sim_vel[17,j,i]=vel_g_los2 
                
                vel_g_tan2=np.nanmean(v_rad_g2_t[nt_g2])
                vel_g_rad2=np.nanmean(v_rad_g2_nc[nt_g2])
                vel_g_los2=np.nanmean(v_rad_g2_tt[nt_g2])-np.nanmean(v_rad_g2c[nt_g2])
                sim_vel[21,j,i]=vel_g_tan2
                sim_vel[22,j,i]=vel_g_rad2
                sim_vel[23,j,i]=vel_g_los2  
                              
    PSF=Gaussian2DKernel(stddev=3)
    for i in range(1,n_vel):
        if i != 3:
            sim_vel[i,:,:]=convolve(sim_vel[i,:,:], PSF)   
                 
    h1=pyf.PrimaryHDU(sim_vel)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_vel
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mook Ilustris IFS Mass"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['KPCSEC2']=dkpcs2
    h['CAMX2']=observer2[0]#0
    h['CAMY2']=observer2[1]#0
    h['CAMZ2']=observer2[2]#cam
    h['REDSHIFT2']=float(red_02)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='km/s'
    h['THETA']=(ang ,'projected angle')
    h['Type0']=('Tangential Stellar Velocity','km/s')
    h['Type1']=('Radial Stellar Velocity','km/s')
    h['Type2']=('Z Stellar Velocity','km/s')
    h['Type3']=('Tangential Gas Velocity','km/s')
    h['Type4']=('Radial Gas Velocity','km/s')
    h['Type5']=('Z Gas Velocity','km/s')
    h['Type6']=('Projected Tangential Stellar Velocity','km/s')
    h['Type7']=('Projected Radial Stellar Velocity','km/s')
    h['Type8']=('Projected LOS Stellar Velocity','km/s')
    h['Type9']=('Projected Tangential Gas Velocity','km/s')
    h['Type10']=('Projected Radial Gas Velocity','km/s')
    h['Type11']=('Projected LOS Gas Velocity','km/s')
    h['Type12']=('Projected Tangential Stellar Velocity withot pm','km/s')
    h['Type13']=('Projected Radial Stellar Velocity withot pm','km/s')
    h['Type14']=('Projected LOS Stellar Velocity withot pm','km/s')
    h['Type15']=('Projected Tangential Gas Velocity withot pm','km/s')
    h['Type16']=('Projected Radial Gas Velocity withot pm','km/s')
    h['Type17']=('Projected LOS Gas Velocity withot pm','km/s')
    h['Type18']=('Projected Tangential Stellar Velocity on LOS','km/s')
    h['Type19']=('Projected Ncircular Stellar Velocity on LOS','km/s')
    h['Type20']=('Recidual Projected LOS Stellar Velocity','km/s')
    h['Type21']=('Projected Tangential Gas Stellar Velocity on LOS','km/s')
    h['Type22']=('Projected Ncircular Gas Stellar Velocity on LOS','km/s')
    h['Type23']=('Recidual Projected Gas LOS Stellar Velocity','km/s')
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    

def sim_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="../../Base_bc03/templete_bc03_2.fits",template="../home/sanchez/ppak/legacy/gsd61_156.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    #radL=10.0*3.08567758e16*100.0
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600#+ran.randn(len(rad))*1.43/2.0
    phi=phi*180/np.pi*3600#+ran.randn(len(rad))*1.43/2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600#+ran.randn(len(rad_g))*1.43/2.0
    phi_g=phi_g*180/np.pi*3600#+ran.randn(len(rad_g))*1.43/2.0
    
    #template="/home/hjibarram/FIT3D_py/Base_bc03/Base_bc17/bc17_salp_Agelin_Metlin_330.fits"
    
    #template2="templete_gas.fits"
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    [ssp_template1,wave1,age_ssp1,met_ssp1,ml_ssp1,crval_w1,cdelt_w1,crpix_w1]=ssp_extract(template2)
    ml_ssp1=1./ml_ssp1
    #for i in range(0, len(age_ssp1)):
    #    print age_ssp1[i],"  ",met_ssp1[i],"  ",ml_ssp1[i]
    #sys.exit()
    n_ages=num_ages(age_ssp)
    ages_r=arg_ages(age_ssp)
    sim_imag=np.zeros([n_ages,nl,nl])
    sim_imag2=np.zeros([n_ages,nl,nl])
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            if len(nt) > 0:
                #if i == 220 & j == 220:
                mass_t=asosiate_ages(age_ssp,age_s[nt],mass_s[nt])
                ligh_t=asosiate_light(age_ssp,age_ssp1,met_ssp1,ml_ssp1,age_s[nt],mass_s[nt],met_s[nt])
                sim_imag[:,j,i]=mass_t   
                sim_imag2[:,j,i]=ligh_t  
    #print np.sum(sim_imag2)   
    #print np.sum(sim_imag)
    #print np.sum(sim_imag)/np.sum(sim_imag2)
    #sys.exit()
    h1=pyf.PrimaryHDU(sim_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mook Ilustris IFS Mass"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
    
    h2=pyf.PrimaryHDU(sim_imag2)
    h=h2.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mook Ilustris IFS Light"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Lsun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist=pyf.HDUList([h2])
    hlist.update_extend()
    out_fit=dir_o+outf+'_L.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_L.fits'
    sycall('gzip -f '+out_fit1)
    
    
def photosim_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    dlam=(1+(v_rad/vel_light+reds*1.0))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    cdelt_w=30.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)
    nw=len(wave_f)
    nw_s=len(wave)
    #import matplotlib.pyplot as plt
    photosim_imag=np.zeros([nw,nl,nl])
    phit=phi
    the1=the
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                for k in range(0, len(nt)):
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)#*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_t=spect_sf+spect_t
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            photosim_imag[:,j,i]=spect
             
    h1=pyf.PrimaryHDU(photosim_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
    
def photosimextgas_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600
    phi_g=phi_g*180/np.pi*3600
    Av=1.1
    #template="../home/sanchez/ppak/legacy/gsd61_156.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=30.0#5.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)#25000
    s_nr=noise_sig(wave_f,25)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    #import matplotlib.pyplot as plt
    #plt.plot(phi,the,'o')
    #plt.ylabel('some numbers')
    #plt.show()
    photo_imag=np.zeros([nw,nl,nl])
    ext_temp=np.zeros([10,nl,nl])
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            nt_g=np.where((phi_gt < xf) & (phi_gt >= xi) & (the_gt < yf) & (the_gt >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                Lt=0
                Ft=0
                met_ligt=0
                met_mas=0
                age_ligt=0
                age_mas=0
                Av_ligt=0
                Av_flux=0
                for k in range(0, len(nt)):
                    nt_e=np.where(rad_g[nt_g] <= rad[nt[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
                    Av_s=Av+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            #ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]#*20.0
                            Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]#/20.0
                            Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                            met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                            met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                            age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                            age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas
                            Av_ligt=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+Av_ligt
                            Av_flux=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)+Av_flux
                if Lt > 0:       
                    ext_temp[0,j,i]=Av_ligt/Lt
                    ext_temp[9,j,i]=Av_flux/Ft
                    ext_temp[5,j,i]=10.0**(met_ligt/Lt)
                    ext_temp[7,j,i]=10.0**(age_ligt/Lt)
                if mass_t > 0:
                    ext_temp[6,j,i]=10.0**(met_mas/mass_t)
                    ext_temp[8,j,i]=10.0**(age_mas/mass_t)
                ext_temp[1,j,i]=mass_t
                ext_temp[2,j,i]=vel_t
                ext_temp[4,j,i]=Lt
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                for k in range(0, len(nt_g)):
                    nt_e=np.where(rad_g[nt_g] <= rad_g[nt_g[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
                    if np.isnan(in_gas[nt_g[k]]):
                        spect_gf=spect_gf
                    else:             
                        if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:   
                            dust=10**(-0.4*Av*dust_rat_gas)
                            spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust*10.0**(-2.18+2.18-1.18)#-3.18)#+0.3+0.6)#0.01#*mass_g[nt_g[k]]
                            spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                            spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                            spect_g=spect_sfg+spect_g  
                ext_temp[3,j,i]=sfr_t
                spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
            photo_imag[:,j,i]=spect+spect_gf
    ext_temp[1,:,:]=np.log10(ext_temp[1,:,:]+1.0)
    ext_temp[4,:,:]=np.log10(ext_temp[4,:,:]+1.0)
    ext_temp[0,:,:]=-2.5*np.log10(ext_temp[0,:,:]+0.0001)
    ext_temp[9,:,:]=-2.5*np.log10(ext_temp[9,:,:]+0.0001)
    h1=pyf.PrimaryHDU(photo_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dus_m=1
    if dus_m==1:
        h2=pyf.PrimaryHDU(ext_temp)
        hf=h2.header
        hf["NAXIS"]=3
        hf["NAXIS1"]=nl
        hf["NAXIS2"]=nl
        hf["NAXIS3"]=10
        hf["CRVAL1"]=0#oap
        hf["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
        hf["CRPIX1"]=nl/2
        hf["CTYPE1"]='RA---TAN'
        hf["CRVAL2"]=0#oap
        hf["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
        hf["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CRPIX2"]=nl/2
        hf["CTYPE2"]='DEC--TAN'
        hf['CUNIT1']='deg     '                                           
        hf['CUNIT2']='deg     '
        hf['RADECSYS']='ICRS    '
        hf['SYSTEM']='FK5     '
        hf['EQUINOX']=2000.00
        #hf['Type1']=('DUST    ','Av')
        #hf['Type2']=('MASS    ','Msun')
        #hf['Type3']=('VEL     ','km/s')
        #hf['Type4']=('SFR     ','Msun/yr')
        hf['Type0']=('Av  ','Mag')
        hf['Type1']=('MASS    ','log10(Msun)')
        hf['Type2']=('VEL     ','km/s')
        hf['Type3']=('SFR     ','Msun/yr')
        hf['Type4']=('LUM     ','log10(Lsun)')
        hf['Type5']=('Z_lw    ','Z/H')
        hf['Type6']=('Z_mw    ','Z/H')
        hf['Type7']=('AGE_lw  ','Gyr')
        hf['Type8']=('AGE_mw  ','Gyr') 
        hf['Type8']=('Av_flx  ','Mag')    
        #hf['PSF']=sig
        hf['FOV']=fov
        hf['KPCSEC']=dkpcs
        hf['CAMX']=observer[0]#0
        hf['CAMY']=observer[1]#0
        hf['CAMZ']=observer[2]#cam
        hf['REDSHIFT']=float(red_0)
        hf['H0']=ho
        hf['Lambda_0']=Lam
        hf['Omega_m']=Om
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        out_fit_d=dir_o+outf+'_val.fits'
        wfits_ext(out_fit_d,hlist)
        dir_od1=dir_o.replace(" ","\ ")
        out_fitd1=dir_od1+outf+'_val.fits'
        sycall('gzip -f '+out_fitd1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
def noise_sig(wave,sn):
    w1=3900
    w2=10100
    dw2=100
    dw1=100
    s=sn/(1.+np.exp(-(wave-w1)/dw1))/(1.+np.exp((wave-w2)/dw2))+1.0
    return s

def inst_disp(pdl_lamb,pdl_flux,sigma_inst):
    vel_light=299792.458
    sigma=sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    rsigma=sigma/dpix
    if sigma > 0.1:
        box=int(3.0*rsigma*2.0)
        if box < 3:
            box=3
        kernel=np.zeros(2*box+1)
        norm=0
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
            kernel[j]=gaus
            norm=norm+gaus
        kernel=kernel/norm
        pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
    else:
        pdl_flux_conv_inst=pdl_flux
    return pdl_flux_conv_inst

def spec_resol(pdl_lamb,pdl_flux,R):
    vel_light=299792.458
    #dl=5500.0/R
    #sigma=dl#sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    pdl_flux_conv_inst=np.zeros(len(pdl_lamb))
    for i in range(0, len(pdl_lamb)):
        sigma=pdl_lamb[i]/R
        rsigma=sigma/dpix
        if sigma > 0.1:
            box=int(3.0*rsigma*2.0*10.0)
            if box < 3:
                box=3
            nk=2*box+1
            kernel=np.zeros(nk)
            norm=0
            for j in range(0, 2*box+1):
                gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
                kernel[j]=gaus
                norm=norm+gaus
            kernel=kernel/norm
            if i < box:
                kernel=kernel[box+i:nk]
            if len(pdl_lamb)-1-i < box+1:
                kernel=kernel[0:len(pdl_lamb)-i+box]
            nlk=len(kernel)
            if i < box:
                pdl_flux_t=pdl_flux[0:nlk]
                pdl_lamb_t=pdl_lamb[0:nlk]
            elif len(pdl_lamb)-1-i < box+1:
                pdl_flux_t=pdl_flux[i-box:len(pdl_lamb)]
                pdl_lamb_t=pdl_lamb[i-box:len(pdl_lamb)]
            else:
                pdl_flux_t=pdl_flux[i-box:i-box+nlk]
                pdl_lamb_t=pdl_lamb[i-box:i-box+nlk]
            #print pdl_lamb_t.shape
            #print kernel.shape
            #print nlk
            #print pdl_lamb_t[nlk-1]
            pdl_flux_t=pdl_flux_t*kernel
            val=simpson_r(pdl_flux_t,pdl_lamb_t,0,nlk-2)/simpson_r(kernel,pdl_lamb_t,0,nlk-2)
            #pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
        else:
            val=pdl_flux[i]
            #pdl_flux_conv_inst=pdl_flux
        pdl_flux_conv_inst[i]=val
    return pdl_flux_conv_inst

def spec_resol_old(pdl_lamb,pdl_flux,R):
    vel_light=299792.458
    dl=5500.0/R
    sigma=dl#sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    rsigma=sigma/dpix
    if sigma > 0.1:
        box=int(3.0*rsigma*2.0)
        if box < 3:
            box=3
        kernel=np.zeros(2*box+1)
        norm=0
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
            kernel[j]=gaus
            norm=norm+gaus
        kernel=kernel/norm
        pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
    else:
        pdl_flux_conv_inst=pdl_flux
    return pdl_flux_conv_inst

def conv_sp(xt,yt,ke=2):
    nsf=len(xt)
    from scipy import signal
    krn=ke/(xt[1]-xt[0])
    ker=signal.gaussian(nsf, krn )
    ker=ker/np.sum(ker)
    ytmf=np.convolve(yt,ker,mode="same")
    return ytmf
     
def shifts(spect_s,wave,dlam):
    wave_f=wave*dlam
    spect_f=interp1d(wave_f, spect_s,bounds_error=False,fill_value=0.)(wave)
    nt=np.where(spect_f == 0)[0]
    #spect_f[nt]=np.nan
    return  spect_f 

def asosiate_gas(ssp_temp,wave,age_s,met_s,den_s,tem_s,ml_s,age,met,den,tem):
    met_s=float_(met_s)*0.02#127
    age_a=[]
    met_a=[]
    den_a=[]
    tem_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    ban2=0
    ban3=0
    age_t=sorted(age_s)
    met_t=sorted(met_s)
    den_t=sorted(den_s)
    tem_t=sorted(tem_s)
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    den_a.extend([den_t[0]])
    tem_a.extend([tem_t[0]])
    ind_ssp=np.zeros((len(met)),dtype=np.int)
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
        if den_t[i-1] > den_t[i]:
                ban2 =1
        if den_t[i-1] < den_t[i] and ban2 == 0:
            den_a.extend([den_t[i]])
        if tem_t[i-1] > tem_t[i]:
                ban3 =1
        if tem_t[i-1] < tem_t[i] and ban3 == 0:
            tem_a.extend([tem_t[i]])
    n_age=len(age_a)
    n_met=len(met_a)
    n_den=len(den_a)
    n_tem=len(tem_a)
    ind_age=ages_definition_l(age,age_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    ind_den=val_definition_l(den[ind_age[i][ind_met[j]]],den_a) 
                    for k in range(0, n_den):
                        if len(ind_den[k]) > 0:
                            ind_tem=val_definition_l(tem[ind_age[i][ind_met[j][ind_den[k]]]],tem_a)
                            for h in range(0, n_tem):
                                if len(ind_tem[h]) >  0:               
                                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]) & (den_s == den_a[k]) & (tem_s == tem_a[h]))[0]
                                    ind_ssp[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]=nt[0]
                                    #print nt[0]
                                    #print "QT=",age_s[nt[0]],";",age_a[i+1],age_a[i-1]
                                    #print age[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "Z=",met_s[nt[0]],";",met_a[j+1],met_a[j-1]
                                    #print met[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "nh=",den_s[nt[0]],";",den_a[k+1],den_a[k-1]
                                    #print den[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "Tem=",tem_s[nt[0]],";",tem_a[h+1],tem_a[h-1]
                                    #print tem[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #sys.exit()
    return ind_ssp

def asosiate_pho2(ssp_temp,wave,age_s,met_s,ml,massp_s,metp_s,agep_s,x_g,y_g,z_g,x_s,y_s,z_s,R_g):  
    mass=np.zeros(len(z_g))
    met=np.zeros(len(z_g))
    age=np.zeros(len(z_g))  
    for i in range(0, len(z_g)):
        r=np.sqrt((x_s-x_g[i])**2.0+(y_s-y_g[i])**2.0+(z_s-z_g[i])**2.0)
        nt=np.where(r <= R_g[i])[0]
        met[i]=np.nanmean(metp_s[nt])
        age[i]=np.nanmean(agep_s[nt])
        mass[i]=np.nansum(massp_s[nt])    
    vel_light=299792458.0
    h_p=6.62607004e-34
    #age=np.ones(len(met))*2.5e6/1e9
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    photo=np.zeros(len(age))
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
                    flux_0=ssp_temp[nt[0],:]/ml[nt[0]]/(h_p*vel_light/wave/1e-10/1e-7)*3.846e33
                    j1=0#int(0.47*n_c)
                    j2=int(0.63*len(wave))
                    norm=simpson_r(flux_0,wave,j1,j2)
                    photo[ind_age[i][ind_met[j]]]=norm*mass[ind_age[i][ind_met[j]]]+1#/(4.0*np.pi*Rs[ind_age[i][ind_met[j]]]**2.0*n_h[ind_age[i][ind_met[j]]])+1
    return np.log10(photo)


def asosiate_pho(ssp_temp,wave,age_s,met_s,ml,mass,met,Rs,n_h):
    vel_light=299792458.0
    h_p=6.62607004e-34
    age=np.ones(len(met))*2.5e6/1e9
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    photo=np.zeros(len(age))
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
                    flux_0=ssp_temp[nt[0],:]/ml[nt[0]]/(h_p*vel_light/wave/1e-10/1e-7)*3.846e33
                    j1=0#int(0.47*n_c)
                    j2=int(0.63*len(wave))
                    norm=simpson_r(flux_0,wave,j1,j2)
                    photo[ind_age[i][ind_met[j]]]=norm*mass[ind_age[i][ind_met[j]]]+1#/(4.0*np.pi*Rs[ind_age[i][ind_met[j]]]**2.0*n_h[ind_age[i][ind_met[j]]])+1
    return np.log10(photo)

def num_ages(age_s):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    n_age=len(age_a)
    return n_age

def reverse_numeric(x, y):
    return y - x

def arg_ages(age_s):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    return age_a[::-1]
  
def asosiate_ages(age_s,age,mass):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    mass_f=np.zeros(n_age)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            mass_f[i]=np.sum(mass[ind_age[i]])
    return mass_f[::-1]

def asosiate_light(age_s,age1_s,met1_s,ml1_s,age,mass,met):
    age_a=[]
    age1_a=[]
    met1_a=[]
    nssp=len(age_s)
    nssp1=len(age1_s)
    ban=0
    ban0=0
    ban1=0
    age_t=sorted(age_s)
    age1_t=sorted(age1_s)
    met1_t=met1_s#sorted(met_s)
    #ml_t=ml_t[np.argsort(age_s)]
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    age1_a.extend([age1_t[0]])
    met1_a.extend([met1_t[0]])
    for i in range(1, nssp1):
        if age1_t[i-1] > age1_t[i]:
            ban0 =1
        if age1_t[i-1] < age1_t[i] and ban0 == 0:
            age1_a.extend([age1_t[i]])
        if met1_t[i-1] > met1_t[i]:
            ban1 =1
        if met1_t[i-1] < met1_t[i] and ban1 == 0:
            met1_a.extend([met1_t[i]])
    ind_age=ages_definition_l(age,age_a)
    ind_age1=ages_definition_l(age1_s,age_a)
    n_age=len(age_a)
    n_age1=len(age1_a)
    n_met1=len(met1_a)
    ligh_f=np.zeros(n_age)
    for i in range(0, n_age):
#        print "Age=",age_a[i]
        if len(ind_age[i]) > 0:
            temp=0.0
            ind_age2=ages_definition_l(age[ind_age[i]],age1_a)
            for j in range(0, n_age1):
                if len(ind_age2[j]) > 0:   
                    ind_met=met_definition_l(met[ind_age[i][ind_age2[j]]],met1_a)
                    for k in range(0, n_met1):
                        if len(ind_met[k]) > 0:
                            nt=np.where((age1_s == age1_a[j]) & (met1_s == met1_a[k]))[0]
                            #temp=np.sum(mass[ind_age[i][ind_age2[j][ind_met[k]]]])+temp
                            temp=np.sum(mass[ind_age[i][ind_age2[j][ind_met[k]]]]/ml1_s[nt[0]])+temp
                            #print "MLr=",ml1_s[nt[0]],age1_s[nt[0]],met1_s[nt[0]]
            ligh_f[i]=temp
 #   print "_________________________________________"
    return ligh_f[::-1]
        
def asosiate_ssp(ssp_temp,wave,age_s,met_s,ml_s,age,met):
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    ind_ssp[:]=np.nan
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
            ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    #print n_age
    #print n_met
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
    return ind_ssp

def gas_extract(template):
    [pdl_flux_c_ini,hdr]=gdata(template, 0, header=True)
    [nf,n_c]=pdl_flux_c_ini.shape
    coeffs=np.zeros([nf,3])
    crpix=hdr["CRPIX1"]
    cdelt=hdr["CDELT1"]
    crval=hdr["CRVAL1"]
    tem_mod=[]
    pht_mod=[]
    den_mod=[]
    met_mod=[]
    Ha=[]
    name=[]
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_gas_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        TEM=data[3]
        PHT=data[2]
        DEN=data[1]
        MET=data[0]
        tem=float_(TEM.replace('t',''))
        pht=float_(PHT.replace('q',''))
        den=float_(DEN.replace('n',''))
        met=float_(MET.replace('z',''))
        tem_mod.extend([tem])    
        pht_mod.extend([pht])
        den_mod.extend([den])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr[header])
        Ha.extend([val_ml])
    wave_c=[]
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c.extend([(crval+cdelt*(j+1-crpix))])
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])
    wave_c=np.array(wave_c)
    return [pdl_flux_c_ini,wave_c,pht_mod,met_mod,den_mod,tem_mod,Ha,crval,cdelt,crpix]

def ssp_extract_arc(wave_i,dir_tem="",col='b'):
    file=dir_tem+"lamphgcdne.txt"
    f=open(file,'r')
    ints=[]
    wave=[]
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=filter(None,data)
            wave.extend([np.float(data[0])])
            ints.extend([np.float(data[1])])
    ints=np.array(ints)
    wave=np.array(wave)
    sgl=0.9
    emi=np.zeros(len(wave_i))
    for i in range(0, len(wave)):
        emi=np.exp(-0.5*((wave_i-wave[i])/sgl)**2.0)*ints[i]+emi
    emi=emi+5.0
    
    file=dir_tem+"output18129.fits"
    [flux, hdr0]=gdata(file, 0, header=True)
    l0=hdr0['CRVAL1']
    dl=hdr0['CDELT1']
    wave=l0+np.arange(len(flux))*dl
    flux_c=interp1d(wave,flux,bounds_error=False,fill_value=0.)(wave_i)
    emi=flux_c*0.0+emi#*10.0
    
    if 'b' in col:
        file_ar=dir_tem+'arc_blue_lib.txt' 
    if 'r' in col:
        file_ar=dir_tem+'arc_red_lib.txt' 
    wave=[]
    arc=[]
    ft=open(file_ar,'r')
    for line in ft:
        data=line.replace('\n','').split(',')
        data=filter(None,data)
        wave.extend([np.float(data[0])])
        arc.extend([np.float(data[1])])
    ft.close()
    arc=np.array(arc)
    wave=np.array(wave)
    flux_t=interp1d(wave,arc,bounds_error=False,fill_value=0.)(wave_i)
    emi=flux_t+emi
    return emi

def ssp_extract_lambdpix(dir_tem="",col="b"):
    if col == "b":
        file=dir_tem+"lamb_blue.txt"
        #file=dir_tem+"pixel_lamb_blue.txt"
    else:
        file=dir_tem+"lamb_red.txt"
        #file=dir_tem+"pixel_lamb_red.txt"
    f=open(file,'r')
    pix=[]
    wave=[]
    dwave=[]
    cont=0
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=filter(None,data)
            wave.extend([10.0**np.float(data[1])])
            pix.extend([np.float(data[0])])
            if cont == 0:
                dwave.extend([0])
            else:
                dwave.extend([wave[cont]-wave[cont-1]])
            cont=cont+1
    pix=np.array(pix)
    wave=np.array(wave)
    dwave=np.array(dwave)
    dwave[0]=dwave[1]
    return [pix,wave,dwave]

def FPS_conf(x,y,r1=7.4,r2=15.0):
    ri=r2-r1
    re=r1+r2
    r0=np.sqrt(x**2.0+y**2.0)
    if r0 <= re and r0 >= ri:
        alpha=(np.arctan2(y,x)-np.arccos((r1**2.0+r0**2.0-r2**2.0)/(2.0*r0*r1)))*180.0/np.pi % 360.0
        beta=(180.0-np.arccos((r1**2.0+r2**2.0-r0**2.0)/(2.0*r2*r1))*180.0/np.pi) % 360.0
    else:
        alpha=0.0
        beta=0.0
    alpha=np.round(alpha,5)
    beta=np.round(beta,5)
    return alpha,beta

def def_FPS(x,y,xs,ys,typ,r1=7.4,r2=15.0,fib2=5.0,Nob=400,Nsky=30,Nstd=70,Nobjt=2000):
    nt3=np.where(typ == 'Fiducial')[0]
    nt4=np.where(typ != 'Fiducial')[0]
    xr=xs[nt4]
    yr=ys[nt4]
    xf=xs[nt3]
    yf=ys[nt3]
    xa=np.copy(x)
    ya=np.copy(y)
    indt=np.arange(len(x))
    xb=[]
    yb=[]
    alp=[]
    bet=[]
    indo=[]
    flavt=[]
    sky=2
    std=2
    obj=2
    sky_c=0
    std_c=0
    obj_c=0
    #print Nobjt
    for i in range(0, len(xr)):
        tt=0
        r=np.sqrt((xr[i]-xa)**2.0+(yr[i]-ya)**2.0)
        nt=np.where((r <= (r1+r2)) & (r >= (r2-r1)))[0]
        if len(nt) > 0:
            #print len(nt),tt
            for j in range(0, len(nt)):
                if tt == 0:
                    x1=xa[nt[j]]
                    y1=ya[nt[j]]
                    r0=np.nanmin(np.sqrt((x1-xf)**2.0+(y1-yf)**2.0))
                    if r0 > (5.0+fib2/2.0):
                        t=ran.rand(100)[10]*100.0
                        #print t
                        #t=np.random.uniform(-0,100,1)[0]
                        if sky != 0:
                            if std != 0 and obj != 0:
                                if t < (Nsky*0.22)/np.float(Nsky+Nstd+Nob+700+100)*100.0:
                                    sky=1
                                    std=2
                                    obj=2
                            if std == 0 and obj != 0:
                                if t < Nsky/np.float(Nsky+Nob+700+100)*100.0:
                                    sky=1
                                    std=0
                                    obj=2
                            if std != 0 and obj == 0:
                                if t < Nsky/np.float(Nsky+Nstd+100)*100.0:
                                    sky=1
                                    std=2
                                    obj=0
                            if std == 0 and obj == 0:
                                sky=1
                                std=0
                                obj=0                                        
                        if std != 0:
                            if sky != 0 and obj != 0:
                                if t < (Nsky*0.22+Nstd*0.45)/np.float(Nsky+Nstd+Nob+700+100)*100.0 and t >= (Nsky*0.22)/np.float(Nsky+Nstd+Nob+700+100)*100.0:
                                    sky=2
                                    std=1
                                    obj=2  
                            if sky == 0 and obj != 0:
                                if t < (Nstd+100)/np.float(Nstd+Nob+700+100)*100.0:
                                    sky=0
                                    std=1
                                    obj=2
                            if sky != 0 and obj == 0:
                                if t >= Nsky/np.float(Nsky+Nstd+100)*100.0:
                                    sky=2
                                    std=1
                                    obj=0
                            if sky == 0 and obj == 0:
                                sky=0
                                std=1
                                obj=0 
                        if obj != 0:  
                            if sky != 0 and std != 0:                           
                                if t >= (Nsky*0.22+Nstd*0.45)/np.float(Nsky+Nstd+Nob+700+100)*100.0:
                                    sky=2
                                    std=2
                                    obj=1  
                            if sky == 0 and std !=0:
                                if t >= (Nstd+100)/np.float(Nstd+Nob+700+100)*100.0:
                                    sky=0
                                    std=2
                                    obj=1
                            if sky != 0 and std ==0:
                                if t >= Nsky/np.float(Nsky+Nob+700+100)*100.0:
                                    sky=2
                                    std=0
                                    obj=1
                            if sky == 0 and std ==0:
                                sky=0
                                std=0
                                obj=1    
                        #print sky,std,obj,j,indt[nt[j]]          
                        if sky == 1:
                            if indt[nt[j]] < Nobjt and indt[nt[j]] >= 0:
                                if i == 0:
                                    #print sky,std,obj,t,sky_c,std_c,obj_c,'sky'
                                    tt=1
                                    ind=j
                                    sky_c=sky_c+1
                                    flav=0
                                else:
                                    r00=np.nanmin(np.sqrt((x1-xb)**2.0+(y1-yb)**2.0))
                                    if r00 > (fib2/2.0):
                                        #print sky,std,obj,t,sky_c,std_c,obj_c,'sky'
                                        tt=1
                                        ind=j
                                        sky_c=sky_c+1
                                        flav=0
                            if sky_c == Nsky:
                                sky=0
                        if std == 1:
                            if indt[nt[j]] < Nobjt*2 and indt[nt[j]] >= Nobjt:
                                if i == 0:
                                    #print sky,std,obj,t,sky_c,std_c,obj_c,'std'
                                    tt=1
                                    ind=j
                                    std_c=std_c+1
                                    flav=1
                                else:
                                    r00=np.nanmin(np.sqrt((x1-xb)**2.0+(y1-yb)**2.0))
                                    #print r00,x1,y1
                                    #print xb
                                    if r00 > (fib2/2.0):
                                        #print sky,std,obj,t,sky_c,std_c,obj_c,'std'
                                        tt=1
                                        ind=j
                                        std_c=std_c+1
                                        flav=1
                            if std_c == Nstd:
                                std=0
                        if obj == 1:
                            if indt[nt[j]] >= Nobjt*2:
                                if i == 0:
                                    #print sky,std,obj,t,sky_c,std_c,obj_c,'obj'
                                    tt=1
                                    ind=j
                                    obj_c=obj_c+1
                                    flav=2
                                else:
                                    r00=np.nanmin(np.sqrt((x1-xb)**2.0+(y1-yb)**2.0))
                                    if r00 > (fib2/2.0):
                                        #print sky,std,obj,t,sky_c,std_c,obj_c,'obj'
                                        tt=1
                                        ind=j
                                        obj_c=obj_c+1
                                        flav=2
                            if obj_c == Nob:
                                obj=0
            if tt == 1:
                x2=xa[nt[ind]]
                y2=ya[nt[ind]]
                idt=indt[nt[ind]]
                alpha,beta=FPS_conf(x2-xr[i],y2-yr[i],r1=r1,r2=r2)
                xa=np.delete(xa,nt[ind])
                ya=np.delete(ya,nt[ind])
                indt=np.delete(indt,nt[ind])  
                indo.extend([idt])  
                xb.extend([x2])
                yb.extend([y2])
                alp.extend([alpha])
                bet.extend([beta])
                flavt.extend([flav])
            else:
                indo.extend([-1])
                xb.extend([np.nan])
                yb.extend([np.nan])
                alp.extend([180.0])
                bet.extend([180.0])
                flavt.extend([-1])
        else:
            indo.extend([-1])
            xb.extend([np.nan])
            yb.extend([np.nan])
            alp.extend([180.0])
            bet.extend([180.0])
            flavt.extend([-1])
            
    if len(xb) > 0:
        flavt=np.array(flavt)
        indo=np.array(indo)
        xb=np.array(xb)
        yb=np.array(yb)
        alp=np.array(alp)
        bet=np.array(bet)
    return xb,yb,alp,bet,indo,flavt

def FPS_conf2(x,y,alp,bet,r1=7.4,r2=15.0):
    x1=r1*np.cos(alp*np.pi/180.0)+x
    y1=r1*np.sin(alp*np.pi/180.0)+y
    x2=r2*np.cos((alp+bet)*np.pi/180.0)+x1
    y2=r2*np.sin((alp+bet)*np.pi/180.0)+y1
    return x2,y2

def plot_armfps(x,y,alp,bet,r1=7.4,r2=15.0,color='grey',lw=1.0):
    import matplotlib.pyplot as plt
    x1=r1*np.cos(alp*np.pi/180.0)+x
    y1=r1*np.sin(alp*np.pi/180.0)+y
    x2=r2*np.cos((alp+bet)*np.pi/180.0)+x1
    y2=r2*np.sin((alp+bet)*np.pi/180.0)+y1
    for i in range(0, len(x)):
        rax=[x[i],x1[i]]
        ray=[y[i],y1[i]]
        rbx=[x1[i],x2[i]]
        rby=[y1[i],y2[i]]
        plt.plot(rax,ray,'-',color=color,lw=lw)
        plt.plot(rbx,rby,'-',color=color,lw=lw)

def plot_fps_tar(xof,yof,alp,bet,x,y,typ,flav,r1=7.4,r2=15.0,name='Target',dir='./'):
    import matplotlib.pyplot as plt
    nt1=np.where(typ == 'BA')[0]
    nt2=np.where(typ == 'BOSS')[0]
    nt3=np.where(typ == 'Fiducial')[0]
    nt4=np.where(typ != 'Fiducial')[0]
    n1=np.where(flav == 0)[0]
    n2=np.where(flav == 1)[0]
    n3=np.where(flav == 2)[0]
    fig, ax = plt.subplots(figsize=(6.3,5.5))
    plt.xlim(-350,350)
    plt.ylim(-300,300)
    plt.plot(xof[n1],yof[n1],'o',markersize=0.2*10,color='firebrick')
    plt.plot(xof[n2],yof[n2],'o',markersize=0.2*10,color='mediumseagreen')
    plt.plot(xof[n3],yof[n3],'o',markersize=0.2*10,color='royalblue')
    plt.plot(x[nt1],y[nt1],'+',color='blue',markersize=0.3*10)
    plt.plot(x[nt2],y[nt2],'+',color='red',markersize=0.3*10)
    plt.xlabel('X [mm]')
    plt.ylabel('Y [mm]')
    plt.title('Target Configuration')
    plot_don(x[nt1],y[nt1],r1=r1,r2=r2)
    plot_don(x[nt2],y[nt2],r1=r1,r2=r2,color='red')
    plot_don(x[nt3],y[nt3],r1=0.0,r2=5.0,color='green',type=2,apl='1.0')
    plot_armfps(x[nt4],y[nt4],alp,bet,r1=r1,r2=r2,color='grey',lw=1.0)
    fig.tight_layout()
    plt.savefig(dir+name+'.pdf',dpi = 1000)
    plt.close()

def fps_layout(file_n='fps_filledHex.txt',dir="libs/"):
    file=dir+file_n
    f=open(file,'r')
    row=[]
    col=[]
    x=[]
    y=[]
    typ=[]
    for line in f:
        if not '#' in line:
            data=line.replace('\n','').split(' ')
            data=filter(None,data)
            row.extend([np.float(data[0])])
            col.extend([np.float(data[1])])
            x.extend([np.float(data[2])])
            y.extend([np.float(data[3])])
            typ.extend([data[4]])
    f.close
    row=np.array(row)
    col=np.array(col)
    x=np.array(x)
    y=np.array(y)
    typ=np.array(typ)
    return row,col,x,y,typ
            
def check_obs_illus(ra0,dec0):
    file='Illustris_obs.csv'
    if ptt.exists(file) == True:
        ra_o=[]
        dec_o=[]
        f=open(file,'r')
        for line in f:
            if not 'RA' in line:
                data=line.replace('\n','').split(',')
                data=filter(None,data)
                ra_o.extend([np.float(data[0])])
                dec_o.extend([np.float(data[1])])
        ra_o=np.array(ra_o)
        dec_o=np.array(dec_o)
        #print ra_o
        ind=[]
        #r_mint=np.copy(ra0)
        for i in range(0, len(ra0)):
            r=np.sqrt((ra_o-ra0[i])**2.0+(dec_o-dec0[i])**2.0)*3600.
            r_min=np.nanmin(r)
            #r_mint[i]=r_min
            if r_min > 3.0:
                ind.extend([i])
        ind=np.array(ind)
        out=True
        #print np.nanmin(r_mint)
    else:
        ind=range(0,len(ra0))
        out=False
    if len(ind) == 0:
        ind=range(0,len(ra0))
        out=False
    #print ind
    #sys.exit()
    return ind,out

def ssp_extract_sky(template):
    f=open(template,'r')
    sky=[]
    wave_t=[]
    wave=[]
    cont=0
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=filter(None,data)
            if cont == 0:
                wave_t.extend([3300.0])
                sky.extend([np.float(data[1])])
            wave_t.extend([10.0**np.float(data[0])])
            wave.extend([10.0**np.float(data[0])])
            sky.extend([np.float(data[1])])
            cont=cont+1
    wave_t.extend([10400.0])
    sky.extend([np.float(data[1])])
    sky=np.array(sky)
    wave_t=np.array(wave_t)
    wave=np.array(wave)
    n_c=len(wave)
    crpix=1.0
    cdelt=(np.amax(wave)-np.amin(wave))/np.float(n_c)
    n_t=np.int(np.round((np.amax(wave_t)-np.amin(wave_t))/cdelt))
    crval=wave_t[0]
    wave_c=crval+(np.arange(n_t)+1.0-crpix)*cdelt
    pdl_flux_c_ini=interp1d(wave_t,sky,bounds_error=False,fill_value=0.)(wave_c)
    return [pdl_flux_c_ini,wave_c,crval,cdelt,crpix]

def ssp_extract_wd(template,wave_min=2500,wave_max=13000):
    f=open(template,'r')
    flux=[]
    wave=[]
    for line in f:
        if not "#" in line:
            data=line.replace('\n','').split(' ')
            data=filter(None,data)
            wave.extend([np.float(data[0])])
            flux.extend([np.float(data[1])])
    flux=np.array(flux)
    wave=np.array(wave)
    n_c=len(wave)
    crpix=1.0
    cdelt=1.5#(np.amax(wave)-np.amin(wave))/np.float(n_c)
    n_t=np.int(np.round((wave_max-wave_min)/cdelt))
    crval=wave_min#wave[0]
    wave_c=crval+(np.arange(n_t)+1.0-crpix)*cdelt
    pdl_flux_c_ini=interp1d(wave,flux,bounds_error=False,fill_value=0.)(wave_c)
    return [pdl_flux_c_ini,wave_c,crval,cdelt,crpix]

def ssp_extract_star(template,star_type="F0II (25291)"):
    [pdl_flux_c_ini,hdr]=gdata(template, 0, header=True)
    [nf,n_c]=pdl_flux_c_ini.shape    
    cdwa=hdr['COEFF1']
    minw=hdr['COEFF0']
    nwav=hdr['NAXIS1']
    wave=10.0**(np.arange(nwav)*cdwa+minw)
    st=1    
    for i in range(0, nf):
        header="NAME"+str(i)
        name=hdr[header]
        if name in star_type:
            pdl_flux=pdl_flux_c_ini[i,:]
            st=0
#    if st == 1:
#        pdl_flux=pdl_flux_c_ini[49,:]
    crpix=1.0
    cdelt=(np.amax(wave)-np.amin(wave))/np.float(n_c)
    #print cdelt,"H"
    crval=wave[0]
    wave_c=crval+(np.arange(n_c)+1.0-crpix)*cdelt
    pdl_flux_c_ini=interp1d(wave,pdl_flux,bounds_error=False,fill_value=0.)(wave_c)
    return [pdl_flux_c_ini,wave_c,crval,cdelt,crpix]

    
def ssp_extract(template):
    [pdl_flux_c_ini,hdr]=gdata(template, 0, header=True)
    [nf,n_c]=pdl_flux_c_ini.shape
    coeffs=np.zeros([nf,3])
    crpix=hdr["CRPIX1"]
    cdelt=hdr["CDELT1"]
    crval=hdr["CRVAL1"]
    age_mod=[]
    met_mod=[]
    ml=[]
    name=[]
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_ssp_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        AGE=data[0]
        MET=data[1]
        if 'Myr' in AGE:
            age=AGE.replace('Myr','')
            age=float_(age)/1000.
        else:
            age=AGE.replace('Gyr','')
            age=float_(age)
        met=float_(MET.replace('z','0.'))
        age_mod.extend([age])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr[header])
        if val_ml != 0:
            ml.extend([1/val_ml])
        else:
            ml.extend([1])
    wave_c=[]
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c.extend([(crval+cdelt*(j+1-crpix))])
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])
    wave_c=np.array(wave_c)
    ml=np.array(ml)
    age_mod=np.array(age_mod)
    met_mod=np.array(met_mod)

    return [pdl_flux_c_ini,wave_c,age_mod,met_mod,ml,crval,cdelt,crpix]

def simpson_r(f,x,i1,i2,typ=0):
    n=(i2-i1)*1.0
    if n % 2:
        n=n+1.0
        i2=i2+1
    b=x[i2]
    a=x[i1]
    h=(b-a)/n
    s= f[i1]+f[i2]
    n=int(n)
    dx=b-a
    for i in range(1, n, 2):
        s += 4 * f[i1+i]
    for i in range(2, n-1, 2):
        s += 2 * f[i1+i]
    if typ == 0:
        return s*h/3.0
    if typ == 1:
        return s*h/3.0/dx
             
def read_data(dir,file1,file2):
#    file1='all_stars_a1.001.out'
#    file1='all_stars_z0.0_r500kpc_centered.out'
    #file2='gas_component_5Rvir_a1.001.out'
    f1=open(dir+file1,'r')
    f2=open(dir+file2,'r')
    vel_str_x=[]
    vel_str_y=[]
    vel_str_z=[]
    pos_str_x=[]
    pos_str_y=[]
    pos_str_z=[]
    massf=[]
    mass0=[]
    age=[]
    metal1=[]
    metal2=[]
    scalf=[]
    for line in f1:
#    for j in range(0, 100000):
#        line=f1.readline()
        if not "#" in line: 
            line=line.replace('\n','')
            data=line.split(',')
            data=filter(None,data)
            pos_str_x.extend([float_(data[5])])
            pos_str_y.extend([float_(data[6])])
            pos_str_z.extend([float_(data[7])])
            vel_str_x.extend([float_(data[8])])
            vel_str_y.extend([float_(data[9])])
            vel_str_z.extend([float_(data[10])])
            massf.extend([float_(data[4])])
            mass0.extend([float_(data[3])])
            metal1.extend([float_(data[12])])
            metal2.extend([float_(data[13])])        
            age.extend([float_(data[1])])
            scalf.extend([float_(data[2])])
    f1.close()
    pos_str_x=np.array(pos_str_x)*1e3
    pos_str_y=np.array(pos_str_y)*1e3
    pos_str_z=np.array(pos_str_z)*1e3
    vel_str_x=np.array(vel_str_x)
    vel_str_y=np.array(vel_str_y)
    vel_str_z=np.array(vel_str_z)
    massf=np.array(massf)
    mass0=np.array(mass0)
    metal1=np.array(metal1)
    metal2=np.array(metal2)    
    #metal=1./((10.**metal*(1.0/0.0127-1.0))+1.)
    metal=(metal1+metal2)
    age=np.array(age)
    scalf=1.0/np.array(scalf)-1
    pos_str=np.zeros([len(pos_str_x),3])
    pos_str[:,0]=pos_str_x
    pos_str[:,1]=pos_str_y
    pos_str[:,2]=pos_str_z
    vel_str=np.zeros([len(pos_str_x),3])
    vel_str[:,0]=vel_str_x
    vel_str[:,1]=vel_str_y
    vel_str[:,2]=vel_str_z
    vel_gas_x=[]
    vel_gas_y=[]
    vel_gas_z=[]
    pos_gas_x=[]
    pos_gas_y=[]
    pos_gas_z=[]
    masses_g=[]
    dens_g=[]
    metal_g=[]
    sfr_g=[]
    Av_g=[]
    T_g=[]
    R_g=[]
    fac=[]
    for line in f2:
#    for j in range(0, 100000):
#        line=f2.readline()
        if not "*" in line:
            line=line.replace('\n','')
            data=line.split(',')
            data=filter(None,data)
            pos_gas_x.extend([float_(data[0])])
            pos_gas_y.extend([float_(data[1])])
            pos_gas_z.extend([float_(data[2])])
            vel_gas_x.extend([float_(data[3])])
            vel_gas_y.extend([float_(data[4])])
            vel_gas_z.extend([float_(data[5])])
            masses_g.extend([float_(data[9])])
            metal_g.extend([float_(data[10])])        
            dens_g.extend([float_(data[7])])
            sfr_g.extend([float_(data[11])])
            Av_g.extend([float_(data[6])])
            T_g.extend([float_(data[8])])
            R_g.extend([float_(data[12])])
            if float_(data[10]) > (10.0**(-0.59)*0.017):
                fac.extend([10.0**(2.21-1.0)])
            else:
                fac.extend([(10.0**(2.21-1.0))/(float_(data[10])/0.017)**(3.1-1.0)])
    f2.close()
    fac=np.array(fac)
    pos_gas_x=np.array(pos_gas_x)*1e3
    pos_gas_y=np.array(pos_gas_y)*1e3
    pos_gas_z=np.array(pos_gas_z)*1e3
    vel_gas_x=np.array(vel_gas_x)
    vel_gas_y=np.array(vel_gas_y)
    vel_gas_z=np.array(vel_gas_z)
    masses_g=np.array(masses_g)
    metal_g=np.array(metal_g)
    dens_g=np.array(dens_g)
    sfr_g=np.array(sfr_g)
    Av_g=np.array(Av_g)/0.017/fac/0.017#*0.5/0.017
    R_g=np.array(R_g)
    T_g=np.array(T_g)
    pos_gas=np.zeros([len(pos_gas_x),3])
    pos_gas[:,0]=pos_gas_x
    pos_gas[:,1]=pos_gas_y
    pos_gas[:,2]=pos_gas_z
    vel_gas=np.zeros([len(pos_gas_x),3])
    vel_gas[:,0]=vel_gas_x
    vel_gas[:,1]=vel_gas_y
    vel_gas[:,2]=vel_gas_z
    #f2.close()
    gas={'Coordinates':pos_gas,'Velocities':vel_gas,'Density':dens_g,'GFM_Metallicity':metal_g,'Radius':R_g,'StarFormationRate':sfr_g,'Extintion':Av_g,'Temperature':T_g,'GFM_Masses':masses_g}
    stars={'Coordinates':pos_str,'Velocities':vel_str,'Masses':massf,'GFM_InitialMass':mass0,'GFM_StellarFormationTime':scalf,'GFM_StellarAge':age,'GFM_Metallicity':metal}
    return [gas, stars]

def mock_halo_test(idh,basename='artsp8-',dir1='',fib_n=7,psf=0,nl=110,fov=30.0,thet=0.0,ifutype="MaNGA"):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
    cubef=basename+id+'.cube_test'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    dir1n=dir1.replace(" ","\ ")
    cube_conv_t(cubef,psfi=psf,dir_o=dir1,nl=fib_n,fov=fov,thet=thet,ifutype=ifutype)
    sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
                 
def mock_halo(idh,sp_res=0.0,sp_samp=1.25,basename='artsp8-',template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",file=file,file2='',dir1='',fib_n=7,psf=0,ho=0.704,Lam=0.7274,SN=15.0,Fluxm=20.0,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="MaNGA"):
    #if not "CALIFA" in ifutype or not "MaNGA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
#    cubef='cube_'+id
    cubef=basename+id+'.cube'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    [gas,stars]=read_data(dir0,file,file2)
    dx_g = gas['Coordinates'][:,0]
    dy_g = gas['Coordinates'][:,1]
    dz_g = gas['Coordinates'][:,2]
    vx_g = gas['Velocities'][:,0]
    vy_g = gas['Velocities'][:,1]
    vz_g = gas['Velocities'][:,2]
    ra_g = gas['Radius'][:]
    dens = gas['Density'][:]
    sfri = gas['StarFormationRate'][:]
    meta_g=gas['GFM_Metallicity'][:]
    mass_g=gas['GFM_Masses'][:]
    temp_g=gas['Temperature'][:]
    Av_g = gas['Extintion'][:]        
    dx = stars['Coordinates'][:,0]
    dy = stars['Coordinates'][:,1]
    dz = stars['Coordinates'][:,2]
    mass=stars['Masses'][:]
    mas0=stars['GFM_InitialMass'][:]
    meta=stars['GFM_Metallicity'][:]
    time=stars['GFM_StellarFormationTime'][:]
    vx = stars['Velocities'][:,0]
    vy = stars['Velocities'][:,1]
    vz = stars['Velocities'][:,2]
    ages=stars['GFM_StellarAge'][:]
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    mass_g=mass_g#/ho
    mass=mass[nt]#/ho
    mas0=mas0[nt]#/ho
    meta=meta[nt]#/0.0127
    ages=ages[nt]
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=ra_g/ho#**3.0
    dens=dens*ho**2.0
    zf=time[nt]
    xo=0#modes(dx,nit=7,n_s=0.8)
    yo=0#modes(dy,nit=7,n_s=0.8)
    zo=0#modes(dz,nit=7,n_s=0.8)
#    print xo,yo,zo,"HOLA",len(dx)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
#    print modes(x,nit=7,n_s=0.8),modes(y,nit=7,n_s=0.8),modes(z,nit=7,n_s=0.8)
#    sys.exit()
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    dir1n=dir1.replace(" ","\ ")
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        cube_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,sp_res=sp_res,sp_samp=sp_samp,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=fib_n,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
        #sys.exit()
        band_cube(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
    else:
        band_cube(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
    sycall('mv *'+basename+id+'* '+dir1n)
        

def mock_halo_s(idh,basename='artsp8-',template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",file=file,file2='',dir1='',psf=0,ho=0.704,Lam=0.7274,SN=15.0,Fluxm=20.0,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="SDSS"):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
#    cubef='cube_'+id
    cubef=basename+id+'.spec'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    [gas,stars]=read_data(dir0,file,file2)
    dx_g = gas['Coordinates'][:,0]
    dy_g = gas['Coordinates'][:,1]
    dz_g = gas['Coordinates'][:,2]
    vx_g = gas['Velocities'][:,0]
    vy_g = gas['Velocities'][:,1]
    vz_g = gas['Velocities'][:,2]
    ra_g = gas['Radius'][:]
    dens = gas['Density'][:]
    sfri = gas['StarFormationRate'][:]
    meta_g=gas['GFM_Metallicity'][:]
    mass_g=gas['GFM_Masses'][:]
    temp_g=gas['Temperature'][:]
    Av_g = gas['Extintion'][:]        
    dx = stars['Coordinates'][:,0]
    dy = stars['Coordinates'][:,1]
    dz = stars['Coordinates'][:,2]
    mass=stars['Masses'][:]
    mas0=stars['GFM_InitialMass'][:]
    meta=stars['GFM_Metallicity'][:]
    time=stars['GFM_StellarFormationTime'][:]
    vx = stars['Velocities'][:,0]
    vy = stars['Velocities'][:,1]
    vz = stars['Velocities'][:,2]
    ages=stars['GFM_StellarAge'][:]
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    mass_g=mass_g#/ho
    mass=mass[nt]#/ho
    mas0=mas0[nt]#/ho
    meta=meta[nt]#/0.0127
    ages=ages[nt]
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=ra_g/ho#**3.0
    dens=dens*ho**2.0
    zf=time[nt]
    xo=0#modes(dx,nit=7,n_s=0.8)
    yo=0#modes(dy,nit=7,n_s=0.8)
    zo=0#modes(dz,nit=7,n_s=0.8)
#    print xo,yo,zo,"HOLA",len(dx)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    dir1n=dir1.replace(" ","\ ")
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        fib_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)

def mock_photo(id,basename='artsp8-',template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photo'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
        #[age_F,ind]=ages_definition(age_s,n_ages=55)
        photo_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
     #   band_photo(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
    #    [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".r_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
    #    sycall("mv "+cubef+".* "+dir1n)
        if ptt.exists(dir1+cubef+".rg_Lc_rad.fits") == True:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=0)
        else:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=[-100,-100,-100,-100,-100,-100]
        hdulist = fits.open(dir1+cubef+".rg_Lc.fits")
        hd=hdulist[0].header
        red_0=hd["REDSHIFT"]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    #sys.exit()
                
    if np.abs(not_fit) > 0:
        R50=-100
    #f4.write(cubef+','+str(R50)+','+str(e_rad)+','+str(red_0)+','+str(ab)+','+str(PA)+','+str(not_fit)+'\n')
    #return [PA,np.sqrt(1.0-ab**2.),e_rad,not_fit]

        
def mock_sim(id,basename='artsp8-',template2="../../Base_bc03/templete_bc03_2.fits",template="../home/sanchez/ppak/legacy/gsd61_156.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,1.5],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='simu-'+id
    cubef=basename+id+'.simu'
    dirf=dir1t+'photo_lin/'
    #dirf2=dir1t+'Ensamble/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
       # nt=np.where(np.abs(z) <= 2.0)[0]
       # x=x[nt]
        #y=y[nt]
        #z=z[nt]
        #vx=vx[nt]
        #vy=vy[nt]
        #vz=vz[nt]
        #ages=ages[nt]
        #meta=meta[nt]
        #mass=mass[nt]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
        #[age_F,ind]=ages_definition(age_s,n_ages=55)
        sim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
        #    ensamble_s(cubef+'.fits.gz',dir1=dir1)
        #    ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            #sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            #sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        #ensamble_int(cubef+'.fits.gz',dir1=dir1)
        #ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        #sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
        #    ensamble_s(cubef+'.fits.gz',dir1=dir1)
        #    ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            #sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            #sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        #ensamble_int(cubef+'.fits.gz',dir1=dir1)
        #ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        #sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
        
def mock_simvel(id,basename='artsp8-',dir1='',file='',file2='',ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,thet=0.0,observer=[0,0,0],observer2=[0,0,0],ang=0.0):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    #outf='image_'+id
    #cubef='simu-'+id
    cubef=basename+id+'.simuvel'
    dirf=dir1t+'photo_lin/'
    #dirf2=dir1t+'Ensamble/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    a=1
    if a==1:#ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)

        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        
        x2=np.copy(x)
        y2=np.copy(y)
        z2=np.copy(z)
        x_g2=np.copy(x_g)
        y_g2=np.copy(y_g)
        z_g2=np.copy(z_g)
        vx2=np.copy(vx)
        vy2=np.copy(vy)
        vz2=np.copy(vz)
        vx_g2=np.copy(vx_g)
        vy_g2=np.copy(vy_g)
        vz_g2=np.copy(vz_g)
        
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        nt=np.where(np.abs(z) <= 4.0)[0]#Particles withi the disk of 2 kpc
        x=x[nt]
        y=y[nt]
        z=z[nt]
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        nt_g=np.where(np.abs(z_g) <= 4.0)[0]#Particles withi the disk of 2 kpc
        x_g=x_g[nt_g]
        y_g=y_g[nt_g]
        z_g=z_g[nt_g]
        vx_g=vx_g[nt_g]
        vy_g=vy_g[nt_g]
        vz_g=vz_g[nt_g]
        
        vxo=modes(vx,nit=7,n_s=0.8)
        vyo=modes(vy,nit=7,n_s=0.8)
        vzo=modes(vz,nit=7,n_s=0.8)
        vx=vx-vxo
        vy=vy-vyo
        vz=vz-vzo
        
        vx_g=vx_g-vxo
        vy_g=vy_g-vyo
        vz_g=vz_g-vzo
        
        x02=observer2[0]-xo
        y02=observer2[1]-yo
        z02=observer2[2]-zo
        Rc2=np.sqrt(x02**2.+y02**2.+z02**2.)
        red_02=reds_cos(Rc2/1e3)
        Ra2=Rc2#/(1+red_0)
        A12=np.arctan2(y02,z02)
        A22=np.arcsin(x02/Ra2)
        R12=np.array([[1,0,0],[0,np.cos(A12),-np.sin(A12)],[0,np.sin(A12),np.cos(A12)]])
        R22=np.array([[np.cos(A22),0,-np.sin(A22)],[0,1,0],[np.sin(A22),0,np.cos(A22)]])
        R32=np.array([[np.cos(A22),-np.sin(A22),0],[np.sin(A22),np.cos(A22),0],[0,0,1]])
        Ve2=np.array([x2,y2,z2])
        Vf2=np.dot(np.dot(R22,R12),Ve2)
        x2=Vf2[0]
        y2=Vf2[1]
        z2=Vf2[2]
        Ve2=np.array([vx2,vy2,vz2])
        Vf2=np.dot(np.dot(R22,R12),Ve2)
        vx2=Vf2[0]
        vy2=Vf2[1]
        vz2=Vf2[2]
        Ve2=np.array([x_g2,y_g2,z_g2])
        Vf2=np.dot(np.dot(R22,R12),Ve2)
        x_g2=Vf2[0]
        y_g2=Vf2[1]
        z_g2=Vf2[2]
        Ve2=np.array([vx_g2,vy_g2,vz_g2])
        Vf2=np.dot(np.dot(R22,R12),Ve2)
        vx_g2=Vf2[0]
        vy_g2=Vf2[1]
        vz_g2=Vf2[2]
        #nt=np.where(np.abs(z) <= 2.0)[0]#Particles withi the disk of 2 kpc
        x2=x2[nt]
        y2=y2[nt]
        z2=z2[nt]
        vx2=vx2[nt]
        vy2=vy2[nt]
        vz2=vz2[nt]
        x_g2=x_g2[nt_g]
        y_g2=y_g2[nt_g]
        z_g2=z_g2[nt_g]
        vx_g2=vx_g2[nt_g]
        vy_g2=vy_g2[nt_g]
        vz_g2=vz_g2[nt_g]

        
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        sim_conv_vel(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,x2,y2,z2,vx2,vy2,vz2,x_g2,y_g2,z_g2,vx_g2,vy_g2,vz_g2,dir_o=dir1,red_0=red_0,red_02=red_02,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,thet=thet,observer=observer,observer2=observer2,ang=ang)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        
        
def mock_photosim(id,basename='artsp8-',template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photosim'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3
        photosim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        #band_photo(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
        
        
def mock_photosimextgas(id,basename='artsp8-',template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photosimextgas'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3
        photosimextgas_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
        sycall("mv "+cubef+".* "+dir1n)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)        
        if ptt.exists(dir1+cubef+".rg_Lc_rad.fits") == True:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=0)
        else:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
        sycall("mv "+cubef+".* "+dir1n)
            
def ensamble_s(name, dir='', dir1='',rx=[0.0,0.5,1.0,1.5],lig=0):
    file2=dir1+name.replace('simu','photo').replace('_L','').replace(".fits.gz",".rg_Lc_rad.fits")
#    if ptt.exists(file2) == True:
    [pdl_rad, hdr2]=gdata(file2,0, header=True)
    ind=[]
    inda=[]
    nr=len(rx)
    rad=1.0
    for ii in range(0, nr-1):
        nt=np.where(pdl_rad< rx[ii+1]*rad)
        nta=np.where(pdl_rad[nt]>= rx[ii]*rad)
        ind.extend([nt])
        inda.extend([nta])
    vel_light=299792.458
    [mass_age_t,hdr]=gdata(dir1+name, 0, header=True)
    [n_age,nx,ny]=mass_age_t.shape
    ages=np.zeros(n_age)
    for i in range(0, n_age):
        ages[i]=np.log10(hdr["AGE"+str(i)])+9.0
    mass_age_t_2=mass_age_t
    mass_age_t_e=np.zeros([n_age,nx,ny])
    MassT=0
    massN=np.zeros(nr)
    massN_e=np.zeros(nr)
    mass=np.zeros([n_age,nr])
    mass_e=np.zeros([n_age,nr])
    sfrt=np.zeros([n_age,nr])
    sfdt=np.zeros([n_age])
    mass_temp_total=0
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
        if i == 0:
            temp1=temp
        else:
            temp1=temp1+temp
        MassT=MassT+np.sum(temp)
        for ii in range(0, nr-1):
            Dt_mass=np.sum(temp[ind[ii]][inda[ii]])
            Dt_mass_e=np.sum(temp_e[ind[ii]][inda[ii]])
            Dt_mass_2=np.sum(temp_2[ind[ii]][inda[ii]])
            massN[ii]=Dt_mass+massN[ii]
            massN_e[ii]=Dt_mass_e+massN_e[ii]
            mass[i,ii]=np.log10(massN[ii])
            mass_e[i,ii]=massN_e[ii]
            sfrt[i,ii]=Dt_mass_2/Dt_age
#    print mass[n_age-1,:]
        #print massN
        #print nr
#    sys.exit()
    mass_temp_total=np.log10(np.sum(10**mass[n_age-1,:]))
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    mass_n=10**(10**(mass-mass[n_age-1,:]))
    mass_n_e=np.sqrt((10**(mass-mass[n_age-1,:]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[n_age-1,:]/10**(2.0*mass[n_age-1,:]))))
    if lig == 1:
        f2=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_L.csv","w")
        f2.write("#  LOG_AGE  N_LIGHR_1  N_LIGHR_2  N_LIGHR_3  N_LIGHR_4  LOG_LIGHR_1  LOG_LIGHR_2  LOG_LIGHR_3  LOG_LIGHR_4 \n")
    else:
        f2=open(dir1+name.replace(".fits.gz","")+"_Ensemble.csv","w")
        f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    for i in range(0, n_age):
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
    f2.close()
    
    
    nrb=1+4*3+1
    ensamb=np.zeros([n_age,nrb])
    co=0
    if lig == 1:
        f=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_L.csv","r")
    else:
        f=open(dir1+name.replace(".fits.gz","")+"_Ensemble.csv","r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1

    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1:4])
    sfh=ensamb[:,7:10]#/10.0**ensamb[:,4:7]
    emgh=(ensamb[:,10:13])#*10
    #print ensamb[:,7:10]
    #print 10.0**ensamb[:,4:7]
    dth=ensamb[:,13]/1e3
    Ddth=np.zeros(len(dth))
    [nx,ny]=mgh.shape
    Dmgh=np.zeros([nx,ny])
    exy1=np.zeros([nx,ny])
    exy2=np.zeros([nx,ny])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i,:]-mgh[i+1,:])
            Ddth[i]=np.abs(dth[i]-dth[i+1])#/2.0
        elif i == len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i-1,:]-mgh[i,:])
            Ddth[i]=np.abs(dth[i-1]-dth[i])#/2.0 
    exy1[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)#np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]
    exy1[:,1]=np.amin([np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1],np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]],axis=0)#np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]
    exy1[:,2]=np.amin([np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2],np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]],axis=0)#np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]
    exy2[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)#np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0]
    exy2[:,1]=np.amin([np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1],np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]],axis=0)#np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1]
    exy2[:,2]=np.amin([np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2],np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]],axis=0)#np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2]

    for i in range(0, 3):
        val1=np.abs(Dmgh[:,i]/(dth+Ddth)*dth)+emgh[:,i]
        val2=np.abs(Dmgh[:,i]/(dth-Ddth)*dth)+emgh[:,i]
        for j in range(0, len(val1)):
            if val1[j] >= 0.15:
                val1[j]=0.1
            if val2[j] >= 0.15:
                val2[j]=0.1
        exy1[:,i]=val1
        exy2[:,i]=val2
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title("SIMULATION LGH INPUT")
    else:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title("SIMULATION MGH INPUT")
    plt.semilogx(age,mgh[:,0],'',color="b",label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=1.5)
    plt.semilogx(age,mgh[:,1],'--',color="g",label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=1.5)
    plt.semilogx(age,mgh[:,2],':',color="r",label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=1.5)
    ax.fill_between(age,mgh[:,0]+exy1[:,0],mgh[:,0]-exy2[:,0],alpha='0.20',color="b")
    ax.fill_between(age,mgh[:,1]+exy1[:,1],mgh[:,1]-exy2[:,1],alpha='0.20',color="g")
    ax.fill_between(age,mgh[:,2]+exy1[:,2],mgh[:,2]-exy2[:,2],alpha='0.20',color="r")

    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'LGH_'+name.replace('_L','')+'.pdf',dpi = 1000)        
    else:
        plt.savefig(dir1+'MGH_'+name+'.pdf',dpi = 1000)
    plt.close()
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)

    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$SFHL(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFHL INPUT")
    else:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFH INPUT")
    plt.semilogx(age,sfh[:,0],'',color="b",drawstyle='steps',label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=2)
    plt.semilogx(age,sfh[:,1],'--',color="g",drawstyle='steps',label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=2)
    plt.semilogx(age,sfh[:,2],':',color="r",drawstyle='steps',label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'SFHL_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'SFH_'+name+'.pdf',dpi = 1000)
    plt.close()
    
def ensamble_int(name, dir='', dir1='',lig=0):
#    file2=dir1+name.replace('simu','photo').replace(".fits.gz",".r_Lc_rad.fits")
#    if ptt.exists(file2) == True:
#    [pdl_rad, hdr2]=gdata(file2,0, header=True)
    ind=[]
    inda=[]
    rad=12.0
    vel_light=299792.458
    [mass_age_t,hdr]=gdata(dir1+name, 0, header=True)
    dpix=hdr['CD1_1']*3600.0
    [n_age,nx,ny]=mass_age_t.shape
    pdl_rad=np.zeros([nx,ny])
    for i in range(0, nx):
        for j in range(0, ny):
            pdl_rad[i,j]=np.sqrt((i-nx/2.)**2.0+(j-ny/2.0)**2.0)*dpix/(2.5*0.5)
    nt=np.where(pdl_rad< rad)
    ages=np.zeros(n_age)
    for i in range(0, n_age):
        ages[i]=np.log10(hdr["AGE"+str(i)])+9.0
    mass_age_t_2=mass_age_t
    mass_age_t_e=np.zeros([n_age,nx,ny])
    MassT=0
    massN=0
    massN_e=0
    mass=np.zeros([n_age])
    mass_e=np.zeros([n_age])
    sfrt=np.zeros([n_age])
    sfdt=np.zeros([n_age])
    mass_temp_total=0
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
        if i == 0:
            temp1=temp
        else:
            temp1=temp1+temp
        MassT=MassT+np.sum(temp)
        
        Dt_mass=np.sum(temp[nt])
        Dt_mass_e=np.sum(temp_e[nt])
        Dt_mass_2=np.sum(temp_2[nt])
        massN=Dt_mass+massN
        massN_e=Dt_mass_e+massN_e
        mass[i]=np.log10(massN)
        mass_e[i]=massN_e
        sfrt[i]=Dt_mass_2/Dt_age
#    print mass[n_age-1,:]
        #print massN
        #print nr
#    sys.exit()
    mass_temp_total=np.log10(10**mass[n_age-1])
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    mass_n=10**(10**(mass-mass[n_age-1]))
    mass_n_e=np.sqrt((10**(mass-mass[n_age-1]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[n_age-1]/10**(2.0*mass[n_age-1]))))
    if lig == 1:
        f2=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_int_L.csv","w")
        f2.write("#  LOG_AGE  N_LIGHR_1  N_LIGHR_2  N_LIGHR_3  N_LIGHR_4  LOG_LIGHR_1  LOG_LIGHR_2  LOG_LIGHR_3  LOG_LIGHR_4 \n")
    else:
        f2=open(dir1+name.replace(".fits.gz","")+"_Ensemble_int.csv","w")
        f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    for i in range(0, n_age):
        line=''
        line=line+str(ages[i])
        line=line+';'+str(mass_n[i])
        line=line+';'+str(mass[i])
        line=line+';'+str(sfrt[i])
        line=line+';'+str(mass_n_e[i])
        line=line+';'+str(sfdt[i])
        line=line+' \n'
        f2.write(line)
    f2.close()
    
    
    nrb=6
    ensamb=np.zeros([n_age,nrb])
    co=0
    if lig == 1:
        f=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_int_L.csv","r")
    else:
        f=open(dir1+name.replace(".fits.gz","")+"_Ensemble_int.csv","r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1

    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1])
    sfh=ensamb[:,3]#/10.0**ensamb[:,4:7]
    emgh=(ensamb[:,4])#*10
    #print ensamb[:,7:10]
    #print 10.0**ensamb[:,4:7]
    dth=ensamb[:,5]/1e3
    Ddth=np.zeros(len(dth))
    nx=len(mgh)
    Dmgh=np.zeros([nx])
    exy1=np.zeros([nx])
    exy2=np.zeros([nx])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i]-mgh[i+1])
            Ddth[i]=np.abs(dth[i]-dth[i+1])#/2.0
        elif i == len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i-1]-mgh[i])
            Ddth[i]=np.abs(dth[i-1]-dth[i])#/2.0 
    exy1=np.amin([np.abs(Dmgh/(dth-Ddth)*dth)+emgh,np.abs(Dmgh/(dth+Ddth)*dth)+emgh],axis=0)#np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]
    exy2=np.amin([np.abs(Dmgh/(dth-Ddth)*dth)+emgh,np.abs(Dmgh/(dth+Ddth)*dth)+emgh],axis=0)#np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0]

    val1=np.abs(Dmgh/(dth+Ddth)*dth)+emgh
    val2=np.abs(Dmgh/(dth-Ddth)*dth)+emgh
    for j in range(0, len(val1)):
        if val1[j] >= 0.15:
            val1[j]=0.1
        if val2[j] >= 0.15:
            val2[j]=0.1
    exy1=val1
    exy2=val2
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    #plt.title(name)
    if lig == 1:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title("SIMULATION LGH INT INPUT")
    else:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title("SIMULATION MGH INT INPUT")
    plt.semilogx(age,mgh,'',color="b",label='$R=15arcsec$',lw=1.5)
    ax.fill_between(age,mgh+exy1,mgh-exy2,alpha='0.20',color="b")
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'LGH_int_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'MGH_int_'+name+'.pdf',dpi = 1000)
    plt.close()
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)

    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$SFHL(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFHL INT INPUT")
    else:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFH INT INPUT")
    plt.semilogx(age,sfh,'',color="b",drawstyle='steps',label='$R=15arcsec$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'SFHL_int_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'SFH_int_'+name+'.pdf',dpi = 1000)
    plt.close()
    
def band_photo_r(pdl_flux, hdr, name='test.fits.gz', dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt']#,'U_Johnson.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt','J_2MASS.txt','H_2MASS.txt','K_2MASS.txt']
    band=['u','g','rg','ig','z']#,'U','B','V','I','R','J','H','K']
    zerop=[3631.0,3730.0,4490.0,4760.0,4810.0]#,1810.0,4260.0,3640.0,2550.0,3080.0,1600.0,1080.0,670.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                flux_t=spec1*trans
                #photo_a[i,j]=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2)/simpson_r(trans,wave,0,len(wave)-2)
                photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                photo_b[i,j]=f_fin*zerop[k]*jans
                photo_c[i,j]=f_fi2
#                if photo_a[i,j] > 0:
#                    print photo_a[i,j]
#                    import matplotlib.pyplot as pl
#                    pl.plot(wave,trans)
#                    pl.show()
#                    pl.plot(wave,flux_t)
#                    pl.show() 
#                    pl.plot(wave,spec1)
#                    pl.plot(wave_s,spec)
#                    pl.show()
#                    print simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)
#                    print simpson_r(trans,wave,0,len(wave)-2,typ=1)
#                    print simpson_r((trans/(wave**2.0)*vel_light/ang),wave,0,len(wave)-2)
#                    print simpson_r(trans/((wave*ang)**2.0)*vel_light,wave*ang,0,len(wave)-2)
#                    print simpson_r(trans,vel_light/(wave*ang),0,len(wave)-2)
#                    sys.exit()
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')
        imag_F2=photo_b#convolve(photo_b, PSF)
        imag_F3=photo_c#convolve(photo_c, PSF)
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr2)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    
def band_photo(name, dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [pdl_flux,hdr]=gdata(dir1+name, 0, header=True)
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','NUV_GALEX.txt','FUV_GALEX.txt','U_Johnson.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt','J_2MASS.txt','H_2MASS.txt','K_2MASS.txt']
    band=['ha','u','g','rg','ig','z','NUV','FUV','U','B','V','I','R','J','H','K']
#    zerop=[3631.0,3631.0,3730.0,4490.0,4760.0,4810.0,3801.4,3619.9,1810.0,4260.0,3640.0,2550.0,3080.0,1600.0,1080.0,670.0]
    zerop=[3631.0,3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                flux_t=spec1*trans*d_wave
                f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                photo_b[i,j]=f_fin*zerop[k]*jans
                photo_c[i,j]=f_fi2
#                if photo_a[i,j] > 0:
#                    print photo_a[i,j]
#                    import matplotlib.pyplot as pl
#                    pl.plot(wave,trans)
#                    pl.show()
#                    pl.plot(wave,flux_t)
#                    pl.show() 
#                    pl.plot(wave,spec1)
#                    pl.plot(wave_s,spec)
#                    pl.show()
#                    print simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)
#                    print simpson_r(trans,wave,0,len(wave)-2,typ=1)
#                    print simpson_r((trans/(wave**2.0)*vel_light/ang),wave,0,len(wave)-2)
#                    print simpson_r(trans/((wave*ang)**2.0)*vel_light,wave*ang,0,len(wave)-2)
#                    print simpson_r(trans,vel_light/(wave*ang),0,len(wave)-2)
#                    sys.exit()
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
            del hdr["CUNIT3"]
            del hdr['RADECSYS']
            hdr['RADECSYSa']='ICRS    '
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2/A'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2=photo_b#convolve(photo_b, PSF)#photo_b#
        imag_F3=photo_c#convolve(photo_c, PSF)#photo_c#
        imag_F1C=convolve(photo_a, PSF,  boundary='symm')#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2C=convolve(photo_b, PSF)#photo_b#
        imag_F3C=convolve(photo_c, PSF)#photo_c#
        imag_F1C[np.where(imag_F1C == 0)]=35.0
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr2)
        wfits(dir1+name_f+band[k]+'_c.fits',imag_F1C,hdr)
        wfits(dir1+name_f+band[k]+'_Fc.fits',imag_F2C,hdr2)
        wfits(dir1+name_f+band[k]+'_Lc.fits',imag_F3C,hdr2)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#5e-15)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)#7e-23)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    from PIL import Image
    max=20.0
    min=35.0
    file00u=dir1+name_f+"u_c.fits"
    file00g=dir1+name_f+"g_c.fits"
    file00r=dir1+name_f+"rg_c.fits"
    file00i=dir1+name_f+"ig_c.fits"
    file00z=dir1+name_f+"z_c.fits"
    file00U=dir1+name_f+"U_c.fits"
    file00B=dir1+name_f+"B_c.fits"
    file00V=dir1+name_f+"V_c.fits"
    file00R=dir1+name_f+"R_c.fits"
    file00I=dir1+name_f+"I_c.fits"    
    [pdl_00u,hdrt]=gdata(file00u, 0, header=True)
    [pdl_00g,hdrt]=gdata(file00g, 0, header=True)
    [pdl_00r,hdrt]=gdata(file00r, 0, header=True)
    [pdl_00i,hdrt]=gdata(file00i, 0, header=True)
    [pdl_00z,hdrt]=gdata(file00z, 0, header=True)
    [pdl_00U,hdrt]=gdata(file00U, 0, header=True)
    [pdl_00B,hdrt]=gdata(file00B, 0, header=True)
    [pdl_00V,hdrt]=gdata(file00V, 0, header=True)
    [pdl_00R,hdrt]=gdata(file00R, 0, header=True)
    [pdl_00I,hdrt]=gdata(file00I, 0, header=True)
    nx,ny=pdl_00g.shape
    pdl_00u=(np.flipud(pdl_00u)-min)/(max-min)*256
    pdl_00g=(np.flipud(pdl_00g)-min)/(max-min)*256
    pdl_00r=(np.flipud(pdl_00r)-min)/(max-min)*256
    pdl_00i=(np.flipud(pdl_00i)-min)/(max-min)*256
    pdl_00z=(np.flipud(pdl_00z)-min)/(max-min)*256
    pdl_00U=(np.flipud(pdl_00U)-min)/(max-min)*256
    pdl_00B=(np.flipud(pdl_00B)-min)/(max-min)*256
    pdl_00V=(np.flipud(pdl_00V)-min)/(max-min)*256
    pdl_00R=(np.flipud(pdl_00R)-min)/(max-min)*256
    pdl_00I=(np.flipud(pdl_00I)-min)/(max-min)*256
    pdl_00u[np.where(pdl_00u < 0)]=0
    pdl_00g[np.where(pdl_00g < 0)]=0
    pdl_00r[np.where(pdl_00r < 0)]=0
    pdl_00i[np.where(pdl_00i < 0)]=0
    pdl_00z[np.where(pdl_00z < 0)]=0
    pdl_00U[np.where(pdl_00U < 0)]=0
    pdl_00B[np.where(pdl_00B < 0)]=0
    pdl_00V[np.where(pdl_00V < 0)]=0
    pdl_00R[np.where(pdl_00R < 0)]=0
    pdl_00I[np.where(pdl_00I < 0)]=0
    pdl_00u[np.where(pdl_00u > 255)]=255
    pdl_00g[np.where(pdl_00g > 255)]=255
    pdl_00r[np.where(pdl_00r > 255)]=255
    pdl_00i[np.where(pdl_00i > 255)]=255
    pdl_00z[np.where(pdl_00z > 255)]=255
    pdl_00U[np.where(pdl_00U > 255)]=255
    pdl_00B[np.where(pdl_00B > 255)]=255
    pdl_00V[np.where(pdl_00V > 255)]=255
    pdl_00R[np.where(pdl_00R > 255)]=255
    pdl_00I[np.where(pdl_00I > 255)]=255
    pdl_00=np.zeros([nx,ny,3],dtype="uint8")
    pdl_00[:,:,0]=pdl_00i
    pdl_00[:,:,1]=pdl_00r
    pdl_00[:,:,2]=pdl_00g
    im = Image.fromarray(pdl_00)
    im.save(dir1+name_f+"gri.jpeg")
    pdl_00[:,:,0]=pdl_00R
    pdl_00[:,:,1]=pdl_00V
    pdl_00[:,:,2]=pdl_00B
    im1 = Image.fromarray(pdl_00)
    im1.save(dir1+name_f+"BVR.jpeg")

    
def band_cube(name, dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [pdl_flux,hdr]=gdata(dir1+name, 0, header=True)
    #pdl_noise=gdata(dir1+name, 1)
    pdl_flux=pdl_flux*1e-16
    #pdl_noise=pdl_noise*1e-16
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt']
    band=['ha','u','g','rg','ig','z','B','V','I','R']
#    zerop=[3631.0,3631.0,3730.0,4490.0,4760.0,4810.0,4260.0,3640.0,2550.0,3080.0]
    zerop=[3631.0,3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]#3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        photo_ae=np.zeros([nx,ny])
        photo_be=np.zeros([nx,ny])
        photo_ce=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                #nois=pdl_noise[:,i,j]
                if np.sum(spec) > 0:
                    spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                    #nois1=interp1d(wave_s, nois,bounds_error=False,fill_value=0.)(wave)
                    flux_t=spec1*trans*d_wave
                    #nois_t=nois1*trans
                    f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                    f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                    #f_fin_e=simpson_r(nois_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                    #f_fi2_e=simpson_r(nois_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                    #photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                    if f_fin <= 0:
                        photo_a[i,j]=-2.5*np.log10(1e-10)#14
                    else:
                        photo_a[i,j]=-2.5*np.log10(f_fin+1e-10)
                    photo_b[i,j]=f_fin*zerop[k]*jans
                    photo_c[i,j]=f_fi2
                    #photo_ae[i,j]=-2.5*np.log10(f_fin_e+1e-12)
                    #photo_be[i,j]=f_fin_e*zerop[k]*jans
                    #photo_ce[i,j]=f_fi2_e
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2/A'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2=photo_b#convolve(photo_b, PSF)#photo_b#
        imag_F3=photo_c#convolve(photo_c, PSF)#photo_c#
        #imag_E1=photo_ae#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        #imag_E2=photo_be#convolve(photo_b, PSF)#photo_b#
        #imag_E3=photo_ce#convolve(photo_c, PSF)#photo_c#
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr3)
        #wfits(dir1+name_f+band[k]+'e.fits',imag_E1,hdr)
        #wfits(dir1+name_f+band[k]+'_Fe.fits',imag_E2,hdr2)
        #wfits(dir1+name_f+band[k]+'_Le.fits',imag_E3,hdr3)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#0.75e-14)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    from PIL import Image
    file00u=dir1+name_f+"u.fits"
    file00g=dir1+name_f+"g.fits"
    file00r=dir1+name_f+"rg.fits"
    file00i=dir1+name_f+"ig.fits"
    file00z=dir1+name_f+"z.fits"
    file00B=dir1+name_f+"B.fits"
    file00V=dir1+name_f+"V.fits"
    file00R=dir1+name_f+"R.fits"
    file00I=dir1+name_f+"I.fits"    
    [pdl_00u,hdrt]=gdata(file00u, 0, header=True)
    [pdl_00g,hdrt]=gdata(file00g, 0, header=True)
    [pdl_00r,hdrt]=gdata(file00r, 0, header=True)
    [pdl_00i,hdrt]=gdata(file00i, 0, header=True)
    [pdl_00z,hdrt]=gdata(file00z, 0, header=True)
    [pdl_00B,hdrt]=gdata(file00B, 0, header=True)
    [pdl_00V,hdrt]=gdata(file00V, 0, header=True)
    [pdl_00R,hdrt]=gdata(file00R, 0, header=True)
    [pdl_00I,hdrt]=gdata(file00I, 0, header=True)
    max=np.amin(pdl_00r[np.where(pdl_00r > 0)])-0.5#24.0#17.5
    min=25.0
    nx,ny=pdl_00g.shape
    pdl_00u=(np.flipud(pdl_00u)-min)/(max-min)*256
    pdl_00g=(np.flipud(pdl_00g)-min)/(max-min)*256
    pdl_00r=(np.flipud(pdl_00r)-min)/(max-min)*256
    pdl_00i=(np.flipud(pdl_00i)-min)/(max-min)*256
    pdl_00z=(np.flipud(pdl_00z)-min)/(max-min)*256
    pdl_00B=(np.flipud(pdl_00B)-min)/(max-min)*256
    pdl_00V=(np.flipud(pdl_00V)-min)/(max-min)*256
    pdl_00R=(np.flipud(pdl_00R)-min)/(max-min)*256
    pdl_00I=(np.flipud(pdl_00I)-min)/(max-min)*256
    pdl_00u[np.where(pdl_00u < 0)]=0
    pdl_00g[np.where(pdl_00g < 0)]=0
    pdl_00r[np.where(pdl_00r < 0)]=0
    pdl_00i[np.where(pdl_00i < 0)]=0
    pdl_00z[np.where(pdl_00z < 0)]=0
    pdl_00B[np.where(pdl_00B < 0)]=0
    pdl_00V[np.where(pdl_00V < 0)]=0
    pdl_00R[np.where(pdl_00R < 0)]=0
    pdl_00I[np.where(pdl_00I < 0)]=0
    pdl_00u[np.where(pdl_00u > 255)]=255
    pdl_00g[np.where(pdl_00g > 255)]=255
    pdl_00r[np.where(pdl_00r > 255)]=255
    pdl_00i[np.where(pdl_00i > 255)]=255
    pdl_00z[np.where(pdl_00z > 255)]=255
    pdl_00B[np.where(pdl_00B > 255)]=255
    pdl_00V[np.where(pdl_00V > 255)]=255
    pdl_00R[np.where(pdl_00R > 255)]=255
    pdl_00I[np.where(pdl_00I > 255)]=255
    pdl_00=np.zeros([nx,ny,3],dtype="uint8")
    pdl_00[:,:,0]=pdl_00i
    pdl_00[:,:,1]=pdl_00r
    pdl_00[:,:,2]=pdl_00g
    im = Image.fromarray(pdl_00)
    im.save(dir1+name_f+"gri.jpeg",quality=100)
    pdl_00[:,:,0]=pdl_00R
    pdl_00[:,:,1]=pdl_00V
    pdl_00[:,:,2]=pdl_00B
    im1 = Image.fromarray(pdl_00)
    im1.save(dir1+name_f+"BVR.jpeg",quality=100)
    
def band_val(value,wave_s,dir='',k=3):
    pdl_flux=np.copy(value)
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt']
    f=open(dir+file[k],'r')
    wave=[]
    trans=[]
    for line in f:
        line=line.replace('\n','')
        data=line.split(' ')
        data=filter(None,data)
        if len(data) > 1:
            wave.extend([float_(data[1])])
            trans.extend([float_(data[2])])
    f.close()
    trans=np.array(trans)
    wave=np.array(wave)
    spec=pdl_flux
    spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
    flux_t=spec1*trans
    f_fin=simpson_r(flux_t,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
    photo_a=f_fin
    return photo_a       
    
    
def band_mag(pdl_flux1,crpix,cdelt,crval,dir='',k=3):
    pdl_flux=np.copy(pdl_flux1)
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    nw=len(pdl_flux)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt']
    band=['ha','u','g','rg','ig','z','B','V','I','R']
    zerop=[3631.0,3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]#3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    f=open(dir+file[k],'r')
    wave=[]
    trans=[]
    for line in f:
        line=line.replace('\n','')
        data=line.split(' ')
        data=filter(None,data)
        if len(data) > 1:
            wave.extend([float_(data[1])])
            trans.extend([float_(data[2])])
    f.close()
    d_wave=np.zeros(len(wave))
    for kk in range(1,len(wave)):
        d_wave[kk]=wave[kk]-wave[kk-1]
    d_wave[0]=d_wave[1]
    trans=np.array(trans)
    wave=np.array(wave)
    spec=pdl_flux*1e-16
    if np.sum(spec) > 0:
        spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
        flux_t=spec1*trans*d_wave
        f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
        if f_fin <= 0:
            photo_a=-1000.0
        else:
            photo_a=-2.5*np.log10(f_fin)-40.0
    else:
        photo_a=-1000.0
    return photo_a   
        
    
def id_str(id,n_z=2):
    id=int(np.float_(id))
    if n_z < 2 or n_z > 8:
        n_z=2
    if n_z == 2:
        if id < 10:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 3:
        if id < 10:
            idt='00'+str(id)
        elif id < 100:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 4:        
        if id < 10:
            idt='000'+str(id)
        elif id < 100:
            idt='00'+str(id)
        elif id < 1000:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 5:        
        if id < 10:
            idt='0000'+str(id)
        elif id < 100:
            idt='000'+str(id)
        elif id < 1000:
            idt='00'+str(id)
        elif id < 10000:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 6:        
        if id < 10:
            idt='00000'+str(id)
        elif id < 100:
            idt='0000'+str(id)
        elif id < 1000:
            idt='000'+str(id)
        elif id < 10000:
            idt='00'+str(id)
        elif id < 100000:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 7:        
        if id < 10:
            idt='000000'+str(id)
        elif id < 100:
            idt='00000'+str(id)
        elif id < 1000:
            idt='0000'+str(id)
        elif id < 10000:
            idt='000'+str(id)
        elif id < 100000:
            idt='00'+str(id)
        elif id < 1000000:
            idt='0'+str(id)
        else:
            idt=str(id)
    elif n_z == 8:        
        if id < 10:
            idt='0000000'+str(id)
        elif id < 100:
            idt='000000'+str(id)
        elif id < 1000:
            idt='00000'+str(id)
        elif id < 10000:
            idt='0000'+str(id)
        elif id < 100000:
            idt='000'+str(id)
        elif id < 1000000:
            idt='00'+str(id)
        elif id < 10000000:
            idt='0'+str(id)
        else:
            idt=str(id)
    return idt


def mock_test(fib_n,typef1="MaNGA",id1='A2-0',psf=0,dir1=''):
    sig=2.5
    thet=0.0
    plots=1
    nl=110
    n_pix=440
    cam=(142135.5*2.0)*7.0/np.float_(fib_n)#/2.34 # 2.34 sp6#/3.0#5.5#/2.2#IF CALIFA 5.5
    fov_p=np.round(142135.5/np.abs(cam)*91.0)
    fov=np.round(142135.5/np.abs(cam)*62.0)#30
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=300.0#150.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name='test-'+id1
    print name
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    Om=0.2726
    Lam=0.7274
    ho=0.704
    mock_halo_test(id1,psf=psf,dir1=dir1,fib_n=fib_n,nl=nl,fov=fov,thet=thet,ifutype=typef1)
    
def mock_vel(id1,Om=0.2726,Lam=0.7274,fov_p=91.0,ang=0.0,ho=0.704,cam=142135.5,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat",npix=440):
    thet=0.0
    name=base_name+id1
    sycall('echo '+name)
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    R3=np.array([[1,0,0],[0,1,0],[0,0,1]])
    ex=np.array([1,0,0])
    Ev=np.transpose(np.array([vx,vy,vz]))
    obs_ang1=np.dot(np.dot(Ev,R3),ex)*cam 
    R3a=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]])
    exa=np.array([1,0,0])
    Eva=np.transpose(np.array([vx,vy,vz]))
    obs_ang2=np.dot(np.dot(Eva,R3a),exa)*cam 
    mock_simvel(id1,basename=base_name,file=file_red,file2=file_gas,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=npix,fov=fov_p,thet=thet,observer=obs_ang1,observer2=obs_ang2,ang=ang)
    

def mock_sp(fib_n,ang,sp_res=2000.0,sp_samp=1.25,modt=0,template_1="libs/gsd61_156.fits",template_2="libs/templete_gas.fits",template_3="libs/templete_bc03_5.fits",template="libs/templete_bc03_2.fits",n_pix=440,fov_p=0,fov=0,rx=142135.5,Om=0.2726,Lam=0.7274,ho=0.704,cam=0,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',typef1="MaNGA",id1='A2-0',psf=0,redo=0,SN=15.0,Fluxm=20.0,dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat",file_out='mock_mass_ill_0.out',file_out_f='mock_mass_ill.out'):
    thet=0.0
    plots=1
    nl=110
    if cam == 0:
        if "CALIFA" in typef1:
            fib_n=14
        cam=(rx*2.0)*7.0/np.float_(fib_n)#/2.34 # 2.34 sp6#/3.0#5.5#/2.2#IF CALIFA 5.5
    cam=-cam
    if fov_p == 0:
        fov_p=np.round(rx/np.abs(cam)*91.0)
    if fov == 0:
        fov=np.round(rx/np.abs(cam)*62.0)#30
    if "MaNGA" in typef1:
        if psf <= 0:
            sig=1.43
        else:
            sig=psf
    elif "CALIFA" in typef1:
        if psf <= 0:
            sig=0.7
        else:
            sig=psf
    elif "MUSE" in typef1:
        if psf <= 0:
            sig=0.6
        else:
            sig=psf
    else:
        if psf == 0:
            sig=1.43
        else:
            sig=psf  
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=300.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name=base_name+id1
    #print name
    sycall('echo '+name)
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    R3=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]])
    ex=np.array([1,0,0])
    Ev=np.transpose(np.array([vx,vy,vz]))
    obs_ang1=np.dot(np.dot(Ev,R3),ex)*cam 
    if modt == 0:
        mock_vel(id1,fov_p=fov_p,Om=Om,Lam=Lam,ho=ho,cam=cam,vx=vx,vy=vy,vz=vz,base_name=base_name,dir1=dir1,file_red=file_red,file_gas=file_gas,npix=n_pix,ang=ang)
        mock_photosim(id1,basename=base_name,file=file_red,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_photosimextgas(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_photo(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_sim(id1,basename=base_name,file=file_red,template2=template,template=template_1,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
        fit3d_only(name,dir_root=dir1,redo=1)#,redo=1)
    if modt == 1:
        mock_photo(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
        fit3d_only(name,dir_root=dir1,redo=1)#,redo=1)
    if modt == 2:
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        fit3d_only(name,dir_root=dir1,redo=1)#,redo=1)
    if modt == 3:
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    if modt == 4:
        mock_vel(id1,fov_p=fov_p,Om=Om,Lam=Lam,ho=ho,cam=cam,vx=vx,vy=vy,vz=vz,base_name=base_name,dir1=dir1,file_red=file_red,file_gas=file_gas,npix=n_pix,ang=ang)
    #fit3d_only(name,dir_root=dir1,redo=1)#,redo=1)
    #fit3d(name,fo,fot_par,dir2=dir1,reens=1)#,redo=1)
    #fit3d_cen(name,fa,reens=0)

    
def mock_cen_re(id1='A2-12700',dir1='',file_out='mock_mass_ill_0.out',n_mt=500,n_mt0=0):
    if ptt.exists(dir1+id1) == False:
        print "The File does not exist"
        sys.exit()
    name='artsp8-'+id1
    print name
    fa=open('FIT3D_c_std'+file_out,'w')                     
    fit3d_cen_stat(name,fa,n_m=n_mt,n_mo=n_mt0)
    fa.close()
    
    
def mock_re(fib_n,ang,template_2="libs/templete_gas.fits",template="libs/templete_bc03_2.fits",ns=0,fov_p=0,fov=0,rx=142135.5,n_pix=440,cam=0,Om=0.2726,Lam=0.7274,ho=0.704,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',typef1="MaNGA",id1='A2-0',redo=0,SN=15.0,dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat"):
    sig=2.5
    thet=0.0
    plots=1
    nl=110
    if cam == 0:
        cam=(rx*2.0)*7.0/np.float_(fib_n)#/2.2
    cam=-cam
    if fov_p == 0:
        fov_p=np.round(rx/np.abs(cam)*91.0)
    if fov == 0:
        fov=np.round(rx/np.abs(cam)*30.0)
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=150.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    if ns == 0:
        ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name=base_name+id1
    print name
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    #Om=0.2726
    #Lam=0.7274
    #ho=0.704
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    #f4=open('Ntemp','w')
    R3=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]])
    ex=np.array([1,0,0])
    Ev=np.array([vx,vy,vz])
    obs_ang1=np.dot(np.dot(Ev,R3),ex)*cam  
    #obs_ang1=[(+0.8972639*np.cos(ang*np.pi/180.0)+0.4414944*np.sin(ang*np.pi/180.0))*cam,(-0.07978683*np.cos(ang*np.pi/180.0)+0.1621534*np.sin(ang*np.pi/180.0))*cam,(-0.4342251*np.cos(ang*np.pi/180.0)+0.8824902*np.sin(ang*np.pi/180.0))*cam] 
    #obs_ang1=[(-0.7156755*np.cos(ang*np.pi/180.0)+0.6984330*np.sin(ang*np.pi/180.0))*cam,(-0.5130859*np.cos(ang*np.pi/180.0)-0.5257526*np.sin(ang*np.pi/180.0))*cam,(0.4738687*np.cos(ang*np.pi/180.0)+0.4855672*np.sin(ang*np.pi/180.0))*cam]
    mock_photosim(id1,file=file_red,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
    mock_photosimextgas(id1,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
    #f4.close()
    
def gen_special(fibf,flav,xr,yr,typ,n_bosr=5,n_apor=5,ra0=180.0,dec0=0.0,scp_s=60.482,config_id='00000',dir='./'):
    nt4=np.where(typ != 'Fiducial')[0]
    xrt=xr[nt4]
    yrt=yr[nt4]
    typf=typ[nt4]
    #flavf=flav[nt4]
    n_rob=len(nt4)
    f=open(dir+'special_target_'+config_id+'.txt','w')
    f.write('HEADER n_bosr '+str(n_bosr)+' \n')
    f.write('HEADER n_apor '+str(n_apor)+' \n')
    #print typf
    cont=500
    for i in range(0, n_rob):
        if typf[i] == 'BOSS' and flav[i] == 2:
            print typf[i],flav[i],i
            xo=xrt[i]
            yo=yrt[i]
            alp=ran.rand(n_bosr)*360.0
            bet=ran.rand(n_bosr)*180.0
            xn,yn=FPS_conf2(xo,yo,alp,bet)
            dec=xn*1e3/3600.0/scp_s+dec0
            ra=yn*1e3/3600.0/scp_s+ra0
            mag=ran.rand(n_bosr)*(22.2-16.0)+16.0
            for j in range(0, n_bosr):
                cont=cont+1
                star_n=star_pool(str='none')
                Av=dust_getval(ra[j],dec[j])*3.1
                #print n_bosr
                f.write(str(fibf[i])+' , '+str(ra[j])+' , '+str(dec[j])+' , '+str(xn[j])+' , '+str(yn[j])+' , '+str(alp[j])+' , '+str(bet[j])+' , '+str(mag[j])+' , '+star_n+' , '+str(Av)+' , '+str(j)+' , BOSS , '+str(cont)+' \n')
    for i in range(0, n_rob):
        if typf[i] == 'BA' and flav[i] == 2:
            xo=xrt[i]
            yo=yrt[i]
            alp=ran.rand(n_apor)*360.0
            bet=ran.rand(n_apor)*180.0
            xn,yn=FPS_conf2(xo,yo,alp,bet)
            dec=xn*1e3/3600.0/scp_s+dec0
            ra=yn*1e3/3600.0/scp_s+ra0
            mag=ran.rand(n_apor)*(13.0-7.0)+7.0
            for j in range(0, n_apor):
                cont=cont+1
                star_n=star_pool(str='none')
                Av=dust_getval(ra[j],dec[j])*3.1
                f.write(str(fibf[i])+' , '+str(ra[j])+' , '+str(dec[j])+' , '+str(xn[j])+' , '+str(yn[j])+' , '+str(alp[j])+' , '+str(bet[j])+' , '+str(mag[j])+' , '+star_n+' , '+str(Av)+' , '+str(j)+' , APOGEE , '+str(cont)+' \n')
            alp=ran.rand(n_bosr)*360.0
            bet=ran.rand(n_bosr)*180.0
            xn,yn=FPS_conf2(xo,yo,alp,bet)
            dec=xn*1e3/3600.0/scp_s+dec0
            ra=yn*1e3/3600.0/scp_s+ra0
            mag=ran.rand(n_bosr)*(22.2-16.0)+16.0
            for j in range(0, n_bosr):
                cont=cont+1
                star_n=star_pool(str='none')
                Av=dust_getval(ra[j],dec[j])*3.1
                f.write(str(fibf[i])+' , '+str(ra[j])+' , '+str(dec[j])+' , '+str(xn[j])+' , '+str(yn[j])+' , '+str(alp[j])+' , '+str(bet[j])+' , '+str(mag[j])+' , '+star_n+' , '+str(Av)+' , '+str(j+n_apor)+' , BOSAPO , '+str(cont)+' \n')   
    f.close()
    
def psf_moffat_flux(psf=1.45,beta=4.765,rfib=1,dt=0.0):
    alpha=psf/(2.0*np.sqrt(2.0**(1.0/beta)-1.0))
    n_r=1000
    r=np.arange(n_r+1)/np.float(n_r)*rfib
    tht=np.arange(n_r+1)/np.float(n_r)*2.0*np.pi
    psf_x=np.zeros(n_r+1)
    for i in range(0, n_r+1):
        #print np.sqrt(r**2.0+dt**2.0-2.0*r*dt*np.cos(tht[i])),tht[i]
        psf_r=(beta-1)/(np.pi*alpha**2.0)*(1.0+((r**2.0+dt**2.0-2.0*r*dt*np.cos(tht[i]))/alpha**2.0))**(-beta)*r
        psf_x[i]=simpson_r(psf_r,r,0,len(r)-2,typ=0)
    flux=simpson_r(psf_x,tht,0,len(tht)-2,typ=0)
    return flux
    
def psf_moffat(nt,psf=1.45,beta=4.765):
    # beta=4.765 see Trujillo et al 2001
    from numpy import random as ran
    #beta=4
    alpha=psf/(2.0*np.sqrt(2.0**(1.0/beta)-1.0))
    n_t=500
    rf=np.zeros(n_t)
    ff=np.zeros(n_t)
    n_r=1000
    rmin=-20.0
    rmax=+20.0
    n_r1=500
    x=alpha*(np.arange(n_r1+1)/np.float(n_r1)*(rmax*1.1-rmin*1.1)+rmin*1.1)
    y=alpha*(np.arange(n_r+1)/np.float(n_r)*(rmax*1.1-rmin*1.1)+rmin*1.1)
    psf_x=np.zeros(n_r1+1)
    for i in range(0, n_r1+1):
        psf_xy=(beta-1)/(np.pi*alpha**2.0)*(1.0+((x[i]**2.0+y**2.0)/(alpha**2.0)))**(-beta)
        psf_x[i]=simpson_r(psf_xy,y,0,len(y)-2,typ=0)
        #psf_x[i]=(beta-1)/(np.pi*alpha**2.0)*(1.0+((x[i]**2.0)/(alpha)**2.0))**(-beta)    
    for i in range(0, n_t):
        ro=alpha*(20.0*i/np.float(n_t))+0.0001
        #ro=alpha*((rmax-rmin)*i/np.float(n_t)+rmin)+0.0001
        r=np.arange(n_r+1)/np.float(n_r)*ro
        #r=np.arange(n_r+1)/np.float(n_r)*(ro-rmin)+rmin
        #psf_r=(beta-1)/(np.pi*alpha**2.0)*(1.0+(r/alpha)**2.0)**(-beta)*2.0*np.pi
        psf_r=interp1d(x, psf_x,bounds_error=False,fill_value=0.)(r)*2.0*np.pi
        f_fin=simpson_r(psf_r,r,0,len(r)-2,typ=0)
        rf[i]=ro
        ff[i]=f_fin
        #print ro, f_fin, np.amax(x)
    f_rand=ran.rand(nt)*(np.amax(ff)-np.amin(ff))+np.amin(ff)
    s_rand=ran.rand(nt)*200.0-100.0
    s_rand[np.where(s_rand <= 0)[0]]=-1.0
    s_rand[np.where(s_rand  > 0)[0]]=+1.0
    r_rand=interp1d(ff, rf,bounds_error=False,fill_value=0.)(f_rand)
    r_rand=r_rand*s_rand
    #r=np.arange(n_r+1)/np.float(n_r)*1
    #psf_r=(beta-1)/(np.pi*alpha**2.0)*(1.0+(r/alpha)**2.0)**(-beta)*2.0*np.pi*r
    #print 1/simpson_r(psf_r,r,0,len(r)-2,typ=0)
    #ro=alpha*20.0+0.0001
    #r=np.arange(n_r+1)/np.float(n_r)*(ro-0.0001)+0.0001
    #psf_r1=psf_x*2.0*np.pi*.49/0.42*0.49/0.833
    #psf_r=(beta-1)/(np.pi*alpha**2.0)*(1.0+(r/alpha)**2.0)**(-beta)*2.0*np.pi*.49/0.42
    #import matplotlib.pyplot as plt
    #plt.xlim(-4,4)
    #plt.plot(r,psf_r/4.0,'-')
    #plt.plot(-r,psf_r/4.0,'-')
    #plt.plot(x,psf_r1/4.0,'-')
    #plt.hist(r_rand, bins=1000,range=(-4,4),normed=True) 
    ##plt.hist(f_rand, bins=100)#,range=(0,5),normed=True) 
    ##plt.plot(rf,ff)
    #plt.show() 
    return r_rand
    
def psf_pdf(nt,psf_m=1.45,d1=10.0,d2=200.0):
    from numpy import random as ran
    from scipy.special import betaincinv#gammaincinv
    psf_0=1.0/1.45*psf_m
    t=(ran.rand(nt))
    t1=betaincinv(d1/2.0,d2/2.0,t)
    t2=psf_0+d2/d1*t1/(1-t1)/(d2/(d2+2.0)*(d1-2.0)/d1)*(psf_m-psf_0)
    return t2
    
def read_special(config_id='0000',dir='./'):
    fib_r=[]
    type_r=[]
    mag_r=[]
    alpha_r=[]
    beta_r=[]
    id_r=[]
    av_r=[]
    x_r=[]
    y_r=[]
    ra_r=[]
    dec_r=[]
    obid_r=[]
    if ptt.exists(dir+'special_target_'+config_id+'.txt') == False:
        print 'No special target file'
    else:
        f=open(dir+'special_target_'+config_id+'.txt','r')
        for line in f:
            if 'HEADER' in line:
                if 'n_bosr' in line:
                    data=line.replace('\n','').split(' ')
                    data=filter(None,data)
                    n_bosr=np.int(data[2])
                if 'n_apor' in line:
                    data=line.replace('\n','').split(' ')
                    data=filter(None,data)
                    n_apor=np.int(data[2])
            if not 'HEADER' in line:
                if 'BOSS' in line:
                    data=line.replace('\n','').split(',')
                    data=filter(None,data)
                    #print data[10]
                    if np.float(data[10]) == 0.0:
                        randc=np.round((ran.rand(1)*n_bosr))[0]
                    randb=ran.rand(1)[0]*100.0
                    if np.float(data[10]) == randc and randb < -100.0:
                        fib_r.extend([np.int(data[0])])
                        ra_r.extend([np.float(data[1])])
                        dec_r.extend([np.float(data[2])])
                        x_r.extend([np.float(data[3])])
                        y_r.extend([np.float(data[4])])
                        alpha_r.extend([np.float(data[5])])
                        beta_r.extend([np.float(data[6])])
                        mag_r.extend([np.float(data[7])])
                        type_r.extend([data[8]])
                        av_r.extend([np.float(data[9])])
                        id_r.extend(['BOSS'])
                        obid_r.extend([np.int(data[12])])
                if 'APOGEE' in line or 'BOSAPO' in line:
                    data=line.replace('\n','').split(',')
                    data=filter(None,data)
                    if np.float(data[10]) == 0.0:
                        randc=np.round((ran.rand(1)*(n_apor)))[0]
                        #print randc
                    randb=ran.rand(1)[0]*100.0
                    if np.float(data[10]) == randc and randb < 1000.0 and 'APOGEE' in line:
                        fib_r.extend([np.int(data[0])])
                        ra_r.extend([np.float(data[1])])
                        dec_r.extend([np.float(data[2])])
                        x_r.extend([np.float(data[3])])
                        y_r.extend([np.float(data[4])])
                        alpha_r.extend([np.float(data[5])])
                        beta_r.extend([np.float(data[6])])
                        mag_r.extend([np.float(data[7])])
                        type_r.extend([data[8]])
                        av_r.extend([np.float(data[9])])
                        id_r.extend(['APOGEE'])
                        obid_r.extend([np.int(data[12])])
                    #else:
                    #    fib_r.extend([np.int(data[0])])
                    #    ra_r.extend([np.float(data[1])])
                    #    dec_r.extend([np.float(data[2])])
                    #    x_r.extend([np.float(data[3])])
                    #    y_r.extend([np.float(data[4])])
                    #    alpha_r.extend([np.float(data[5])])
                    #    beta_r.extend([np.float(data[6])])
                    #    mag_r.extend([np.float(data[7])])
                    #    type_r.extend([data[8]])
                    #    av_r.extend([np.float(data[9])])
                    #    id_r.extend(['APOGEE'])
                    #    obid_r.extend([np.int(data[12])])                        
                    if np.float(data[10]) == randc and randb < 1000.0 and 'BOSAPO' in line:
                        fib_r.extend([np.int(data[0])])
                        ra_r.extend([np.float(data[1])])
                        dec_r.extend([np.float(data[2])])
                        x_r.extend([np.float(data[3])])
                        y_r.extend([np.float(data[4])])
                        alpha_r.extend([np.float(data[5])])
                        beta_r.extend([np.float(data[6])])
                        mag_r.extend([np.float(data[7])])
                        type_r.extend([data[8]])
                        av_r.extend([np.float(data[9])])
                        id_r.extend(['BOSS'])
                        obid_r.extend([np.int(data[12])])
        f.close()
        fib_r=np.array(fib_r)
        ra_r=np.array(ra_r)
        dec_r=np.array(dec_r)
        x_r=np.array(x_r)
        y_r=np.array(y_r)
        type_r=np.array(type_r)
        mag_r=np.array(mag_r)
        alpha_r=np.array(alpha_r)
        beta_r=np.array(beta_r)
        id_r=np.array(id_r)
        av_r=np.array(av_r)
        obid_r=np.array(obid_r)
    return fib_r,ra_r,dec_r,x_r,y_r,type_r,mag_r,alpha_r,beta_r,id_r,av_r,obid_r 

def read_obssumary(file):
    f=open(file,'r')
    fiberid=[]
    fibertype=[]
    objectype=[]
    mag=[]
    for line in f:
        if 'FIBERMAP' in line:
            if len(line.split(' ')) > 10:
                data=line.replace('/n','').split(' ')
                data=filter(None,data)
                if np.int(data[1]) >= 0:
                    fiberid.extend([np.int(data[1])])
                    fibertype.extend([data[3]])
                    objectype.extend([data[5]])
                    mag.extend([np.array([np.float(data[13]),np.float(data[14]),np.float(data[15]),np.float(data[16]),np.float(data[17])])])
    return fiberid,fibertype,objectype,mag     
    
def create_photo_files(fieldt='0000',dir1='./',mjd='56008',conf='000000'):
    from pyfits import Column
    from pyfits import ColDefs
    from pyfits import BinTableHDU
    logfile='Log_obs-'+conf+'-'+mjd+'.log'
    obsumary='obsSummary-'+conf+'-'+mjd+'-01.par'
    print logfile
    ra_f,dec_f,x_f,y_f,alp_f,bet_f,fibf,typ,mago,Avo,obj_id,obj_cor,obj_typ,ra_0,dec_0,Ntar,Nobj,Nstd,Nsky,survey,junk=read_logreobsfile_bhm(dir1+'/bhm/'+logfile)
    fiberid,fibertype,objectype,mag=read_obssumary(dir1+'/bhm/'+mjd+'/raw_mock/'+obsumary)
    cdelt_w=1.25
    crval_w=1522.0
    crpix_w=1
    wave=np.arange(crval_w,10352.0,cdelt_w)
    dust_rat_ssp=A_l(3.1,wave)
    objid=[]
    skyver=np.ones(len(Avo),dtype='int8')
    mode=skyver
    clean=skyver
    run=skyver*2964
    rerut=skyver*300
    rerun=np.array([str(a) for a in rerut])
    camcol=skyver
    field=skyver*200
    thing=field
    thing[:]=-1
    parent=skyver
    nchild=skyver
    obj_type=skyver*0
    obj_prob=np.zeros(len(Avo))
    obj_flags=skyver
    obj_row=np.zeros(len(Avo))
    obj_rowerr=np.zeros(len(Avo))
    obj_col=np.zeros(len(Avo))
    obj_colerr=np.zeros(len(Avo))
    row_deg=np.zeros(len(Avo))
    row_degerr=np.zeros(len(Avo))
        
    flags2=[]
    
    Extin=[]
    airm=[]
    phi_offset=[]
    phi_dev_deg=[]
    phi_exp_deg=[]
    fiberflux=[]
    fiber2flux=[]
    fiberflux_ivar=[]
    fiber2flux_ivar=[]
    fibermag=[]
    fiber2mag=[]
    fibermag_err=[]
    fiber2mag_err=[]
    for i in range(0, len(Avo)):
        if typ[i] == 'NA':
            objid.extend([str(obj_id[i])])
        else:
            objid.extend([str(fiberid[i])])
        flags2.extend([np.array([0,0,0,0,0])])
        dust=10**(-0.4*Avo[i]*dust_rat_ssp)
        flux1=np.ones(len(wave))*10.0
        flux=flux1*dust
        M_0=band_mag(flux*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=1)
        M_1=band_mag(flux*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=2)
        M_2=band_mag(flux*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
        M_3=band_mag(flux*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=4)
        M_4=band_mag(flux*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=5)
        Mi_0=band_mag(flux1*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=1)
        Mi_1=band_mag(flux1*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=2)
        Mi_2=band_mag(flux1*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=3)
        Mi_3=band_mag(flux1*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=4)
        Mi_4=band_mag(flux1*1e-16,crpix_w,cdelt_w,crval_w,dir='legacy/',k=5)
        Ext1=M_0-Mi_0
        Ext2=M_1-Mi_1
        Ext3=M_2-Mi_2
        Ext4=M_3-Mi_3
        Ext5=M_4-Mi_4
        Extin.extend([np.array([Ext1,Ext2,Ext3,Ext4,Ext5])])
        airm.extend([np.array([1,1,1,1,1])])
        phi_offset.extend([np.array([0,0,0,0,0])])
        phi_dev_deg.extend([np.array([0,0,0,0,0])])
        phi_exp_deg.extend([np.array([0,0,0,0,0])])
        fiberflux.extend([10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[1.343, 1.336, 1.354, 1.363, 1.367]])
        fiber2flux.extend([10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[2.085, 2.085, 2.116, 2.134, 2.135]])
        fiberflux_ivar.extend([10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[1.343**2.0, 1.336**2.0, 1.354**2.0, 1.363**2.0, 1.367**2.0]*[0.1,0.1,0.1,0.1,0.1]])
        fiber2flux_ivar.extend([10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[2.085**2.0, 2.085**2.0, 2.116**2.0, 2.134**2.0, 2.135**2.0]*[0.1,0.1,0.1,0.1,0.1]])
        fibermag.extend([22.5-2.5*np.log10(10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[1.343, 1.336, 1.354, 1.363, 1.367])])
        fibermag_err.extend([22.5-2.5*np.log10(10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[1.343**2.0, 1.336**2.0, 1.354**2.0, 1.363**2.0, 1.367**2.0]*[0.1,0.1,0.1,0.1,0.1])])
        fiber2mag.extend([22.5-2.5*np.log10(10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[2.085, 2.085, 2.116, 2.134, 2.135])])
        fiber2mag_err.extend([22.5-2.5*np.log10(10.0**((22.5 - mag[i]) / 2.5)*[1.25, 1.0, 1.0, 1.0, 1.0]/[2.085**2.0, 2.085**2.0, 2.116**2.0, 2.134**2.0, 2.135**2.0]*[0.1,0.1,0.1,0.1,0.1])])
    objid=np.array(objid)
    Extin=np.array(Extin)
    airm=np.array(airm)
    phi_offset=np.array(phi_offset)
    phi_dev_deg=np.array(phi_dev_deg)
    phi_exp_deg=np.array(phi_exp_deg)
    fiberflux=np.array(fiberflux)
    fiber2flux=np.array(fiber2flux)
    fiberflux_ivar=np.array(fiberflux_ivar)
    fiber2flux_ivar=np.array(fiber2flux_ivar)
    fibermag=np.array(fibermag)
    fiber2mag=np.array(fiber2mag)
    fibermag_err=np.array(fibermag_err)
    fiber2mag_err=np.array(fiber2mag_err)
    flags2=np.array(flags2)
    
    t01 = Column(name='FLUX_OBJID', format='19A', array=objid) 
    t02 = Column(name='RUN', format='J', array=run)    
    t03 = Column(name='CAMCOL', format='J', array=camcol)
    t04 = Column(name='FIELD', format='J', array=field)
    t05 = Column(name='FLUX_ID', format='J', array=field)
    t06 = Column(name='THING_ID', format='J', array=thing)
    t07 = Column(name='RERUN', format='3A', array=rerun)
    t08 = Column(name='NOBSERVE', format='J', array=field)
    t09 = Column(name='NODETECT', format='J', array=field)
    t10 = Column(name='FLUX_RA', format='D', array=ra_f)
    t11 = Column(name='FLUX_DEC', format='D', array=dec_f)
    t12 = Column(name='MATCH_RA', format='D', array=ra_f)
    t13 = Column(name='MATCH_DEC', format='D', array=dec_f)
    t14 = Column(name='ORIG_OBJID', format='19A', array=objid) 
    t15 = Column(name='ORIG_RA', format='D', array=ra_f)
    t16 = Column(name='ORIG_DEC', format='D', array=dec_f)
    t17 = Column(name='APERFLUX3_R', format='E', array=ra_f*0)
    t18 = Column(name='FIBERFLUX_R', format='E', array=dec_f*0)
    t19 = Column(name='APERFLUX3_R_TOTAL', format='E', array=ra_f*0)
    t20 = Column(name='APERFLUX3_R_MATCH', format='E', array=dec_f*0)
    t21 = Column(name='FLUXMATCH_PARENT', format='J', array=field*0)
    t22 = Column(name='FLUXMATCH_STATUS', format='J', array=field*0)
    
    a01 = Column(name='OBJID', format='19A', array=objid)
    a02 = Column(name='PARENTID', format='19A', array=objid)
    a03 = Column(name='FIELDID', format='19A', array=objid)
    a04 = Column(name='SKYVERSION', format='B', array=skyver)
    a05 = Column(name='MODE', format='B', array=mode)
    a06 = Column(name='CLEAN', format='B', array=clean)
    a07 = Column(name='RUN', format='I', array=run)    
    a08 = Column(name='RERUN', format='3A', array=rerun)
    a09 = Column(name='CAMCOL', format='B', array=camcol)
    a10 = Column(name='FIELD', format='I', array=field)
    a11 = Column(name='PARENT', format='I', array=parent)
    a12 = Column(name='NCHILD', format='I', array=nchild)
    a13 = Column(name='OBJC_TYPE', format='J', array=obj_type)
    a14 = Column(name='OBJC_PROB_PSF', format='E', array=obj_prob)
    a15 = Column(name='OBJC_FLAGS', format='J', array=obj_flags)
    a16 = Column(name='OBJC_FLAGS2', format='J', array=obj_flags)
    a17 = Column(name='OBJC_ROWC', format='E', array=obj_row)
    a18 = Column(name='OBJC_ROWCERR', format='E', array=obj_rowerr)
    a19 = Column(name='OBJC_COLC', format='E', array=obj_col)
    a20 = Column(name='OBJC_COLCERR', format='E', array=obj_colerr)
    a21 = Column(name='ROWVDEG', format='E', array=row_deg)
    a22 = Column(name='ROWVDEGERR', format='E', array=row_degerr)
    a23 = Column(name='COLVDEG', format='E', array=row_deg)
    a24 = Column(name='COLVDEGERR', format='E', array=row_degerr)
    b01 = Column(name='FLAGS', format='5J', array=flags2)
    b02 = Column(name='FLAGS2', format='5J', array=flags2)
    b03 = Column(name='RA', format='D', array=ra_f)
    b04 = Column(name='DEC', format='D', array=dec_f)
    c01 = Column(name='ARIMASS', format='5E', array=airm)
    c02 = Column(name='PHI_OFFSET', format='5E', array=phi_offset)
    c03 = Column(name='PHI_DEV_DEG', format='5E', array=phi_dev_deg)
    c04 = Column(name='PHI_EXP_DEG', format='5E', array=phi_exp_deg)
    c05 = Column(name='EXTINTION', format='5E', array=Extin)
    c06 = Column(name='SKYFLUX', format='5E', array=fiberflux)
    c07 = Column(name='SKYFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c08 = Column(name='PSFFLUX', format='5E', array=fiberflux)
    c09 = Column(name='PSFFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c10 = Column(name='PSFMAG', format='5E', array=fibermag)
    c11 = Column(name='PSFMAGERR', format='5E', array=fibermag_err)
    c12 = Column(name='FIBERFLUX', format='5E', array=fiberflux)
    c13 = Column(name='FIBERFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c14 = Column(name='FIBERMAG', format='5E', array=fibermag)
    c15 = Column(name='FIBERMAGERR', format='5E', array=fibermag_err)
    c16 = Column(name='FIBER2FLUX', format='5E', array=fiber2flux)
    c17 = Column(name='FIBER2FLUX_IVAR', format='5E', array=fiber2flux_ivar)
    c18 = Column(name='FIBER2MAG', format='5E', array=fiber2mag)
    c19 = Column(name='FIBER2MAGERR', format='5E', array=fiber2mag_err)
    c20 = Column(name='CMODELFLUX', format='5E', array=fiberflux)
    c21 = Column(name='CMODELFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c22 = Column(name='CMODELMAG', format='5E', array=fibermag)
    c23 = Column(name='CMODELMAGERR', format='5E', array=fibermag_err)
    c24 = Column(name='MODELFLUX', format='5E', array=fiberflux)
    c25 = Column(name='MODELFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c26 = Column(name='MODELMAG', format='5E', array=fibermag)
    c27 = Column(name='MODELMAGERR', format='5E', array=fibermag_err)
    c28 = Column(name='PETROFLUX', format='5E', array=fiberflux)
    c29 = Column(name='PETROFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c30 = Column(name='PETROMAG', format='5E', array=fibermag)
    c31 = Column(name='PETROMAGERR', format='5E', array=fibermag_err)
    c32 = Column(name='DEVFLUX', format='5E', array=fiberflux)
    c33 = Column(name='DEVFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c34 = Column(name='DEVMAG', format='5E', array=fibermag)
    c35 = Column(name='DEVMAGERR', format='5E', array=fibermag_err)
    c36 = Column(name='EXPFLUX', format='5E', array=fiberflux)
    c37 = Column(name='EXPFLUX_IVAR', format='5E', array=fiberflux_ivar)
    c38 = Column(name='EXPMAG', format='5E', array=fibermag)
    c39 = Column(name='EXPMAGERR', format='5E', array=fibermag_err)
    d01 = Column(name='THING_ID', format='J', array=thing)
    coldefs = ColDefs([a01,a02,a03,a04,a05,a06,a07,a08,a09,a10, 
                       a11,a12,a13,a14,a15,a16,a17,a18,a19,a20, 
                       a21,a22,a23,a24,b01,b02,b03,b04,c01,c02,c03,c04,
                       c05,c06,c07,c08,c09,c10,c11,c12,c13,c14,
                       c15,c16,c17,c18,c19,c20,c21,c22,c23,c24,
                       c25,c26,c27,c28,c29,c30,c31,c32,c33,c34,
                       c35,c36,c37,c38,c39,d01])
    tbhdu = BinTableHDU.from_columns(coldefs)
    sycall('rm '+dir1+'/bhm/'+mjd+'/raw_mock/'+'photoField-'+fieldt+'.fits')
    tbhdu.writeto(dir1+'/bhm/'+mjd+'/raw_mock/'+'photoField-'+fieldt+'.fits')
    sycall('cp '+dir1+'/bhm/'+mjd+'/raw_mock/'+'photoField-'+fieldt+'.fits '+dir1+'/bhm/'+mjd+'/raw_mock/'+'photoPosField-'+fieldt+'.fits ')
    coldefs1 = ColDefs([t01,t02,t03,t04,t05,t06,t07,t08,t09,t10, 
                        t11,t12,t13,t14,t15,t16,t17,t18,t19,t20, 
                        t21,t22])
    tbhdu1 = BinTableHDU.from_columns(coldefs1)
    sycall('rm '+dir1+'/bhm/'+mjd+'/raw_mock/'+'photoMatchField-'+fieldt+'.fits')
    tbhdu1.writeto(dir1+'/bhm/'+mjd+'/raw_mock/'+'photoMatchField-'+fieldt+'.fits')    
        
def deg_sexagesimal(ra,dec):
    ra_1=np.copy(ra)*24./360.
    ra_h=np.int(ra_1)    
    ra_2=(ra_1-ra_h)*60.0
    ra_m=np.int(ra_2)
    ra_s=np.round((ra_2-ra_m)*60.0)
    if dec < 0:
        sig='-'
    else:
        sig='+'
    dec_1=np.abs(np.copy(dec))
    dec_d=np.int(dec_1)    
    dec_2=(dec_1-dec_d)*60.0
    dec_m=np.int(dec_2)
    dec_s=np.round((dec_2-dec_m)*60.0)
    str_id=id_str(ra_h,n_z=3)+id_str(ra_m,n_z=2)+id_str(ra_s,n_z=2)+sig+id_str(dec_d,n_z=3)+id_str(dec_m,n_z=2)+id_str(dec_s,n_z=2)
    return str_id
    
def target_idform(ra,dec,typ=1):
    if typ == 1:
        typ_n='BHM'
    elif typ ==2:
        typ_n='MWM'
    elif typ ==3:
        typ_n='STD'
    elif typ ==4:
        typ_n='SKY'
    elif typ ==5:
        typ_n='BLK'
    else:
        typ_n='BHM'
    string_id=deg_sexagesimal(ra,dec)
    tar_fin=typ_n+string_id
    return tar_fin

def mock_bhm(dirtemp="libs/",reobs_f='none',special_t='none',dsep=1.0,mjdid=0,lshf=5070.0,prs=0,obsn=0,tfield=0,field=0,ra_0=180.0,dec_0=0.0,ra_0i=180.0,dec_0i=0.0,sp_res=2000.0,Av=0.00,sp_samp=1.25,basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1',template2="libs/sky_model.txt",template_1="libs/spEigenStar-55734.fits",template_wd="libs/da012500_800.dat",wd=False,base_name='sdR',config_name=0,typef1="APO",id1='A2-0',psf=0,SN=15.0,Fluxm=20.0,dir1='./',t_exp=0,all_boss=0,nsp=4):
    plots=0
    if "APO" in typef1:
        Nst=70
        Nsky=30
        Ntar=500
        if psf <= 0:
            sig=1.43
        else:
            sig=psf
    if "LCO" in typef1:
        Nst=70
        Nsky=30
        Ntar=500#330
        if psf <= 0:
            sig=0.63
        else:
            sig=psf
    else:
        if psf == 0:
            sig=1.43
        else:
            sig=psf  
    calib=False
    #leter=['A','B','C','D','E','F','G','H','I','J','K','L','M']
    #plate_name=np.str((config_name))#+'-'+
    if all_boss == 1:
        plate_name=id_str(np.float(config_name+tfield),n_z=6)
    else:
        plate_name=id_str(obsn+np.float(config_name),n_z=6)
    field_id=id_str(field,n_z=5)
    bhm_config_create(Nstd=Nst,Nsky=Nsky,Ntar=Ntar,dsep=dsep,special_t=special_t,mjdid=mjdid,lshf=lshf,ra_0=ra_0,obsn=obsn,dec_0=dec_0,ra_0i=ra_0i,dec_0i=dec_0i,prs=prs,reobs_f=reobs_f,dirtemp=dirtemp,basePath=basePath,sp_res=sp_res,Av=Av,sp_samp=sp_samp,template2=template2,template_1=template_1,template_wd=template_wd,wd=wd,base_name=base_name,plate_name=plate_name,field_id=field_id,typef1=typef1,id0=id1,psf=sig,SN=SN,Fluxm=Fluxm,dir1=dir1,t_exp=t_exp,all_boss=all_boss,nsp=nsp)
    create_photo_files(fieldt=field_id,dir1=dir1,mjd=id1,conf=plate_name)
    sys.exit()
    copy_data(id1,typef1,field_id)
    #sys.exit()

def bhm_config_create(Nstd=480,dsep=1.0,reobs_f='none',special_t='none',lshf=5070.0,Nsky=220,ra_0=180.0,obsn=0,dec_0=0.0,ra_0i=180.0,dec_0i=0.0,rw=1,prs=0,Ntar=500,basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1',dirtemp="libs/",sp_res=2000.0,Av=0.00,sp_samp=1.25,template2="libs/sky_model.txt",template_1="libs/spEigenStar-55734.fits",template_0="libs/gsd61_156.fits",template_2="libs/templete_gas.fits",template_3="libs/templete_bc03_5.fits",template_wd="libs/da012500_800.dat",wd=False,base_name='sdR',plate_name='0000',field_id='0000',typef1="APO",id0='A2-0',psf=0,SN=15.0,Fluxm=20.0,dir1='./',t_exp=0,all_boss=0,mjdid=0,nsp=4):
    sycall('mkdir -p '+dir1+'bhm')
    if not 'none' in reobs_f:
        print reobs_f
        if ptt.exists(dir1+'/bhm/'+reobs_f) == False:
            print 'No reobservation log file'
            re_obs = 0
        else:
            re_obs = 1
    else:
        re_obs =0
    if re_obs == 0:
        Nobj=Ntar-Nstd-Nsky
        Ntarc=Nstd+Nsky
        if "APO" in typef1:
            Nobjt=2000
            scp_s=60.482#microns per arcsec
            file_fps='fps_filledHex.txt'
            fibB=120.0
            Dp=3.0
            obs_l=32.78
        elif "LCO" in typef1:
            Nobjt=2000
            scp_s=91.475#microns per arcsec
            file_fps='fps_LCO330.txt'
            file_fps='fps_filledHex.txt'
            fibB=180.0
            Dp=2.0
            obs_l=-29.0146
        if prs == 0:
            typ_o=2
            obs0,ra_0,dec_0,raf,decf,xf,yf,ind0,ho,Om,Lam=cone_mock_ill(fov_p=Dp,rw=rw,dir1=dir1+'bhm/',basePath=basePath,scp_s=scp_s,ra_o=ra_0i,dec_o=dec_0i,the_o=ra_0,phi_o=dec_0)
            ind_t,out=check_obs_illus(raf,decf)
            if out == True:
                obs0=obs0[:,ind_t]  
                raf=raf[ind_t]
                decf=decf[ind_t]
                xf=xf[ind_t]
                yf=yf[ind_t]
                ind0=ind0[ind_t]
        else:
            ho=0.74
            Om=0.2726
            Lam=0.7274
            Ntarc=Ntar
            typ_o=3
        Rp=Dp/2.0
        indf=np.zeros(Nobjt*typ_o,dtype='int16')
        obsf=np.zeros([3,Nobjt*typ_o])
        rr=np.sqrt(ran.rand(Nobjt*typ_o))*Rp
        rt=ran.rand(Nobjt*typ_o)*2.0*np.pi
        ra=np.array(rr*np.sin(rt))+ra_0
        dec=np.array(rr*np.cos(rt))+dec_0
        y=(dec-dec_0)*3600.0*scp_s/1e3
        x=(ra-ra_0)*3600.0*scp_s/1e3
        #print obsf.shape
        if prs == 0:
            x=np.concatenate((x,xf))
            y=np.concatenate((y,yf))
            ra=np.concatenate((ra,raf))
            dec=np.concatenate((dec,decf))
            indf=np.concatenate((indf,ind0))
            obsf=np.concatenate((obsf,obs0),axis=1)
        #print obs0.shape
        #print obsf.shape      
        row,col,xr,yr,typ=fps_layout(file_n=file_fps)
        #nt1=np.where(typ == 'BA')[0]
        #nt2=np.where(typ == 'BOSS')[0]
        #nt3=np.where(typ == 'Fiducial')[0]
        #nt4=np.where(typ != 'Fiducial')[0]
        #print nt1
        #print xr[nt2]
        #print yr[nt2]
        fps_try=0
        while True:
            fps_try=fps_try+1
            xof,yof,alp,bet,ind,flav=def_FPS(x,y,xr,yr,typ,r1=7.4,r2=15.0,Nob=Nobj,Nsky=Nsky,Nstd=Nstd,Nobjt=Nobjt)
            if len(np.where(flav == -1)[0]) == 0:
                break 
            else:
                print len(np.where(flav == -1)[0]),fps_try,Nsky,Nobj
            if fps_try >= 500:
                tm_std=np.where(flav == 1)[0]
                tm_sky=np.where(flav == 0)[0]
                if len (tm_std) == Nstd and len (tm_sky) == Nsky:
                    break
                else:
                    fps_try=499
                #Nobj=Nobj-5
                #Nsky=Nsky+5
                #fps_try=0
        plot_fps_tar(xof,yof,alp,bet,xr,yr,typ,flav,r1=7.4,r2=15.0,name='Target0',dir=dir1+'bhm/')
        #print "DONE"
        #sys.exit()
        mjdf=np.float(id0)
        sycall('mkdir -p '+dir1+'bhm/plot_hole')
        plots_hol_type(mjdf,name=dir1+'bhm/plot_hole/'+typef1+'_plot_hole',ra_0=ra_0,dec_0=dec_0,phi=obs_l,scp_s=scp_s,fibB=fibB)
        #sys.exit()
        fibf=np.array(range(len(xof)))+1
        apo_t=np.zeros(len(xof))
        indf_f=indf[ind]
        obsf_f=obsf[:,ind]
        #sys.exit()  
        dec_f=dec[ind]
        ra_f=ra[ind]
        x_f=x[ind]
        y_f=y[ind]
        alp_f=alp
        bet_f=bet
        ind_std=np.where(flav == 1)[0]
        ind_tar=np.array(range(Ntar))
        ind_sky=np.where(flav == 0)[0]
        #Nstd=len(ind_std)
        #Nsky=len(ind_sky)
        #print Nstd,Nsky
        if special_t == 'none' and all_boss == 0:      
            gen_special(fibf,flav,xr,yr,typ,n_bosr=2,n_apor=6,ra0=ra_0,dec0=dec_0,scp_s=scp_s,config_id=plate_name,dir=dir1+'/bhm/')
        #else:
        #    gen_s=np.zeros(len(fibf))
        #    mago=np.zeros(len(fibf))
        #    fib_s,ra_s,dec_s,x_s,y_s,type_s,mag_s,alpha_s,beta_s,id_s,av_s=read_special(config_id=special_t,dir=dir1+'/bhm/')
        #    apo_t=np.zeros(len(fibf))
        #    if len(ra_s) > 0:
        #        for kk in range(0, len(fib_s)):
        #            for jj in range(0, len(fibf)):
        #                if fibf[jj] == fib_s[kk]:
        #                    if id_s[kk] == 'APOGEE':
        #                        apo_t[jj]=1.0
        #                    gen_s[jj]=1.0
        #                    ra_f[jj]=ra_s[kk]
        #                    dec_f[jj]=dec_s[kk]
        #                    x_f[jj]=x_s[kk]
        #                    y_f[jj]=y_s[kk]
        #                    #Avo[jj]=av_s[kk]
        #                    mago[jj]=mag_s[kk]
        #                    obj_id[jj]=type_s[kk]
        #                    obj_typ[jj]='star'
        #print len(np.where(flav == -1)[0])
        #print len(np.where(flav == 0)[0]),Nsky
        #sys.exit()
        #fibif[np.where(flav == -1)[0]]=-1
    else:
        ho=0.74
        Om=0.2726
        Lam=0.7274
        ra_f,dec_f,x_f,y_f,alp_f,bet_f,fibf,typ,mago,Avo,obj_id,obj_cor,obj_typ,ra_0,dec_0,Ntar,Nobj,Nstd,Nsky,junk,flav=read_logreobsfile_bhm(dir1+'/bhm/'+reobs_f)
        fib_s,ra_s,dec_s,x_s,y_s,type_s,mag_s,alpha_s,beta_s,id_s,av_s,obid_s=read_special(config_id=special_t,dir=dir1+'/bhm/')
        apo_t=np.zeros(len(fibf))
        objf=np.array(range(len(fibf)))+1
        #print len(ra_s)
        #sys.exit()
        if len(ra_s) > 0 and all_boss == 0:
            for kk in range(0, len(fib_s)):
                for jj in range(0, len(fibf)):
                    if not 'SKY' in typ[jj]:
                        if not 'SPECTROPHOTO_STD' in typ[jj]:
                            if fibf[jj] == fib_s[kk]:
                                if id_s[kk] == 'APOGEE':
                                    apo_t[jj]=1.0
                                ra_f[jj]=ra_s[kk]
                                dec_f[jj]=dec_s[kk]
                                x_f[jj]=x_s[kk]
                                y_f[jj]=y_s[kk]
                                alp_f[jj]=alpha_s[kk]
                                bet_f[jj]=beta_s[kk]
                                Avo[jj]=av_s[kk]
                                mago[jj]=mag_s[kk]
                                obj_id[jj]=type_s[kk]
                                obj_typ[jj]='star'
                                objf[jj]=obid_s[kk]
                                flav[jj]=2
        #print obj_typ
        #sys.exit()
    exptime_sc=900.0
    exptime_ar=4.0
    exptime_fl=25.0
    cmr=60/exptime_fl
    redot=55.6
    FPS_rect=120.0
    ttime_sc=exptime_sc+redot
    ttime_ar=exptime_ar+redot
    ttime_fl=exptime_fl+redot
    #ttime=4.0*ttime_sc+ttime_fl+ttime_ar
    obs_t=0
    if (obsn % nsp) > 0:# the previus setup use 4 exposures, then change 8 by 4
        obs_t=1
        calib=False
    else:
        calib=True
        
    if calib == True:
        nexp=3
        ttime=ttime_sc+ttime_fl+ttime_ar
        fex=1
    else:
        fex=0
        nexp=1
        ttime=ttime_sc
    ttime2=(nsp*ttime_sc+nsp*FPS_rect+ttime_fl+ttime_ar)
    #toffs=(obsn-3.5*4.0)*(ttime+FPS_rect)
    toffs=-13869.1+obsn*(ttime_sc+FPS_rect)+(obs_t+int(obsn/nsp))*(ttime_fl+ttime_ar)
    toffs2=-13869.1+int(obsn/nsp)*(nsp*ttime_sc+nsp*FPS_rect+ttime_fl+ttime_ar)
    #exp_r=1000+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6*8+6*10+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+6*9+nexp*(obsn+1)
    #exp_r=2506+42+42+42+42+42+42+42+42+42+42+42+42+42+42+42*(mjdid+1)
    #exp_r=2506+42+42+42+42+42+42+42+42+42+42+42+42+42+42+42+48*(mjdid)#+54
    exp_r=2506+42+42+42+42+42+42+42+42+42+42+42+42+42+42+42+48+324+480+(60*70)+(12*5)+10*(mjdid)#60*(mjdid)#+3*6#+54
    
    
    expoff=ttime_ar*fex+ttime_fl*fex+toffs
    ha=lst_c(tai=((np.float(id0)+.25)*3600.0*24.0+expoff+exptime_sc/2.0))-ra_0/180.0*12.0 
    
    dir0='bhm/'+id0#+'/'+plate_name.split('-')[0]#+'-'+mjd  
    sycall('mkdir -p '+dir1+'bhm') 
    sycall('mkdir -p '+dir1+'bhm/'+id0)
    #sycall('mkdir -p '+dir1+dir0)        
    dir0f=dir0+'/raw_mock'
    sycall('mkdir -p '+dir1+dir0f)  
    dir0f=dir0f+'/'
    while True:
        out_fit=dir1+dir0f+base_name+'-'+'b1-'+id_str(exp_r,n_z=8)+'.fit.gz'
        if ptt.exists(out_fit) == False:
            break
        else:
            exp_r=exp_r+1
    
    cot1=0
    cot2=0
    cot3=0
    if Ntar > 0:
        if calib == True:
            for kk in range(0,Ntar):
                #print kk
                #print x_f[kk],y_f[kk],xof[kk],yof[kk]
                spec_b,spec_r=mock_arc(dirtemp=dirtemp)
                if kk == 0:
                    blu_sp=np.zeros([len(spec_b),1])
                    red_sp=np.zeros([len(spec_r),1])
                    blu_sp[:,0]=spec_b
                    red_sp[:,0]=spec_r
                else:
                    blu_sp1=np.zeros([len(spec_b),1])
                    red_sp1=np.zeros([len(spec_r),1])
                    blu_sp1[:,0]=spec_b
                    red_sp1[:,0]=spec_r
                    blu_sp=np.concatenate((blu_sp,blu_sp1),axis=1)
                    red_sp=np.concatenate((red_sp,red_sp1),axis=1)
            #sys.exit()
#            raw_exp_bhm(blu_sp,fibf,base_name,fc=[0.88,0.94],n_cr=20,d_cr=2,type="blue",dir1=dir1,mjd=id0,plate=plate_name,flb='a',exp=exp_r+0,ra0=ra_0,dec0=dec_0,expof=toffs)
#            raw_exp_bhm(red_sp,fibf,base_name,fc=[0.83,0.89],n_cr=20,d_cr=2,type="red",dir1=dir1,mjd=id0,plate=plate_name,flb='a',exp=exp_r+0,ra0=ra_0,dec0=dec_0,expof=toffs)        
            for kk in range(0,Ntar):
                spec_b,spec_r=mock_flat(dirtemp=dirtemp)
                if kk == 0:
                    blu_sp=np.zeros([len(spec_b),1])
                    red_sp=np.zeros([len(spec_r),1])
                    blu_sp[:,0]=spec_b
                    red_sp[:,0]=spec_r
                else:
                    blu_sp1=np.zeros([len(spec_b),1])
                    red_sp1=np.zeros([len(spec_r),1])
                    blu_sp1[:,0]=spec_b
                    red_sp1[:,0]=spec_r
                    blu_sp=np.concatenate((blu_sp,blu_sp1),axis=1)
                    red_sp=np.concatenate((red_sp,red_sp1),axis=1)
#            raw_exp_bhm(blu_sp,fibf,base_name,fc=[0.88,0.94],n_cr=60,type="blue",dir1=dir1,mjd=id0,plate=plate_name,flb='f',exp=exp_r+1,ra0=ra_0,dec0=dec_0,expof=ttime_ar+toffs)
#            raw_exp_bhm(red_sp,fibf,base_name,fc=[0.83,0.89],n_cr=60,type="red",dir1=dir1,mjd=id0,plate=plate_name,flb='f',exp=exp_r+1,ra0=ra_0,dec0=dec_0,expof=ttime_ar+toffs)        
       # sys.exit()
        dir0='bhm/'+plate_name+"-"+id0        
        sycall('mkdir -p '+dir1+dir0)
        dir00='bhm/'+id0#+'/'+plate_name.split('-')[0]        
        dir0f=dir00+'/raw_mock'
        sycall('mkdir -p '+dir1+dir0f)        
        ft=open(dir1+dir0f+"/obsSummary-"+plate_name+"-"+id0+"-01.par","w")
        fl=open(dir1+"/bhm/Log_obs-"+plate_name+"-"+id0+".log","w")
        fl.write('$ RA  DEC  X  Y ALP  BET  FIB  TYPE MAG AV OBJ_ID  OBJ_COR OBJ_TYP OBS_TYP \n')
        fl.write('# '+str(ra_0)+' '+str(dec_0)+' '+str(Ntar)+' '+str(Nobj)+' '+str(Nstd)+' '+str(Nsky)+'\n')
        head_map_bhm(ft,plate_name,field_id,Nobj,Nsky,Nstd,id0,ha=ha,ra=ra_0,dec=dec_0,obs=typef1)
        for kk in range(0,Ntar):
            if re_obs == 0:
                kt=0
                for jj in range(0, Nstd):  
                    if ind_tar[kk] == ind_std[jj]:
                        kt=1
                        indt=ind_tar[kk]
                        type='SPECTROPHOTO_STD'
                        cot1=cot1+1
                for jj in range(0, Nsky):  
                    if ind_tar[kk] == ind_sky[jj]:
                        kt=2
                        indt=ind_tar[kk]
                        type='SKY'
                        cot2=cot2+1
                if kt == 0:
                    indt=ind_tar[kk]
                    type='NA'
                    cot3=cot3+1
            else:
                type=typ[kk]
                indt=kk
                #print type
            if fibf[indt] >= 500:
                samp=2
            else:
                samp=1
            id1=id0+'-'+id_str(fibf[indt],n_z=4)
            #if fibf[indt] % 20 == 1:
            #    t=ran.rand(1)[0]*100.0
            #    if t  < 5:
            if re_obs == 0:
                if flav[indt] == -1:
                    fib_t="-1"
                else:
                    fib_t=str(fibf[indt])
            else:
                if flav[indt] == -1:
                    fib_t="-1"
                else:
                    fib_t=str(fibf[indt])
            rob_id=fib_t
            if typef1 == 'APO':
                fib_si='120'
            else:
                fib_si='120'
            if type == 'NA':
                #print re_obs
                if re_obs == 0:
                    idh=indf_f[indt]
                    obs=obsf_f[:,indt]
                    star_n=star_pool(str='none')
                    Av=dust_getval(ra_f[indt],dec_f[indt])*3.1
                    #id_obj=id_str(fibf[indt],n_z=6)
                    id_obj_log=id_str(fibf[indt],n_z=18)
                    id_obj=target_idform(ra_f[indt],dec_f[indt],typ=1)
                else:
                    Av=Avo[kk]
                    #Av=dust_getval(ra_f[indt],dec_f[indt])*3.1
                    #print obj_typ[kk]
                    if obj_typ[kk] == 'galaxy':
                        idh=np.int(obj_id[kk])
                        prs=0
                        obs=obj_cor[kk]
                        #id_obj=id_str(idh,n_z=6)
                    if obj_typ[kk] == 'star':
                        prs=1
                        star_n=obj_id[kk]
                        #id_obj=id_str(objf[kk],n_z=6)
                    id_obj_log=id_str(objf[kk],n_z=18)
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                if apo_t[kk] == 0:
                    bostyp='BOSS'
                    if prs == 0:
                        id_obj=target_idform(ra_f[indt],dec_f[indt],typ=1)
                    else:                       
                        id_obj=target_idform(ra_f[indt],dec_f[indt],typ=2)
                else:
                    bostyp='APOGEE'
                    id_obj=target_idform(ra_f[indt],dec_f[indt],typ=5)
                if indt < 20:
                    plots=1
                else:
                    plots=0
                if prs == 0 and bostyp == 'BOSS':
                    spec_b1,spec_r1,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",ho=ho,Lam=Lam,Om=Om,plots=plots,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_f[indt],dec_0=dec_f[indt],expof=0.0*ttime_sc+ttime_ar*fex+ttime_fl*fex+toffs,expot=ttime2,toffs=toffs2,t_exp=t_exp)
                    #spec_b2,spec_r2,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    #spec_b3,spec_r3,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    #spec_b4,spec_r4,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(np.round(alp_f[indt],5))+' , '+str(np.round(bet_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(-1)+' , '+str(Av)+' , '+str(idh)+' , ['+str(obs[0])+'_'+str(obs[1])+'_'+str(obs[2])+'] , galaxy , '+bostyp+' , '+id_obj_log+' , '+str(flav[indt])+' \n'
                else: 
                    if re_obs == 0:  
                        mag_ob=ran.rand(1)*(22.2-17.0)+17.0
                        mag_ob=mag_ob[0]
                    else:
                        mag_ob=mago[kk]
                    #print star_n
                    spec_b1,spec_r1,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,alp=alp_f[indt],bet=bet_f[indt],star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,template_wd=template_wd,wd=wd,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=plots,ifutype=typef1,expt=exptime_sc,ra_0=ra_f[indt],dec_0=dec_f[indt],expof=0.0*ttime_sc+ttime_ar*fex+ttime_fl*fex+toffs,expot=ttime2,toffs=toffs2,t_exp=t_exp,apot=apo_t[kk])
                    #spec_b2,spec_r2,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp,apot=apo_t[kk])
                    #spec_b3,spec_r3,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp,apot=apo_t[kk])        
                    #spec_b4,spec_r4,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp,apot=apo_t[kk])
                    line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(np.round(alp_f[indt],5))+' , '+str(np.round(bet_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(mag_ob)+' , '+str(Av)+' , '+star_n+' , [0_0_0] , star , '+bostyp+' , '+id_obj_log+' , '+str(flav[indt])+' \n' 
                line1='FIBERMAP '+fib_t+' '+rob_id+' '+bostyp+' '+fib_si+' NA '+str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' '+str(samp)+' '+str(np.int(11959+ran.rand(1)*1000.0))+' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+ ' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 '
                line3=id_obj+' 000000 000000' 
                #line1='ROBOMAPOBJ { 0 0 0 0 0 } '+bostyp+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                #line2=' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+ ' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 0.00000 0.00000 NA '
                #line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' '+id_obj+' '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
            if type == 'SPECTROPHOTO_STD':
                if apo_t[kk] == 0:
                    bostyp='BOSS'
                else:
                    bostyp='APOGEE'
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                #if indt < 30:
                plots=1
                #else:
                #    plots=0
                #id_obj=id_str(fibf[indt],n_z=18)#6)
                id_obj=target_idform(ra_f[indt],dec_f[indt],typ=3)
                if re_obs == 0:
                    Mag=ran.rand(1)*(18.2-16.0)+16.0
                    Mag=Mag[0]
                    Avt=ran.rand(1)*(Av-0.0)+0.0
                    Avt=Avt[0]
                    star_n=star_pool(str='none')
                    Avt=dust_getval(ra_f[indt],dec_f[indt])*3.1
                else:
                    Mag=mago[kk]
                    Avt=Avo[kk]
                    #Avt=dust_getval(ra_f[indt],dec_f[indt])*3.1
                    star_n=obj_id[kk]
                spec_b1,spec_r1,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],dsep=dsep,lshf=lshf,alp=alp_f[indt],bet=bet_f[indt],star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=plots,ifutype=typef1,expt=exptime_sc,ra_0=ra_f[indt],dec_0=dec_f[indt],expof=0.0*ttime_sc+ttime_ar*fex+ttime_fl*fex+toffs,expot=ttime2,toffs=toffs2,t_exp=t_exp)        
                #spec_b2,spec_r2,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                #spec_b3,spec_r3,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                #spec_b4,spec_r4,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)                       
                line1='FIBERMAP '+fib_t+' '+rob_id+' '+bostyp+' '+fib_si+' SPECTROPHOTO_STD '+str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' '+str(samp)+' '+str(np.int(11959+ran.rand(1)*1000.0))+' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 '
                #line3='000000 000000 000000'
                line3=id_obj+' 000000 000000'
                #line1='ROBOMAPOBJ { 3180 301 6 54 65 } '+bostyp+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                #line2=' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+ ' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 0.00000 0.00000 SPECTROPHOTO_STD '
                #line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' 000000 '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
                line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(np.round(alp_f[indt],5))+' , '+str(np.round(bet_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(Mag)+' , '+str(Avt)+' , '+star_n+' , [0_0_0] , star , '+bostyp+' , '+id_obj+' , '+str(flav[indt])+' \n'
            if type == 'SKY':
                if apo_t[kk] == 0:
                    bostyp='BOSS'
                else:
                    bostyp='APOGEE'
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                #id_obj=id_str(fibf[indt],n_z=18)#6)
                id_obj=target_idform(ra_f[indt],dec_f[indt],typ=4)
                if indt < 30:
                    plots=1
                else:
                    plots=0
                spec_b1,spec_r1=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=plots,ifutype=typef1,expt=exptime_sc)
                #spec_b2,spec_r2=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,ifutype=typef1,expt=exptime_sc)
                #spec_b3,spec_r3=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,ifutype=typef1,expt=exptime_sc)
                #spec_b4,spec_r4=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1+"/bhm/",plots=0,ifutype=typef1,expt=exptime_sc)
                line1='FIBERMAP '+fib_t+' '+rob_id+' '+bostyp+' '+fib_si+' SKY '+str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' '+str(samp)+' '+str(np.int(11959+ran.rand(1)*1000.0))+' { 0 0 0 0 0 } 0.00000 '
                #line3='000000 000000 000000'
                line3=id_obj+' 000000 000000' 
                #line1='ROBOMAPOBJ { 0 0 0 0 0 } '+bostyp+' '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                #line2=' { 0.00000 0.00000 0.00000 0.00000 0.00000 } 0.00000 0.00000 0.00000 SKY '
                #line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' 000000 '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
                line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(np.round(alp_f[indt],5))+' , '+str(np.round(bet_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(0)+' , '+str(0)+' , sky , [0_0_0] , sky , '+bostyp+' , '+id_obj+' , '+str(flav[indt])+' \n'
            if kk == 0:
                blu_sp1=np.zeros([len(spec_b1),1])
                red_sp1=np.zeros([len(spec_r1),1])
                #blu_sp2=np.zeros([len(spec_b2),1])
                #red_sp2=np.zeros([len(spec_r2),1])
                #blu_sp3=np.zeros([len(spec_b3),1])
                #red_sp3=np.zeros([len(spec_r3),1])
                #blu_sp4=np.zeros([len(spec_b4),1])
                #red_sp4=np.zeros([len(spec_r4),1])
                blu_sp1[:,0]=spec_b1
                red_sp1[:,0]=spec_r1
                #blu_sp2[:,0]=spec_b2
                #red_sp2[:,0]=spec_r2
                #blu_sp3[:,0]=spec_b3
                #red_sp3[:,0]=spec_r3
                #blu_sp4[:,0]=spec_b4
                #red_sp4[:,0]=spec_r4
            else:
                blu_sp11=np.zeros([len(spec_b1),1])
                red_sp11=np.zeros([len(spec_r1),1])
                #blu_sp12=np.zeros([len(spec_b2),1])
                #red_sp12=np.zeros([len(spec_r2),1])
                #blu_sp13=np.zeros([len(spec_b3),1])
                #red_sp13=np.zeros([len(spec_r3),1])
                #blu_sp14=np.zeros([len(spec_b4),1])
                #red_sp14=np.zeros([len(spec_r4),1])
                blu_sp11[:,0]=spec_b1
                red_sp11[:,0]=spec_r1
                #blu_sp12[:,0]=spec_b2
                #red_sp12[:,0]=spec_r2
                #blu_sp13[:,0]=spec_b3
                #red_sp13[:,0]=spec_r3
                #blu_sp14[:,0]=spec_b4
                #red_sp14[:,0]=spec_r4
                blu_sp1=np.concatenate((blu_sp1,blu_sp11),axis=1)
                red_sp1=np.concatenate((red_sp1,red_sp11),axis=1)
                #blu_sp2=np.concatenate((blu_sp2,blu_sp12),axis=1)
                #red_sp2=np.concatenate((red_sp2,red_sp12),axis=1)
                #blu_sp3=np.concatenate((blu_sp3,blu_sp13),axis=1)
                #red_sp3=np.concatenate((red_sp3,red_sp13),axis=1)
                #blu_sp4=np.concatenate((blu_sp4,blu_sp14),axis=1)
                #red_sp4=np.concatenate((red_sp4,red_sp14),axis=1)
            ft.write(linef)
            fl.write(line_log)
        linef='FIBERMAP -9999 -9999 QUALITY -9999 NA 0.00000 0.00000 '+str(np.round(ra_0,5))+' '+str(np.round(dec_0,5))+' 0 0 { 0.00000 0.00000 0.00000 0.00000 0.00000 } 0.00000 000000 000000 000000\n'          
        #linef='ROBOMAPOBJ { 0 0 0 0 0 } QUALITY '+str(np.round(ra_0,5))+' '+str(np.round(dec_0,5))+' { 0.00000 0.00000 0.00000 0.00000 0.00000 } 0.00000 0.00000 0.00000 NA 0.0000000 0.0000000 0 -9999 -9999 0 0 0\n'
        ft.write(linef)
        ft.close()
        fl.close()
        sys.exit()
#        raw_exp_bhm(blu_sp1,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+2*fex,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=0.0*ttime_sc+ttime_ar*fex+ttime_fl*fex+toffs)
#        raw_exp_bhm(red_sp1,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+2*fex,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=0.0*ttime_sc+ttime_ar*fex+ttime_fl*fex+toffs)
        #raw_exp_bhm(blu_sp2,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+3,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        #raw_exp_bhm(red_sp2,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+3,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        #raw_exp_bhm(blu_sp3,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+4,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        #raw_exp_bhm(red_sp3,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+4,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        #raw_exp_bhm(blu_sp4,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+5,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        #raw_exp_bhm(red_sp4,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+5,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs)

    
def mock_fit_star(dirtemp="libs/",Mag=15.4,lshf=5070.0,reobs_f='none',prs=0,obsn=0,ra_0=180.0,dec_0=0.0,sp_res=2000.0,Av=0.00,sp_samp=1.25,basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1',template2="libs/sky_model.txt",template_1="libs/spEigenStar-55734.fits",base_name='sdR',plate_name='0000',typef1="BOSS",id1='A2-0',psf=0,SN=15.0,Fluxm=20.0,dir1='./',t_exp=0):
    plots=0
    if "BOSS" in typef1:
        if psf <= 0:
            sig=1.43
        else:
            sig=psf
    else:
        if psf == 0:
            sig=1.43
        else:
            sig=psf  
    plate_name=np.str(np.int(np.float(plate_name)+obsn))
    poit_plate_create(Nstd=120,Nsky=80,Ntar=1000,lshf=lshf,ra_0=ra_0,obsn=obsn,dec_0=dec_0,prs=prs,reobs_f=reobs_f,dirtemp=dirtemp,basePath=basePath,Mag=Mag,sp_res=sp_res,Av=Av,sp_samp=sp_samp,template2=template2,template_1=template_1,base_name=base_name,plate_name=plate_name,typef1=typef1,id0=id1,psf=sig,SN=SN,Fluxm=Fluxm,dir1=dir1,t_exp=t_exp)

def head_map_bhm(f,plate,field_id,Nobj,Nsky,Nstd,mjd,obs='APO',ha=0.0,ra=0.0,dec=0.0):
    Ebv=dust_getval(ra, dec)
    Ext=A_l(3.1,np.array([3551.0,4686.0,6165.0,7481.0,8931.0]))*Ebv*3.1
    f1=open('Robo_header.txt','r')
    for line in f1:
        if 'nboss_science' in line:
            line=line.replace('850',str(Nobj))
        if 'nboss_sky' in line:
            line=line.replace('80',str(Nsky))
        if 'nboss_standard' in line:
            line=line.replace('70',str(Nstd))
        if 'configuration_id' in line:
            line=line.replace('0000',str(plate))
        if 'bhmfield_id'in line:
            line=line.replace('0000',str(field_id))
        if 'raCen' in line:
            line=line.replace('213.70400000',str(ra))
        if 'decCen' in line:
            line=line.replace('53.08300000',str(dec))
        if 'haMed' in line:
            line=line.replace('0.000',str(ha))
        if 'observatory' in line:
            line=line.replace('APO',obs)    
        if 'reddeningMed' in line:
            line=line.replace('0.0526',str(np.round(Ext[0],4))).replace('0.0387',str(np.round(Ext[1],4))).replace('0.0281',str(np.round(Ext[2],4))).replace('0.0213',str(np.round(Ext[3],4))).replace('0.0151',str(np.round(Ext[4],4)))
        f.write(line)
            
def head_map_plug(f,plate,Nobj,Nsky,Nstd,mjd,ra=0.0,dec=0.0):
    Ebv=dust_getval(ra, dec)
    Ext=A_l(3.1,np.array([3551.0,4686.0,6165.0,7481.0,8931.0]))*Ebv*3.1
    f1=open('Plug_header.txt','r')
    for line in f1:
        if 'nboss_science' in line:
            line=line.replace('850',str(Nobj))
        if 'nboss_sky' in line:
            line=line.replace('80',str(Nsky))
        if 'nboss_standard' in line:
            line=line.replace('70',str(Nstd))
        if 'mjdDesign' in line:
            line=line.replace('56631',str(mjd))
        if 'plateId' in line:
            line=line.replace('7339',str(plate))
        if 'plateId' in line:
            line=line.replace('7339',str(plate))
        if 'fscanMJD' in line:
            line=line.replace('57518',str(mjd))
        if 'raCen' in line:
            line=line.replace('213.70400000',str(ra))
        if 'decCen' in line:
            line=line.replace('53.08300000',str(dec))
        if 'reddeningMed' in line:
            line=line.replace('0.0526',str(np.round(Ext[0],4))).replace('0.0387',str(np.round(Ext[1],4))).replace('0.0281',str(np.round(Ext[2],4))).replace('0.0213',str(np.round(Ext[3],4))).replace('0.0151',str(np.round(Ext[4],4)))
        f.write(line)

def read_logreobsfile_bhm(file):
    ra_f=[]
    dec_f=[]
    x_f=[]
    y_f=[]
    alp_f=[]
    bet_f=[]
    fib=[]
    typ=[]
    mag=[]
    Av=[]
    obj_id=[]
    obj_cor=[]
    obj_typ=[]
    surv_typ=[]
    flav=[]
    f1=open(file,'r')
    for line in f1:
        if '#' in line:
            data=line.replace('\n','').replace('#','').split(' ')
            data=filter(None,data)
            ra_0=np.float(data[0])
            dec_0=np.float(data[1])
            Ntar=np.int(data[2])
            Nobj=np.int(data[3])
            Nstd=np.int(data[4])
            Nsky=np.int(data[5])
        if not '#' in line:
            if not '$' in line:
                data=line.replace('\n','').split(',')
                data=filter(None,data)
                ra_f.extend([np.float(data[0])])
                dec_f.extend([np.float(data[1])])
                x_f.extend([np.float(data[2])])
                y_f.extend([np.float(data[3])])
                alp_f.extend([np.float(data[4])])
                bet_f.extend([np.float(data[5])])
                fib.extend([np.int(data[6])])
                typ.extend([data[7].replace(' ','')])
                mag.extend([np.float(data[8])])
                Av.extend([np.float(data[9])])
                obj_id.extend([data[10]])
                obj_cor.extend([[np.float(i) for i in data[11].replace('[','').replace(']','').split('_')]])
                obj_typ.extend([data[12].replace(' ','')])
                surv_typ.extend([data[13].replace(' ','')])
                flav.extend([np.int(data[15])])
    ra_f=np.array(ra_f)
    dec_f=np.array(dec_f)
    x_f=np.array(x_f)
    y_f=np.array(y_f)
    alp_f=np.array(alp_f)
    bet_f=np.array(bet_f)
    fib=np.array(fib)
    Av=np.array(Av)
    mag=np.array(mag)
    flav=np.array(flav)
    return ra_f,dec_f,x_f,y_f,alp_f,bet_f,fib,typ,mag,Av,obj_id,obj_cor,obj_typ,ra_0,dec_0,Ntar,Nobj,Nstd,Nsky,surv_typ,flav
        
def read_logreobsfile(file):
    ra_f=[]
    dec_f=[]
    x_f=[]
    y_f=[]
    fib=[]
    typ=[]
    mag=[]
    Av=[]
    obj_id=[]
    obj_cor=[]
    obj_typ=[]
    f1=open(file,'r')
    for line in f1:
        if '#' in line:
            data=line.replace('\n','').replace('#','').split(' ')
            data=filter(None,data)
            ra_0=np.float(data[0])
            dec_0=np.float(data[1])
            Ntar=np.int(data[2])
            Nobj=np.int(data[3])
            Nstd=np.int(data[4])
            Nsky=np.int(data[5])
        if not '#' in line:
            if not '$' in line:
                data=line.replace('\n','').split(',')
                data=filter(None,data)
                ra_f.extend([np.float(data[0])])
                dec_f.extend([np.float(data[1])])
                x_f.extend([np.float(data[2])])
                y_f.extend([np.float(data[3])])
                fib.extend([np.int(data[4])])
                typ.extend([data[5].replace(' ','')])
                mag.extend([np.float(data[6])])
                Av.extend([np.float(data[7])])
                obj_id.extend([data[8]])
                obj_cor.extend([[np.float(i) for i in data[9].replace('[','').replace(']','').split('_')]])
                obj_typ.extend([data[10].replace(' ','')])
    ra_f=np.array(ra_f)
    dec_f=np.array(dec_f)
    x_f=np.array(x_f)
    y_f=np.array(y_f)
    fib=np.array(fib)
    Av=np.array(Av)
    mag=np.array(mag)
    return ra_f,dec_f,x_f,y_f,fib,typ,mag,Av,obj_id,obj_cor,obj_typ,ra_0,dec_0,Ntar,Nobj,Nstd,Nsky

def ecuatorial_galactin(ra,dec):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    c_icrs = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs') 
    l=c_icrs.galactic.l.deg
    b=c_icrs.galactic.b.deg
    return l,b

def dust_getval(ra,dec,not_idl=1):
    l,b=ecuatorial_galactin(ra, dec)
    if not_idl == 0:
        return  0.0
    else:
        call="echo 'print, dust_getval("+str(l)+", "+str(b)+", /interp, ipath="+'"/home/hjibarram/FIT3D_py/idl_files/idlutils/pro/dust/maps/",outfile="dust"'+")' | /usr/local/bin/idl"
        sycall(call)
        if ptt.exists('dust'):
            f=open('dust','r')
            for line in f:
                data=line.replace('\n','').split(' ')
                data=filter(None,data)
                EBV=np.float(data[2])
            sycall('rm dust')
            return EBV
        else:
            return 0.00
        
def poit_plate_create(Nstd=480,reobs_f='none',Nsky=220,lshf=5070.0,ra_0=180.0,obsn=0,dec_0=0.0,rw=1,prs=0,Dp=3.0,Ntar=500,basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1',dirtemp="libs/",Mag=15.4,sp_res=2000.0,Av=0.00,sp_samp=1.25,template2="libs/sky_model.txt",template_1="libs/spEigenStar-55734.fits",template_0="libs/gsd61_156.fits",template_2="libs/templete_gas.fits",template_3="libs/templete_bc03_5.fits",base_name='sdR',plate_name='0000',typef1="BOSS",id0='A2-0',psf=0,SN=15.0,Fluxm=20.0,dir1='./',t_exp=0):
    if not 'none' in reobs_f:
        print reobs_f
        if ptt.exists(reobs_f) == False:
            print 'No reobservation log file'
            re_obs = 0
        else:
            re_obs = 1
    else:
        re_obs =0
    if re_obs == 0:
        Nobj=Ntar-Nstd-Nsky
        Ntarc=Nstd+Nsky
        if "SDSS" in typef1:
            scp_s=60.4#microns per arcsec
        if "BOSS" in typef1:
            scp_s=60.4#microns per arcsec
        if prs == 0:
            obs,ra_0,dec_0,raf,decf,xf,yf,ind0,ho,Om,Lam=mock_sill(fov_p=Dp,rw=rw,ntar=Nobj,dir1=dir1,basePath=basePath,scp_s=scp_s,ra_o=ra_0,dec_o=dec_0)
        else:
            Ntarc=Ntar
        Rp=Dp/2.0
        indf=np.zeros(Ntarc,dtype='int16')
        rr=np.sqrt(ran.rand(Ntarc))*Rp
        rt=ran.rand(Ntarc)*2.0*np.pi
        ra=np.array(rr*np.sin(rt))+ra_0
        dec=np.array(rr*np.cos(rt))+dec_0
        x=(dec-dec_0)*3600.0*scp_s/1e3
        y=(ra-ra_0)*3600.0*scp_s/1e3
        if prs == 0:
            x=np.concatenate((x,xf))
            y=np.concatenate((y,yf))
            ra=np.concatenate((ra,raf))
            dec=np.concatenate((dec,decf))
            indf=np.concatenate((indf,ind0))
        nt2=np.where(y > 0)[0]
        nt1=np.where(y < 0)[0]
        ntp1=np.argsort(x[nt1])
        ntp2=np.argsort(x[nt2])
        fib1=np.array(range(len(nt1)))
        fib2=np.array(range(len(nt1),len(nt1)+len(nt2)))
        fib1f=fib1[ntp1]
        fib2f=fib2[ntp2]
        indf_f1=indf[nt1][ntp1]
        indf_f2=indf[nt2][ntp2]
        x_f1=x[nt1][ntp1]
        y_f1=y[nt1][ntp1]
        x_f2=x[nt2][ntp2]
        y_f2=y[nt2][ntp2]
        dec_f1=dec[nt1][ntp1]
        ra_f1=ra[nt1][ntp1]
        dec_f2=dec[nt2][ntp2]
        ra_f2=ra[nt2][ntp2]
        indf_f=np.concatenate((indf_f1,indf_f2))
        ra_f=np.concatenate((ra_f1,ra_f2))
        dec_f=np.concatenate((dec_f1,dec_f2))
        x_f=np.concatenate((x_f1,x_f2))
        y_f=np.concatenate((y_f1,y_f2))
        fibf=np.concatenate((fib1f,fib2f))+1
        import random as ram    
        ind_std=np.array(ram.sample(range(Ntar), Nstd))  
        ind_tar=np.array(range(Ntar))
        ind_obj1=np.delete(ind_tar,ind_std)
        ind_sky=np.array(ram.sample(ind_obj1, Nsky)) 
    else:
        ra_f,dec_f,x_f,y_f,fibf,typ,mago,Avo,obj_id,obj_cor,obj_typ,ra_0,dec_0,Ntar,Nobj,Nstd,Nsky=read_logreobsfile(dir1+'/'+reobs_f)
        #print obj_typ
        #sys.exit()
    #if len(id1t) == 1:
    #    id=0
    #    while True:
    #        if ptt.exists(dir1+'/'+id1t[0]+'/'+base_name+'-'+id_str(id,n_z=4)+'.spec.fits.gz') == False:
    #            break
    #        else:
    #            id=id+1
    #    id1=id_str(id,n_z=4)
    #else:
    #    id1=id_str(id1t[1],n_z=4)
    exptime_sc=900.0
    exptime_ar=4.0
    exptime_fl=25.0
    cmr=60/exptime_fl
    redot=55.6
    ttime_sc=exptime_sc+redot
    ttime_ar=exptime_ar+redot
    ttime_fl=exptime_fl+redot
    ttime=4.0*ttime_sc+ttime_fl+ttime_ar
    toffs=(obsn-3.5)*ttime
    exp_r=1000+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6+6*(obsn+1)+6*8+6*10+6*9+6*9+6*9+6*9
    cot1=0
    cot2=0
    cot3=0
    if Ntar > 0:
        for kk in range(0,Ntar):
            #print kk
            spec_b,spec_r=mock_arc(dirtemp=dirtemp)
            if kk == 0:
                blu_sp=np.zeros([len(spec_b),1])
                red_sp=np.zeros([len(spec_r),1])
                blu_sp[:,0]=spec_b
                red_sp[:,0]=spec_r
            else:
                blu_sp1=np.zeros([len(spec_b),1])
                red_sp1=np.zeros([len(spec_r),1])
                blu_sp1[:,0]=spec_b
                red_sp1[:,0]=spec_r
                blu_sp=np.concatenate((blu_sp,blu_sp1),axis=1)
                red_sp=np.concatenate((red_sp,red_sp1),axis=1)
        raw_exp(blu_sp,fibf,base_name,fc=[0.88,0.94],n_cr=20,d_cr=2,type="blue",dir1=dir1,mjd=id0,plate=plate_name,flb='a',exp=exp_r+0,ra0=ra_0,dec0=dec_0,expof=toffs)
        raw_exp(red_sp,fibf,base_name,fc=[0.83,0.89],n_cr=20,d_cr=2,type="red",dir1=dir1,mjd=id0,plate=plate_name,flb='a',exp=exp_r+0,ra0=ra_0,dec0=dec_0,expof=toffs)        
        for kk in range(0,Ntar):
            spec_b,spec_r=mock_flat(dirtemp=dirtemp)
            if kk == 0:
                blu_sp=np.zeros([len(spec_b),1])
                red_sp=np.zeros([len(spec_r),1])
                blu_sp[:,0]=spec_b
                red_sp[:,0]=spec_r
            else:
                blu_sp1=np.zeros([len(spec_b),1])
                red_sp1=np.zeros([len(spec_r),1])
                blu_sp1[:,0]=spec_b
                red_sp1[:,0]=spec_r
                blu_sp=np.concatenate((blu_sp,blu_sp1),axis=1)
                red_sp=np.concatenate((red_sp,red_sp1),axis=1)
        raw_exp(blu_sp,fibf,base_name,fc=[0.88,0.94],n_cr=60,type="blue",dir1=dir1,mjd=id0,plate=plate_name,flb='f',exp=exp_r+1,ra0=ra_0,dec0=dec_0,expof=ttime_ar+toffs)
        raw_exp(red_sp,fibf,base_name,fc=[0.83,0.89],n_cr=60,type="red",dir1=dir1,mjd=id0,plate=plate_name,flb='f',exp=exp_r+1,ra0=ra_0,dec0=dec_0,expof=ttime_ar+toffs)        
       # sys.exit()
        dir0=plate_name+"-"+id0        
        sycall('mkdir -p '+dir1+dir0)        
        dir0f=dir0+'/raw_mock'
        sycall('mkdir -p '+dir1+dir0f)        
        ft=open(dir1+dir0f+"/plPlugMapM-"+plate_name+"-"+id0+"-01.par","w")
        fl=open(dir1+"/Log_obs-"+plate_name+"-"+id0+".log","w")
        fl.write('$ RA  DEC  X  Y  FIB  TYPE MAG AV OBJ_ID  OBJ_COR OBJ_TYP \n')
        fl.write('# '+str(ra_0)+' '+str(dec_0)+' '+str(Ntar)+' '+str(Nobj)+' '+str(Nstd)+' '+str(Nsky)+'\n')
        head_map_plug(ft,plate_name,Nobj,Nsky,Nstd,id0,ra=ra_0,dec=dec_0)
        for kk in range(0,Ntar):
            if re_obs == 0:
                kt=0
                for jj in range(0, Nstd):  
                    if ind_tar[kk] == ind_std[jj]:
                        kt=1
                        indt=ind_tar[kk]
                        type='SPECTROPHOTO_STD'
                        cot1=cot1+1
                for jj in range(0, Nsky):  
                    if ind_tar[kk] == ind_sky[jj]:
                        kt=2
                        indt=ind_tar[kk]
                        type='SKY'
                        cot2=cot2+1
                if kt == 0:
                    indt=ind_tar[kk]
                    type='NA'
                    cot3=cot3+1
            else:
                type=typ[kk]
                indt=kk
                #print type
            if fibf[indt] >= 500:
                samp=2
            else:
                samp=1
            id1=id0+'-'+id_str(fibf[indt],n_z=4)
            if fibf[indt] % 20 == 1:
                t=ran.rand(1)[0]*100.0
                if t  < 5:
                    fib_t="-1"
                else:
                    fib_t=str(fibf[indt])
            else:
                fib_t=str(fibf[indt])
            if type == 'NA':
                #print re_obs
                if re_obs == 0:
                    idh=indf_f[indt]
                    star_n=star_pool(str='none')
                    Av=dust_getval(ra_f[indt],dec_f[indt])*3.1
                else:
                    Av=Avo[kk]
                    #Av=dust_getval(ra_f[indt],dec_f[indt])*3.1
                    #print obj_typ[kk]
                    if obj_typ[kk] == 'galaxy':
                        idh=np.int(obj_id[kk])
                        prs=0
                        obs=obj_cor[kk]
                    if obj_typ[kk] == 'star':
                        prs=1
                        star_n=obj_id[kk]
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                if indt < 20:
                    plots=1
                else:
                    plots=0
                if prs == 0:
                    spec_b1,spec_r1,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,plots=plots,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=0.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    spec_b2,spec_r2,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    spec_b3,spec_r3,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    spec_b4,spec_r4,mag_v=mock_gali(idh,plate_name+"-"+id1t[0],id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,mjd=id1t[0],Av_gal=Av,dirtemp=dirtemp,template4=template2,template3=template_0,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,plots=0,outs=0,observer=obs,ifutype=typef1,basePath=basePath,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)
                    line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(-1)+' , '+str(Av)+' , '+idh+' , ['+str(obs[0])+'_'+str(obs[1])+'_'+str(obs[2])+'] , galaxy\n'
                else: 
                    if re_obs == 0:  
                        mag_ob=ran.rand(1)*(22.2-17.0)+17.0
                        mag_ob=mag_ob[0]
                    else:
                        mag_ob=mago[kk]
                    #print star_n
                    spec_b1,spec_r1,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=plots,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=0.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                    spec_b2,spec_r2,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                    spec_b3,spec_r3,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                    spec_b4,spec_r4,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=mag_ob,dirtemp=dirtemp,sp_res=sp_res,Av_g=Av,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)    
                    line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(mag_ob)+' , '+str(Av)+' , '+star_n+' , [0_0_0] , star\n' 
                line1='PLUGMAPOBJ { 0 0 0 0 0 } OBJECT '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+ ' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 0.00000 0.00000 NA '
                line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
            if type == 'SPECTROPHOTO_STD':
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                #if indt < 30:
                plots=1
                #else:
                #    plots=0
                if re_obs == 0:
                    Mag=ran.rand(1)*(18.2-16.0)+16.0
                    Mag=Mag[0]
                    Avt=ran.rand(1)*(Av-0.0)+0.0
                    Avt=Avt[0]
                    star_n=star_pool(str='none')
                    Avt=dust_getval(ra_f[indt],dec_f[indt])*3.1
                else:
                    Mag=mago[kk]
                    Avt=Avo[kk]
                    #Avt=dust_getval(ra_f[indt],dec_f[indt])*3.1
                    star_n=obj_id[kk]
                spec_b1,spec_r1,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=plots,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=0.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                spec_b2,spec_r2,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                spec_b3,spec_r3,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)        
                spec_b4,spec_r4,mag_v=mock_star(id1,x_f[indt],y_f[indt],fibf[indt],lshf=lshf,star_t=star_n,mjd=id1t[0],idn=plate_name+"-"+id1t[0],Mag_s=Mag,dirtemp=dirtemp,sp_res=sp_res,Av_g=Avt,basename=base_name+'-',template2=template2,template=template_1,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,plots=0,outs=0,ifutype=typef1,expt=exptime_sc,ra_0=ra_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs,expot=ttime,toffs=toffs,t_exp=t_exp)                       
                line1='PLUGMAPOBJ { 3180 301 6 54 65 } OBJECT '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' { '+str(np.round(mag_v[0],4))+' '+str(np.round(mag_v[1],4))+ ' '+str(np.round(mag_v[2],4))+' '+str(np.round(mag_v[3],4))+' '+str(np.round(mag_v[4],4))+' } 0.00000 0.00000 0.00000 SPECTROPHOTO_STD '
                line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
                line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(Mag)+' , '+str(Avt)+' , '+star_n+' , [0_0_0] , star\n'
            if type == 'SKY':
                id1t=id1.split('-')
                id1=id_str(id1t[1],n_z=4)                
                name=base_name+'-'+id1
                if indt < 30:
                    plots=1
                else:
                    plots=0
                spec_b1,spec_r1=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1,plots=plots,ifutype=typef1,expt=exptime_sc)
                spec_b2,spec_r2=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1,plots=0,ifutype=typef1,expt=exptime_sc)
                spec_b3,spec_r3=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1,plots=0,ifutype=typef1,expt=exptime_sc)
                spec_b4,spec_r4=mock_sky(id1,x_f[indt],y_f[indt],fibf[indt],idn=plate_name+"-"+id1t[0],dirtemp=dirtemp,sp_res=sp_res,basename=base_name+'-',template2=template2,SN=SN,Fluxm=Fluxm,dir1=dir1,plots=0,ifutype=typef1,expt=exptime_sc)
                line1='PLUGMAPOBJ { 0 0 0 0 0 } OBJECT '+str(np.round(ra_f[indt],5))+' '+str(np.round(dec_f[indt],5))
                line2=' { 0.00000 0.00000 0.00000 0.00000 0.00000 } 0.00000 0.00000 0.00000 SKY '
                line3=str(np.round(x_f[indt],5))+' '+str(np.round(y_f[indt],5))+' '+str(samp)+' '+fib_t+' '+str(np.int(11959+ran.rand(1)*1000.0))+' 0 0' 
                linef=line1+line2+line3+'\n'
                line_log=str(np.round(ra_f[indt],5))+' , '+str(np.round(dec_f[indt],5))+' , '+str(np.round(x_f[indt],5))+' , '+str(np.round(y_f[indt],5))+' , '+str(fibf[indt])+' , '+type+' , '+str(0)+' , '+str(0)+' , sky , [0_0_0] , sky\n'
            if kk == 0:
                blu_sp1=np.zeros([len(spec_b1),1])
                red_sp1=np.zeros([len(spec_r1),1])
                blu_sp2=np.zeros([len(spec_b2),1])
                red_sp2=np.zeros([len(spec_r2),1])
                blu_sp3=np.zeros([len(spec_b3),1])
                red_sp3=np.zeros([len(spec_r3),1])
                blu_sp4=np.zeros([len(spec_b4),1])
                red_sp4=np.zeros([len(spec_r4),1])
                blu_sp1[:,0]=spec_b1
                red_sp1[:,0]=spec_r1
                blu_sp2[:,0]=spec_b2
                red_sp2[:,0]=spec_r2
                blu_sp3[:,0]=spec_b3
                red_sp3[:,0]=spec_r3
                blu_sp4[:,0]=spec_b4
                red_sp4[:,0]=spec_r4
            else:
                blu_sp11=np.zeros([len(spec_b1),1])
                red_sp11=np.zeros([len(spec_r1),1])
                blu_sp12=np.zeros([len(spec_b2),1])
                red_sp12=np.zeros([len(spec_r2),1])
                blu_sp13=np.zeros([len(spec_b3),1])
                red_sp13=np.zeros([len(spec_r3),1])
                blu_sp14=np.zeros([len(spec_b4),1])
                red_sp14=np.zeros([len(spec_r4),1])
                blu_sp11[:,0]=spec_b1
                red_sp11[:,0]=spec_r1
                blu_sp12[:,0]=spec_b2
                red_sp12[:,0]=spec_r2
                blu_sp13[:,0]=spec_b3
                red_sp13[:,0]=spec_r3
                blu_sp14[:,0]=spec_b4
                red_sp14[:,0]=spec_r4
                blu_sp1=np.concatenate((blu_sp1,blu_sp11),axis=1)
                red_sp1=np.concatenate((red_sp1,red_sp11),axis=1)
                blu_sp2=np.concatenate((blu_sp2,blu_sp12),axis=1)
                red_sp2=np.concatenate((red_sp2,red_sp12),axis=1)
                blu_sp3=np.concatenate((blu_sp3,blu_sp13),axis=1)
                red_sp3=np.concatenate((red_sp3,red_sp13),axis=1)
                blu_sp4=np.concatenate((blu_sp4,blu_sp14),axis=1)
                red_sp4=np.concatenate((red_sp4,red_sp14),axis=1)
            ft.write(linef)
            fl.write(line_log)
        linef='PLUGMAPOBJ { 0 0 0 0 0 } QUALITY '+str(np.round(ra_0,5))+' '+str(np.round(dec_0,5))+' { 0.00000 0.00000 0.00000 0.00000 0.00000 } 0.00000 0.00000 0.00000 NA 0.0000000 0.0000000 0 -9999 0 0 0\n'
        ft.write(linef)
        ft.close()
        fl.close()
        #ran_b=ran.randn(24)*0.7
        #ran_r=ran.randn(24)*0.7  fc=[0.50,0.35] BLUE
        raw_exp(blu_sp1,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+2,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=0.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(red_sp1,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+2,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=0.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(blu_sp2,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+3,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(red_sp2,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+3,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=1.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(blu_sp3,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+4,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(red_sp3,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+4,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=2.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(blu_sp4,fibf,base_name,fc=[0.88,0.94],n_cr=600,type="blue",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+5,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs)
        raw_exp(red_sp4,fibf,base_name,fc=[0.83,0.89],n_cr=600,type="red",dir1=dir1,mjd=id1t[0],plate=plate_name,flb='s',exp=exp_r+5,expt=exptime_sc,ra0=ra_0,dec0=dec_0,expof=3.0*ttime_sc+ttime_ar+ttime_fl+toffs)
            
def cosmic_rays(arrayf,n_cr=100,d_cr=5):
    nf,ng=arrayf.shape
    array1=np.zeros([nf,ng])    
    nc=n_cr+np.int(ran.randn(1)[0]*d_cr) 
    xof=ran.rand(nc)*nf
    yof=ran.rand(nc)*ng
    thet=ran.rand(nc)*(90.0+90.0)-90.0
    phi=ran.rand(nc)*(90.0-5.0)+5.0
    deep=10.0
    lent=deep/np.sin(phi*np.pi/180.0)
    for k in range(0, nc):
        lx=np.int(lent[k]*np.cos(thet[k]*np.pi/180.0))+1        
        x_tc=np.arange(lx)+xof[k]
        cof=yof[k]-np.tan(thet[k]*np.pi/180.0)*xof[k]       
        y_tc=np.tan(thet[k]*np.pi/180.0)*x_tc+cof
        for i in range(0, lx):
            xt1=np.int(x_tc[i])-1
            xt2=np.int(x_tc[i])#+1
            yt1=np.int(y_tc[i])-1
            yt2=np.int(y_tc[i])#+1
            if yt1 < 0:
                yt1=0
            if yt2 < 0:
                yt2 =0
            if yt1 > ng:
                yt1=ng
            if yt2 > ng:
                yt2=ng
            if xt1 < 0:
                xt1=0
            if xt2 < 0:
                xt2=0
            if xt1 > nf:
                xt1=nf
            if xt2 > nf:
                xt2=nf
            array1[xt1:xt2,yt1:yt2]=100000.0#60000.0
    dv=0.5
    PSF=Gaussian2DKernel(stddev=dv)
    array1=convolve(array1, PSF)#, mode='full')#,  boundary='symm')#photo_a#
    arrayf=arrayf+array1
    arrayf[np.where(arrayf >= 60000.0)]=60000.0
    return arrayf

def read_op_fib(cart,cam,dir='libs/'):
    if cart < 1 and cart > 18:
        cart=16
    if not 'b1' in cam:
        if not 'b2' in cam:
            if not 'r1' in cam:
                if not 'r2' in cam:
                    cam='b1' 
    file=dir+'opFibers.par'  
    f=open(file,'r')
    for line in f:
        if 'FIBERPARAM' in line:
            dat=line.replace('\n','').split('{')
            data1=dat[0].split(' ')
            data1=filter(None,data1)
            if len(data1) > 2:
                car=np.int(data1[1])
                lap=data1[2]
                if car == cart and cam == lap:
                    print car,lap
                    data2=dat[1].replace('}','').split(' ')
                    data2=filter(None,data2)
                    space=np.array([np.float(val) for val in data2])
                    data3=dat[2].replace('}','').split(' ')
                    data3=filter(None,data3)
                    bspa=np.array([np.float(val) for val in data3])
                    bspa[0]=bspa[0]-279.0
    return space,bspa


def raw_exp_bhm(spec,fibf,base_name,fc=[1.0,1.0],n_cr=130,d_cr=5,type="blue",dir1='./',l=500,mjd='00000',plate='0000',exp=0,flb='s',expt=900.0,ra0=0.0,dec0=0.0,expof=0.0):       
    nx,ny=spec.shape
    nf=4224
    ng=4352
    arrayf_1=np.zeros([nf,ng])
    b_arrayf_1=np.zeros([nf,ng])
    if type == "blue":
        p_arrayf_1=np.zeros([4112,4096])
        f_arrayf_1=np.ones([4112,4096])
        b1_arrayf_1=np.zeros([4112,4096])
    if type == "red":
        p_arrayf_1=np.zeros([4128,4114])    
        f_arrayf_1=np.ones([4128,4114])
        b1_arrayf_1=np.zeros([4128,4114])
    nt=np.argsort(fibf)
    if type == "blue":
        let=800
        let2=300
        ty='b'
        fibs1,bunds1=read_op_fib(16,'b1')
    if type == "red":
        let=620
        let2=300
        ty='r'
        fibs1,bunds1=read_op_fib(16,'r1')
    for i in range(0, len(fibf)):
        if i < 500:
            it=i
        else:
            it=i-500
        r=let*(1.3+0.0)
        then=np.arcsin((it-250)/r)                 
        dr=np.cos(then)*r-let*(0.3+0.0)
        dx=np.int(np.round(dr))
        spect=spec[:,nt[i]]
        npix1=100
        npix2=200
        y_e=spect[0]/np.float(npix1)*np.arange(npix1)
        y_r=spect[::-1][0]/np.float(npix2)*np.arange(npix2)
        spect=np.concatenate((y_e,spect))
        spect=np.concatenate((spect,y_r[::-1]))
        dx=dx-npix1        
        nx=len(spect)
        dsi=np.ones(nx)
        indt=it % 20
        itt=it/20
        if indt == 0:
            if it > 0:
                dt=16.4541#+np.abs(ran.randn(1)*0.7)-0.1
                dt=bunds1[itt]+fibs1[itt]
            else:
                if type == 'blue':
                    dt=403.0+bunds1[0]+36.6-2.0#+6.578*4.5#+280
                if type == 'red':
                    dt=403.0+bunds1[0]+36.6-18.3+9.15
                dg=0.0
        else:
            dt=6.578
            dt=fibs1[itt]
        dg=dg+dt

        dy=np.int(np.round(dg))
        off=dy-dg
        dc=10.0
        nxt=np.int(dc*2+1)
        x_t=np.arange(nxt)*1.0-dc+off
        x_t=np.array([x_t]*nx)
        ds_t=np.array([dsi]*nxt).T
        At=np.array([spect]*nxt).T
        spec_t=np.exp(-0.5*(x_t/ds_t)**2.0)/(np.sqrt(2.0*np.pi)*ds_t)*At
        y1=0
        y2=nxt
        x1=0
        x2=nx
        arrayf_1[dx+x1:dx+x2,dy+y1:dy+y2]=spec_t+arrayf_1[dx+x1:dx+x2,dy+y1:dy+y2]    
    
    if type == "blue":
        bias_1=2171.0
        sig_1=2.0*fc[0]*0.56
        gain_1=[1.048, 1.048, 1.018, 1.006]
        arrayf_1[0:2111,0:2175]=arrayf_1[0:2111,0:2175]/gain_1[0]+ran.randn(2111,2175)*sig_1+bias_1+35.0#2206
        arrayf_1[2111:4224,0:2175]=arrayf_1[2111:4224,0:2175]/gain_1[1]+ran.randn(2113,2175)*sig_1+bias_1+0.0#2171       
        arrayf_1[0:2111,2175:4352]=arrayf_1[0:2111,2175:4352]/gain_1[2]+ran.randn(2111,2177)*sig_1+bias_1-30.0#2141
        arrayf_1[2111:4224,2175:4352]=arrayf_1[2111:4224,2175:4352]/gain_1[3]+ran.randn(2113,2177)*sig_1+bias_1+129.0#2300
        b_arrayf_1[0:2111,0:2175]=ran.randn(2111,2175)*sig_1+bias_1+35.0#2206
        b_arrayf_1[2111:4224,0:2175]=ran.randn(2113,2175)*sig_1+bias_1+0.0#2171       
        b_arrayf_1[0:2111,2175:4352]=ran.randn(2111,2177)*sig_1+bias_1-30.0#2141
        b_arrayf_1[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_1+bias_1+129.0#2300
        b1_arrayf_1=ran.randn(4112,4096)*sig_1+2.0
    if type == "red":
        bias_1=2120.0#2490.0
        sig_1=2.0*fc[0]*0.56
        gain_1=[1.9253, 1.5122, 1.4738, 1.5053]
        arrayf_1[0:2111,0:2175]=arrayf_1[0:2111,0:2175]/gain_1[0]+ran.randn(2111,2175)*sig_1+bias_1+5.0#+2495-4.0
        arrayf_1[2111:4224,0:2175]=arrayf_1[2111:4224,0:2175]/gain_1[1]+ran.randn(2113,2175)*sig_1+bias_1-4.0#2490-4.0
        arrayf_1[0:2111,2175:4352]=arrayf_1[0:2111,2175:4352]/gain_1[2]+ran.randn(2111,2177)*sig_1+bias_1+55.0#+2545-4.0
        arrayf_1[2111:4224,2175:4352]=arrayf_1[2111:4224,2175:4352]/gain_1[3]+ran.randn(2113,2177)*sig_1+bias_1+100.0#246.0#2740-4.0
        b_arrayf_1[0:2111,0:2175]=ran.randn(2111,2175)*sig_1+bias_1+5.0#+2495-4.0
        b_arrayf_1[2111:4224,0:2175]=ran.randn(2113,2175)*sig_1+bias_1-4.0#2490-4.0
        b_arrayf_1[0:2111,2175:4352]=ran.randn(2111,2177)*sig_1+bias_1+55.0#+2545-4.0
        b_arrayf_1[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_1+bias_1+100.0#246.0#2740-4.0
        b1_arrayf_1=ran.randn(4128,4114)*sig_1+2.0
    
    sycall('mkdir -p '+dir1+'bhm')        
    dir0='bhm/'+str(mjd)#+'/'+plate.split('-')[0]#+'-'+mjd        
    sycall('mkdir -p '+dir1+'bhm/'+str(mjd))
    #sycall('mkdir -p '+dir1+dir0)        
    dir0f=dir0+'/raw_mock'
    sycall('mkdir -p '+dir1+dir0f)  
    dir0f=dir0f+'/'
    arrayf_1=cosmic_rays(arrayf_1,n_cr=n_cr,d_cr=d_cr)-32768.0
    arrayf_1[np.where(arrayf_1 > 32767)]=32767.0
    arrayf_1=np.array(arrayf_1,dtype='int16')
    h1=pyf.PrimaryHDU(arrayf_1)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=ng
    h["NAXIS2"]=nf
    h=row_data_header_bhm(h,plate,mjd,exp,ty+'',flb=flb,expt=expt,ra=ra0,dec=dec0,expof=expof)
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir1+dir0f+base_name+'-'+ty+'1-'+id_str(exp,n_z=8)+'.fit'
    wfits_ext(out_fit,hlist)
    sycall('gzip -f '+out_fit)


    
    if flb == 'f':
        sycall('mkdir -p '+dir1+'bhm/flats')
        sycall('mkdir -p '+dir1+'bhm/biases')
        h1=pyf.PrimaryHDU(b_arrayf_1)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=ng
        h["NAXIS2"]=nf
        h=row_data_header2(h,mjd)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'bhm/biases/boss_pixbias-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
    

        h1=pyf.PrimaryHDU(f_arrayf_1)
        h=h1.header
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'bhm/flats/pixflatave-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit) 
        
        h1=pyf.PrimaryHDU(p_arrayf_1)
        h=h1.header
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'bhm/flats/badpixels-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
    
        
        h1=pyf.PrimaryHDU(b1_arrayf_1)
        h=h1.header
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'bhm/biases/pixbiasave-00001'+'-'+ty+'1.fits'
        wfits_ext(out_fit,hlist)
        #sycall('gzip -f '+out_fit)
    
        
        
def raw_exp(spec,fibf,base_name,fc=[1.0,1.0],n_cr=130,d_cr=5,type="blue",dir1='./',l=500,mjd='00000',plate='0000',exp=0,flb='s',expt=900.0,ra0=0.0,dec0=0.0,expof=0.0):       
    nx,ny=spec.shape
    nf=4224
    ng=4352
    arrayf_1=np.zeros([nf,ng])
    arrayf_2=np.zeros([nf,ng])
    b_arrayf_1=np.zeros([nf,ng])
    b_arrayf_2=np.zeros([nf,ng])
    if type == "blue":
        p_arrayf_1=np.zeros([4112,4096])
        p_arrayf_2=np.zeros([4112,4096])
        f_arrayf_1=np.ones([4112,4096])
        f_arrayf_2=np.ones([4112,4096])
        b1_arrayf_1=np.zeros([4112,4096])
        b1_arrayf_2=np.zeros([4112,4096])
    if type == "red":
        p_arrayf_1=np.zeros([4128,4114])
        p_arrayf_2=np.zeros([4128,4114])    
        f_arrayf_1=np.ones([4128,4114])
        f_arrayf_2=np.ones([4128,4114])
        b1_arrayf_1=np.zeros([4128,4114])
        b1_arrayf_2=np.zeros([4128,4114])
    nt=np.argsort(fibf)
    #*0.6
    if type == "blue":
        let=800
        let2=300
        ty='b'
        fibs1,bunds1=read_op_fib(16,'b1')
        fibs2,bunds2=read_op_fib(16,'b2')
    if type == "red":
        let=620
        let2=300
        ty='r'
        fibs1,bunds1=read_op_fib(16,'r1')
        fibs2,bunds2=read_op_fib(16,'r2')
    for i in range(0, len(fibf)):
        if i < 500:
            it=i
        else:
            it=i-500
        r=let*(1.3+0.0)
        then=np.arcsin((it-250)/r)                 
        dr=np.cos(then)*r-let*(0.3+0.0)
        dx=np.int(np.round(dr))
        spect=spec[:,nt[i]]
        npix1=100
        npix2=200
        y_e=spect[0]/np.float(npix1)*np.arange(npix1)
        y_r=spect[::-1][0]/np.float(npix2)*np.arange(npix2)
        spect=np.concatenate((y_e,spect))
        spect=np.concatenate((spect,y_r[::-1]))
        #nxf=nf-dx-let2#-npix1#-npix2
        #nxi=len(spect)+npix2
        #fact=np.float(nxf)/np.float(nxi)
        #npix1f=np.int(np.round(npix1*fact))
        #npix2f=np.int(np.round(npix2*fact))
        #nlen=np.int(np.round(len(spect)*fact))
        #nxf=nlen+npix1f+npix2f
        #dx=dx-npix1f
        dx=dx-npix1
        #xi=np.linspace(0,nxf-1,len(spect))  
        #xf=np.arange(nxf)      
        #spect=interp1d(xi,spect,bounds_error=False,fill_value=0.)(xf)
        
        
        nx=len(spect)
        dsi=np.ones(nx)
        #con=0
        if i < 500:
            indt=it % 20
            itt=it/20
            if indt == 0:
                if it > 0:
                    #if len(ran_bs) >= 24:
                    #    dt=16.4541+np.abs(ran_bs[con])-0.1
                    #else:
                    dt=16.4541#+np.abs(ran.randn(1)*0.7)-0.1
                    dt=bunds1[itt]+fibs1[itt]
                    #con=con+1
                    #if dt < 16.4:
                    #    dt = 16.4
                    #print dt
                else:
                    if type == 'blue':
                        dt=403.0+bunds1[0]+36.6-2.0#+6.578*4.5#+280
                    if type == 'red':
                        dt=403.0+bunds1[0]+36.6-18.3+9.15
                    dg=0.0
            else:
                dt=6.578
                dt=fibs1[itt]
            dg=dg+dt
        else:
            indt=it % 20
            itt=it/20
            if indt == 0:
                if it > 0:
                    dt=16.4541
                    dt=bunds2[itt]+fibs2[itt]
                    #print dt
                else:
                    if type == 'blue':
                        dt=403.0+bunds2[0]+36.6-18.3+9.15-2.0
                    if type == 'red':
                        dt=403.0+bunds2[0]+36.6-18.3-18.3+18.3+9.15
                    dg=0.0
            else:
                dt=6.578
                dt=fibs2[itt]
            dg=dg+dt
        dy=np.int(np.round(dg))
        off=dy-dg
        dc=10.0
        nxt=np.int(dc*2+1)
        x_t=np.arange(nxt)*1.0-dc+off
        x_t=np.array([x_t]*nx)
        ds_t=np.array([dsi]*nxt).T
        At=np.array([spect]*nxt).T
        spec_t=np.exp(-0.5*(x_t/ds_t)**2.0)/(np.sqrt(2.0*np.pi)*ds_t)*At
        y1=0
        y2=nxt
        x1=0
        x2=nx
        if i < 500:
            arrayf_1[dx+x1:dx+x2,dy+y1:dy+y2]=spec_t+arrayf_1[dx+x1:dx+x2,dy+y1:dy+y2]
        else:
            arrayf_2[dx+x1:dx+x2,dy+y1:dy+y2]=spec_t+arrayf_2[dx+x1:dx+x2,dy+y1:dy+y2]
    
    
    if type == "blue":
        bias_1=2171.0
        bias_2=2171.0
        sig_1=2.0*fc[0]*0.56
        sig_2=2.0*fc[1]*0.56
        gain_1=[1.048, 1.048, 1.018, 1.006]
        gain_2=[1.040, 0.994, 1.002, 1.010]
        arrayf_1[0:2111,0:2175]=arrayf_1[0:2111,0:2175]/gain_1[0]+ran.randn(2111,2175)*sig_1+bias_1+35.0#2206
        arrayf_1[2111:4224,0:2175]=arrayf_1[2111:4224,0:2175]/gain_1[1]+ran.randn(2113,2175)*sig_1+bias_1+0.0#2171       
        arrayf_1[0:2111,2175:4352]=arrayf_1[0:2111,2175:4352]/gain_1[2]+ran.randn(2111,2177)*sig_1+bias_1-30.0#2141
        arrayf_1[2111:4224,2175:4352]=arrayf_1[2111:4224,2175:4352]/gain_1[3]+ran.randn(2113,2177)*sig_1+bias_1+129.0#2300
        arrayf_2[0:2111,0:2175]=arrayf_2[0:2111,0:2175]/gain_2[0]+ran.randn(2111,2175)*sig_2+bias_2+60.0#2170-5.0
        arrayf_2[2111:4224,0:2175]=arrayf_2[2111:4224,0:2175]/gain_2[1]+ran.randn(2113,2175)*sig_2+bias_2-5.0#2110-5.0
        arrayf_2[0:2111,2175:4352]=arrayf_2[0:2111,2175:4352]/gain_2[2]+ran.randn(2111,2177)*sig_2+bias_2+42.0#2150-5.0
        arrayf_2[2111:4224,2175:4352]=arrayf_2[2111:4224,2175:4352]/gain_2[3]+ran.randn(2113,2177)*sig_2+bias_2+11.0#2121-5.0
        b_arrayf_1[0:2111,0:2175]=ran.randn(2111,2175)*sig_1+bias_1+35.0#2206
        b_arrayf_1[2111:4224,0:2175]=ran.randn(2113,2175)*sig_1+bias_1+0.0#2171       
        b_arrayf_1[0:2111,2175:4352]=ran.randn(2111,2177)*sig_1+bias_1-30.0#2141
        b_arrayf_1[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_1+bias_1+129.0#2300
        b_arrayf_2[0:2111,0:2175]=ran.randn(2111,2175)*sig_2+bias_2+60.0#2170-5.0
        b_arrayf_2[2111:4224,0:2175]=ran.randn(2113,2175)*sig_2+bias_2-5.0#2110-5.0
        b_arrayf_2[0:2111,2175:4352]=ran.randn(2111,2177)*sig_2+bias_2+42.0#2150-5.0
        b_arrayf_2[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_2+bias_2+11.0#2121-5.0
        #b_arrayf_1=ran.randn(nf,ng)*sig_1+bias_1
        #b_arrayf_2=ran.randn(nf,ng)*sig_2+bias_2
        b1_arrayf_1=ran.randn(4112,4096)*sig_1+2.0
        b1_arrayf_2=ran.randn(4112,4096)*sig_2+2.0
    if type == "red":
        bias_1=2120.0#2490.0
        bias_2=2112.0#2670.0
        sig_1=2.0*fc[0]*0.56
        sig_2=2.0*fc[1]*0.56
        gain_1=[1.9253, 1.5122, 1.4738, 1.5053]
        gain_2=[1.5980, 1.6560, 1.5820, 1.5940]
        arrayf_1[0:2111,0:2175]=arrayf_1[0:2111,0:2175]/gain_1[0]+ran.randn(2111,2175)*sig_1+bias_1+5.0#+2495-4.0
        arrayf_1[2111:4224,0:2175]=arrayf_1[2111:4224,0:2175]/gain_1[1]+ran.randn(2113,2175)*sig_1+bias_1-4.0#2490-4.0
        arrayf_1[0:2111,2175:4352]=arrayf_1[0:2111,2175:4352]/gain_1[2]+ran.randn(2111,2177)*sig_1+bias_1+55.0#+2545-4.0
        arrayf_1[2111:4224,2175:4352]=arrayf_1[2111:4224,2175:4352]/gain_1[3]+ran.randn(2113,2177)*sig_1+bias_1+100.0#246.0#2740-4.0
        arrayf_2[0:2111,0:2175]=arrayf_2[0:2111,0:2175]/gain_2[0]+ran.randn(2111,2175)*sig_2+bias_2-75.0#+2595-8.0
        arrayf_2[2111:4224,0:2175]=arrayf_2[2111:4224,0:2175]/gain_2[1]+ran.randn(2113,2175)*sig_2+bias_2-8.0#+2670-8.0
        arrayf_2[0:2111,2175:4352]=arrayf_2[0:2111,2175:4352]/gain_2[2]+ran.randn(2111,2177)*sig_2+bias_2-100.0#245.0#+2425-8.0
        arrayf_2[2111:4224,2175:4352]=arrayf_2[2111:4224,2175:4352]/gain_2[3]+ran.randn(2113,2177)*sig_2+bias_2-100.0#235.0#+2435-8.0
        b_arrayf_1[0:2111,0:2175]=ran.randn(2111,2175)*sig_1+bias_1+5.0#+2495-4.0
        b_arrayf_1[2111:4224,0:2175]=ran.randn(2113,2175)*sig_1+bias_1-4.0#2490-4.0
        b_arrayf_1[0:2111,2175:4352]=ran.randn(2111,2177)*sig_1+bias_1+55.0#+2545-4.0
        b_arrayf_1[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_1+bias_1+100.0#246.0#2740-4.0
        b_arrayf_2[0:2111,0:2175]=ran.randn(2111,2175)*sig_2+bias_2-75.0#+2595-8.0
        b_arrayf_2[2111:4224,0:2175]=ran.randn(2113,2175)*sig_2+bias_2-8.0#+2670-8.0
        b_arrayf_2[0:2111,2175:4352]=ran.randn(2111,2177)*sig_2+bias_2-100.0#245.0#+2425-8.0
        b_arrayf_2[2111:4224,2175:4352]=ran.randn(2113,2177)*sig_2+bias_2-100.0#235.0#+2435-8.0
        #b_arrayf_1=ran.randn(nf,ng)*sig_1+bias_1
        #b_arrayf_2=ran.randn(nf,ng)*sig_2+bias_2
        b1_arrayf_1=ran.randn(4128,4114)*sig_1+2.0
        b1_arrayf_2=ran.randn(4128,4114)*sig_2+2.0
           
    dir0=plate+'-'+mjd        
    sycall('mkdir -p '+dir1+dir0)        
    dir0f=dir0+'/raw_mock'
    sycall('mkdir -p '+dir1+dir0f)  
    dir0f=dir0f+'/'
    arrayf_1=cosmic_rays(arrayf_1,n_cr=n_cr,d_cr=d_cr)-32768.0
    arrayf_2=cosmic_rays(arrayf_2,n_cr=n_cr,d_cr=d_cr)-32768.0
    arrayf_1[np.where(arrayf_1 > 32767)]=32767.0
    arrayf_2[np.where(arrayf_2 > 32767)]=32767.0
    arrayf_1=np.array(arrayf_1,dtype='int16')
    arrayf_2=np.array(arrayf_2,dtype='int16')
    h1=pyf.PrimaryHDU(arrayf_1)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=ng
    h["NAXIS2"]=nf
    h=row_data_header(h,plate,mjd,exp,ty+'1',flb=flb,expt=expt,ra=ra0,dec=dec0,expof=expof)
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir1+dir0f+base_name+'-'+ty+'1-'+id_str(exp,n_z=8)+'.fit'
    wfits_ext(out_fit,hlist)
    sycall('gzip -f '+out_fit)

    h1=pyf.PrimaryHDU(arrayf_2)
    h=h1.header
    h["NAXIS"]=2 
    h["NAXIS1"]=ng
    h["NAXIS2"]=nf
    h=row_data_header(h,plate,mjd,exp,ty+'2',flb=flb,expt=expt,ra=ra0,dec=dec0,expof=expof)
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir1+dir0f+base_name+'-'+ty+'2-'+id_str(exp,n_z=8)+'.fit'
    wfits_ext(out_fit,hlist)
    sycall('gzip -f '+out_fit)
    
    if flb == 'f':
        sycall('mkdir -p '+dir1+'flats')
        sycall('mkdir -p '+dir1+'biases')
        h1=pyf.PrimaryHDU(b_arrayf_1)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=ng
        h["NAXIS2"]=nf
        h=row_data_header2(h,mjd)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'biases/boss_pixbias-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
    
        h1=pyf.PrimaryHDU(b_arrayf_2)
        h=h1.header
        h["NAXIS"]=2 
        h["NAXIS1"]=ng
        h["NAXIS2"]=nf
        h=row_data_header2(h,mjd)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'biases/boss_pixbias-'+str(mjd)+'-'+ty+'2.fits'
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)        


        h1=pyf.PrimaryHDU(f_arrayf_1)
        h=h1.header
        #h["NAXIS"]=2 
        #h["NAXIS1"]=ng
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'1',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'flats/pixflatave-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
    
        h1=pyf.PrimaryHDU(f_arrayf_2)
        h=h1.header
        #h["NAXIS1"]=ng
        #h["NAXIS"]=2 
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'2',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'flats/pixflatave-'+str(mjd)+'-'+ty+'2.fits'
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)    
        
        h1=pyf.PrimaryHDU(p_arrayf_1)
        h=h1.header
        #h["NAXIS"]=2 
        #h["NAXIS1"]=ng
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'1',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'flats/badpixels-'+str(mjd)+'-'+ty+'1.fits' 
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
    
        h1=pyf.PrimaryHDU(p_arrayf_2)
        h=h1.header
        #h["NAXIS"]=2 
        #h["NAXIS1"]=ng
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'2',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'flats/badpixels-'+str(mjd)+'-'+ty+'2.fits'
        wfits_ext(out_fit,hlist)
        sycall('gzip -f '+out_fit)
        
        h1=pyf.PrimaryHDU(b1_arrayf_1)
        h=h1.header
        #h["NAXIS"]=2 
        #h["NAXIS1"]=ng
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'1',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'biases/pixbiasave-00001'+'-'+ty+'1.fits'
        wfits_ext(out_fit,hlist)
        #sycall('gzip -f '+out_fit)
    
        h1=pyf.PrimaryHDU(b1_arrayf_2)
        h=h1.header
        #h["NAXIS"]=2 
        #h["NAXIS1"]=ng
        #h["NAXIS2"]=nf
        #h=row_data_header2(h,plate,mjd,exp,ty+'2',flb=flb)
        hlist=pyf.HDUList([h1])
        hlist.update_extend()
        out_fit=dir1+'biases/pixbiasave-00001'+'-'+ty+'2.fits'
        wfits_ext(out_fit,hlist)
        #sycall('gzip -f '+out_fit)
        #sys.exit()
        
def row_data_header2(h,mjd):                                                
    h['BSCALE']  = 1                                                  
    h['BZERO']   = 32768                                                  
    h['EXTEND']  =  True                                                 
    h['TELESCOP']= 'SDSS 2-5m'                                                           
    h['FILENAME']= 'sdR-b1-00104337.fit'                                                 
    h['CAMERAS'] = 'b1      '                                                            
    h['EXPOSURE']=  104337                                                  
    h['DAQVER']  = '1.2.7   '                                                            
    h['CAMDAQ']  = '1.2.0:28'                                                            
    h['ERRCNT']  = 'NONE    '                                                            
    h['SYNCERR'] = 'NONE    '                                                            
    h['SLINES']  = 'NONE    '                                                            
    h['PIXERR']  = 'NONE    '                                                            
    h['PLINES']  = 'NONE    '                                                            
    h['FLAVOR']  = ('bias    ','exposure type, SDSS spectro style')              
    h['BOSSVER'] = ('branch_jme-rewrite+svn105840M','ICC version')                         
    h['MJD']     = (np.int(mjd),'APO MJD day at start of exposure')  
    h['TAI-BEG'] = ((np.float(mjd)+0.25)*24.0*3600.0,'MJD(TAI) seconds at start of exposure')        
    h['DATE-OBS']= ('2012-03-20T06:00:00','TAI date at start of exposure')               
    h['FF']      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                           
    h['NE']      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                         
    h['HGCD']    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')       
    h['FFS']     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open')         
    h['OBJSYS']  = ('Mount   ','The TCC objSys')             
    h['RA']      = ('NaN     ','Telescope is not tracking the sky')   
    h['DEC']     = ('NaN     ','Telescope is not tracking the sky')            
    h['RADEG']   = ('NaN     ','Telescope is not tracking the sky')            
    h['DECDEG']  = ('NaN     ','Telescope is not tracking the sky')            
    h['ROTTYPE'] = ('Mount   ','Rotator request type')        
    h['ROTPOS']  = (0.0 ,'Rotator request position (deg)')                 
    h['BOREOFFX']= (0.0 ,'TCC Boresight offset, deg')               
    h['BOREOFFY']= (0.0 ,'TCC Boresight offset, deg')                     
    h['ARCOFFX'] = (0.0 ,'TCC ObjArcOff, deg')                    
    h['ARCOFFY'] = (0.0 ,'TCC ObjArcOff, deg')                            
    h['OBJOFFX'] = (0.0 ,'TCC ObjOff, deg')                           
    h['OBJOFFY'] = (0.0 ,'TCC ObjOff, deg')                            
    h['CALOFFX'] = (0.0 ,'TCC CalibOff, deg')                             
    h['CALOFFY'] = (0.0 ,'TCC CalibOff, deg')                             
    h['CALOFFR'] = (0.0 ,'TCC CalibOff, deg')                             
    h['GUIDOFFX']= (0.0 ,'TCC GuideOff, deg')                             
    h['GUIDOFFY']= (0.0 ,'TCC GuideOff, deg')                             
    h['GUIDOFFR']= (0.0 ,'TCC GuideOff, deg')             
    h['AZ']      = (121.0 ,'Azimuth axis pos. (approx, deg)')               
    h['ALT']     = (30.0 ,'Altitude axis pos. (approx, deg)')               
    h['IPA']     = (0.0 ,'Rotator axis pos. (approx, deg)')              
    h['FOCUS']   = (0.0 ,'User-specified focus offset (um)')             
    h['M2PISTON']= (1256.77 ,'TCC SecOrient')                  
    h['M2XTILT'] = (-3.21 ,'TCC SecOrient')                               
    h['M2YTILT'] = (-10.45 ,'TCC SecOrient')                                
    h['M2XTRAN'] = (9.24 ,'TCC SecOrient')                              
    h['M2YTRAN'] = (234.17 ,'TCC SecOrient')                                 
    h['M1PISTON']= (0.0 ,'TCC PrimOrient')                           
    h['M1XTILT'] = (-3.04 ,'TCC PrimOrient')                                 
    h['M1YTILT'] = (5.26 ,'TCC PrimOrient')                             
    h['M1XTRAN'] = (277.7 ,'TCC PrimOrient')                          
    h['M1YTRAN'] = (178.07 ,'TCC PrimOrient')                                 
    h['SCALE']   = (1.0 ,'User-specified scale factor')                   
    h['NAME']    = ('3521-55170-01','The name of the currently loaded plate')   
    h['PLATEID'] = (3521 ,'The currently loaded plate')
    h['CARTID']  = (11 ,'The currently loaded cartridge')                 
    h['MAPID']   = (1 ,'The mapping version of the loaded plate')      
    h['POINTING']= ('A       ','The currently specified pointing')               
    h['GUIDER1'] = ('proc-gimg-0135.fits','The first guider image')                        
    h['GUIDERN'] = ('proc-gimg-0135.fits','The last guider image')                     
    h['EXPTIME'] = (0.08755397796630859,'exposure time')              
    h['DATASUM'] = '0000000000000000'                                                    
    h['COMMENT'] = 'failed to make SPA card from None'                                     
    h['EXPID00'] = ('b1-00104337','ID string for exposure 00')                      
    h['EXPID01'] = ('b1-00104338','ID string for exposure 01')                      
    h['EXPID02'] = ('b1-00104339','ID string for exposure 02')                      
    h['EXPID03'] = ('b1-00104340','ID string for exposure 03')                      
    h['EXPID04'] = ('b1-00104341','ID string for exposure 04')                      
    h['EXPID05'] = ('b1-00104342','ID string for exposure 05')                      
    h['EXPID06'] = ('b1-00104343','ID string for exposure 06')                      
    h['EXPID07'] = ('b1-00104344','ID string for exposure 07')                       
    h['EXPID08'] = ('b1-00104345','ID string for exposure 08')                       
    h['EXPID09'] = ('b1-00104346','ID string for exposure 09')                       
    h['EXPID10'] = ('b1-00104347','ID string for exposure 10')                       
    h['EXPID11'] = ('b1-00104348','ID string for exposure 11')                       
    h['EXPID12'] = ('b1-00104349','ID string for exposure 12')                       
    h['EXPID13'] = ('b1-00104350','ID string for exposure 13')                       
    h['EXPID14'] = ('b1-00104351','ID string for exposure 14')                       
    h['EXPID15'] = ('b1-00104352','ID string for exposure 15')                       
    h['EXPID16'] = ('b1-00104353','ID string for exposure 16')                       
    h['EXPID17'] = ('b1-00104354','ID string for exposure 17')                       
    h['EXPID18'] = ('b1-00104355','ID string for exposure 18')                       
    h['EXPID19'] = ('b1-00104356','ID string for exposure 19')                       
    h['EXPID20'] = ('b1-00104357','ID string for exposure 20')                       
    h['EXPID21'] = ('b1-00104358','ID string for exposure 21')                       
    h['EXPID22'] = ('b1-00104359','ID string for exposure 22') 
    return h

def row_data_header2_bhm(h,mjd):                                                
    h['BSCALE']  = 1                                                  
    h['BZERO']   = 32768                                                  
    h['EXTEND']  =  True                                                 
    h['TELESCOP']= 'SDSS 2-5m'                                                           
    h['FILENAME']= 'sdR-b1-00104337.fit'                                                 
    h['CAMERAS'] = 'b1      '                                                            
    h['EXPOSURE']=  104337                                                  
    h['DAQVER']  = '1.2.7   '                                                            
    h['CAMDAQ']  = '1.2.0:28'                                                            
    h['ERRCNT']  = 'NONE    '                                                            
    h['SYNCERR'] = 'NONE    '                                                            
    h['SLINES']  = 'NONE    '                                                            
    h['PIXERR']  = 'NONE    '                                                            
    h['PLINES']  = 'NONE    '                                                            
    h['FLAVOR']  = ('bias    ','exposure type, SDSS spectro style')              
    h['BOSSVER'] = ('branch_jme-rewrite+svn105840M','ICC version')                         
    h['MJD']     = (np.int(mjd),'APO MJD day at start of exposure')  
    h['TAI-BEG'] = ((np.float(mjd)+0.25)*24.0*3600.0,'MJD(TAI) seconds at start of exposure')        
    h['DATE-OBS']= ('2012-03-20T06:00:00','TAI date at start of exposure')               
    h['FF']      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                           
    h['NE']      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                         
    h['HGCD']    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')       
    h['FFS']     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open')         
    h['OBJSYS']  = ('Mount   ','The TCC objSys')             
    h['RA']      = ('NaN     ','Telescope is not tracking the sky')   
    h['DEC']     = ('NaN     ','Telescope is not tracking the sky')            
    h['RADEG']   = ('NaN     ','Telescope is not tracking the sky')            
    h['DECDEG']  = ('NaN     ','Telescope is not tracking the sky')            
    h['ROTTYPE'] = ('Mount   ','Rotator request type')        
    h['ROTPOS']  = (0.0 ,'Rotator request position (deg)')                 
    h['BOREOFFX']= (0.0 ,'TCC Boresight offset, deg')               
    h['BOREOFFY']= (0.0 ,'TCC Boresight offset, deg')                     
    h['ARCOFFX'] = (0.0 ,'TCC ObjArcOff, deg')                    
    h['ARCOFFY'] = (0.0 ,'TCC ObjArcOff, deg')                            
    h['OBJOFFX'] = (0.0 ,'TCC ObjOff, deg')                           
    h['OBJOFFY'] = (0.0 ,'TCC ObjOff, deg')                            
    h['CALOFFX'] = (0.0 ,'TCC CalibOff, deg')                             
    h['CALOFFY'] = (0.0 ,'TCC CalibOff, deg')                             
    h['CALOFFR'] = (0.0 ,'TCC CalibOff, deg')                             
    h['GUIDOFFX']= (0.0 ,'TCC GuideOff, deg')                             
    h['GUIDOFFY']= (0.0 ,'TCC GuideOff, deg')                             
    h['GUIDOFFR']= (0.0 ,'TCC GuideOff, deg')             
    h['AZ']      = (121.0 ,'Azimuth axis pos. (approx, deg)')               
    h['ALT']     = (30.0 ,'Altitude axis pos. (approx, deg)')               
    h['IPA']     = (0.0 ,'Rotator axis pos. (approx, deg)')              
    h['FOCUS']   = (0.0 ,'User-specified focus offset (um)')             
    h['M2PISTON']= (1256.77 ,'TCC SecOrient')                  
    h['M2XTILT'] = (-3.21 ,'TCC SecOrient')                               
    h['M2YTILT'] = (-10.45 ,'TCC SecOrient')                                
    h['M2XTRAN'] = (9.24 ,'TCC SecOrient')                              
    h['M2YTRAN'] = (234.17 ,'TCC SecOrient')                                 
    h['M1PISTON']= (0.0 ,'TCC PrimOrient')                           
    h['M1XTILT'] = (-3.04 ,'TCC PrimOrient')                                 
    h['M1YTILT'] = (5.26 ,'TCC PrimOrient')                             
    h['M1XTRAN'] = (277.7 ,'TCC PrimOrient')                          
    h['M1YTRAN'] = (178.07 ,'TCC PrimOrient')                                 
    h['SCALE']   = (1.0 ,'User-specified scale factor')                   
    h['NAME']    = ('3521-55170-01','The name of the currently observation')   
    h['CONFIID'] = (3521 ,'The currently loaded plate')
    h['CARTID']  = (11 ,'The currently loaded cartridge')                 
    h['MAPID']   = (1 ,'The mapping version of the loaded plate')      
    h['POINTING']= ('A       ','The currently specified pointing')               
    h['GUIDER1'] = ('proc-gimg-0135.fits','The first guider image')                        
    h['GUIDERN'] = ('proc-gimg-0135.fits','The last guider image')                     
    h['EXPTIME'] = (0.08755397796630859,'exposure time')              
    h['DATASUM'] = '0000000000000000'                                                    
    h['COMMENT'] = 'failed to make SPA card from None'                                     
    h['EXPID00'] = ('b1-00104337','ID string for exposure 00')                      
    h['EXPID01'] = ('b1-00104338','ID string for exposure 01')                      
    h['EXPID02'] = ('b1-00104339','ID string for exposure 02')                      
    h['EXPID03'] = ('b1-00104340','ID string for exposure 03')                      
    h['EXPID04'] = ('b1-00104341','ID string for exposure 04')                      
    h['EXPID05'] = ('b1-00104342','ID string for exposure 05')                      
    h['EXPID06'] = ('b1-00104343','ID string for exposure 06')                      
    h['EXPID07'] = ('b1-00104344','ID string for exposure 07')                       
    h['EXPID08'] = ('b1-00104345','ID string for exposure 08')                       
    h['EXPID09'] = ('b1-00104346','ID string for exposure 09')                       
    h['EXPID10'] = ('b1-00104347','ID string for exposure 10')                       
    h['EXPID11'] = ('b1-00104348','ID string for exposure 11')                       
    h['EXPID12'] = ('b1-00104349','ID string for exposure 12')                       
    h['EXPID13'] = ('b1-00104350','ID string for exposure 13')                       
    h['EXPID14'] = ('b1-00104351','ID string for exposure 14')                       
    h['EXPID15'] = ('b1-00104352','ID string for exposure 15')                       
    h['EXPID16'] = ('b1-00104353','ID string for exposure 16')                       
    h['EXPID17'] = ('b1-00104354','ID string for exposure 17')                       
    h['EXPID18'] = ('b1-00104355','ID string for exposure 18')                       
    h['EXPID19'] = ('b1-00104356','ID string for exposure 19')                       
    h['EXPID20'] = ('b1-00104357','ID string for exposure 20')                       
    h['EXPID21'] = ('b1-00104358','ID string for exposure 21')                       
    h['EXPID22'] = ('b1-00104359','ID string for exposure 22') 
    return h

def row_data_header_bhm(h,plate,mjd,exp,typ,flb='s',ra=0.0,dec=0.0,azim=180.0,alt=90.0,expt=900.0,expof=0.0):  
    if flb == 's':
        flab='science '
    elif flb == 'a':
        flab='arc     '
        expt=4.0
    elif flb == 'f':
        flab='flat    '  
        expt=25.0                                        
    h["TELESCOP"]= 'SDSS 2-5m'                                                           
    h["FILENAME"]= 'sdR-'+typ+'-'+id_str(exp,n_z=8)+'.fit'                                                 
    h["CAMERAS"] = typ+'      '                                                            
    h["EXPOSURE"]=  exp                                                  
    h["V_BOSS"]  = ('v4_0    ' ,'Active version of the BOSS ICC')              
    h["CAMDAQ"]  = '1.5.0:37'                                                            
    h["SUBFRAME"]= ('' ,'the subframe readout command')                                    
    h["ERRCNT"]  = 'NONE    '                                                            
    h["SYNCERR"] = 'NONE    '                                                            
    h["SLINES"]  = 'NONE    '                                                            
    h["PIXERR"]  = 'NONE    '                                                            
    h["PLINES"]  = 'NONE    '                                                            
    h["PFERR"]   = 'NONE    '                                                            
    h["DIDFLUSH"]= (True ,'CCD was flushed before integration')          
    h["FLAVOR"]  = (flab,'exposure type, SDSS spectro style')            
    h["MJD"]     = (np.int(mjd) ,'APO fMJD day at start of exposure')          
    h["TAI-BEG"] = ((np.float(mjd)+0.25)*24.0*3600.0+expof,'MJD(TAI) seconds at start of integration')  
    h["DATE-OBS"]= ('2012-03-20T06:00:00','TAI date at start of integration')           
    h["V_GUIDER"]= ('v3_4    ','version of the current guiderActor')            
    h["V_SOP"]   = ('v3_8_1  ','version of the current sopActor')           
    h["NAME"]    = (plate+'-'+mjd+'-01','The name of the currently loaded plate')     
    h["CONFIID"] = (plate,'The currently FPS configuration')              
    h["CARTID"]  = (16,'The currently loaded cartridge')                 
    h["MAPID"]   = (1,'The mapping version of the loaded plate') 
    h["POINTING"]= ('A       ','The currently specified pointing')               
    h["CONFTYP"]= ('BOSS    ','Type of plate (e.g. BOSS, APOGEE, BA')
    h["SRVYMODE"]= ('None    ','Survey leading this observation and its mode')
    h["OBJSYS"]  = ('ICRS    ','The TCC objSys')             
    h["RA"]      = (ra,'RA of telescope boresight (deg)')             
    h["DEC"]     = (dec,'Dec of telescope boresight (deg)')              
    h["RADEG"]   = (ra+0.704,'RA of telescope pointing(deg)')                 
    h["DECDEG"]  = (dec+0.083,'Dec of telescope pointing (deg)')                
    h["SPA"]     = (-158.0698343797722,'TCC SpiderInstAng')                              
    h["ROTTYPE"] = ('Obj     ','Rotator request type')                           
    h["ROTPOS"]  = (0.0,'Rotator request position (deg)')            
    h["BOREOFFX"]= (0.0,'TCC Boresight offset, deg')               
    h["BOREOFFY"]= (0.0,'TCC Boresight offset, deg')                    
    h["ARCOFFX"] = (-8.8999999999999E-05,'TCC ObjArcOff, deg')                    
    h["ARCOFFY"] = (-0.000807,'TCC ObjArcOff, deg')                           
    h["CALOFFX"] = (0.0,'TCC CalibOff, deg')                           
    h["CALOFFY"] = (0.0,'TCC CalibOff, deg')                            
    h["CALOFFR"] = (0.0,'TCC CalibOff, deg')                            
    h["GUIDOFFX"]= (0.0,'TCC GuideOff, deg')                            
    h["GUIDOFFY"]= (0.0,'TCC GuideOff, deg')                            
    h["GUIDOFFR"]= (0.052684,'TCC GuideOff, deg')                            
    h["AZ"]      = (azim,'Azimuth axis pos. (approx, deg)')           
    h["ALT"]     = (alt,'Altitude axis pos. (approx, deg)')              
    h["IPA"]     = (21.60392,'Rotator axis pos. (approx, deg)')             
    h["FOCUS"]   = (10.7512,'User-specified focus offset (um)')              
    h["M2PISTON"]= (357.36,'TCC SecOrient')             
    h["M2XTILT"] = (7.19,'TCC SecOrient')                                
    h["M2YTILT"] = (-18.2,'TCC SecOrient')                                
    h["M2XTRAN"] = (24.89,'TCC SecOrient')                                
    h["M2YTRAN"] = (-110.34,'TCC SecOrient')                                
    h["M2ZROT"]  = (-19.77,'TCC SecOrient')                                
    h["M1PISTON"]= (-949.28,'TCC PrimOrient')                                
    h["M1XTILT"] = (-24.31,'TCC PrimOrient')                               
    h["M1YTILT"] = (6.14,'TCC PrimOrient')                               
    h["M1XTRAN"] = (356.01,'TCC PrimOrient')                               
    h["M1YTRAN"] = (1322.6,'TCC PrimOrient')                               
    h["M1ZROT"]  = (0.03,'TCC PrimOrient')                  
    h["SCALE"]   = (1.000096,'User-specified scale factor')              
    h["V_APO"]   = ('trunk+svn158476M','version of the current apoActor')              
    h["PRESSURE"]=21.413                                                  
    h["WINDD"]   =286.0                                                  
    h["WINDS"]   =18.6                                                  
    h["GUSTD"]   =295.6                                                  
    h["GUSTS"]   =25.1                                                  
    h["AIRTEMP2"]=8.1                                                  
    h["DEWPOINT"]=-4.2                                                  
    h["TRUSTEMP"]=7.79                                                  
    h["HUMIDITY"]=39.9                                                  
    h["DUSTA"]   =16084.0                                                  
    h["DUSTB"]   =1020.0                                                  
    h["WINDD25M"]=318.3                                                  
    h["WINDS25M"]=1.4    
    if 'flat    ' in flab:                                              
        h["FF"]      = ('1 1 1 1 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open')    
    elif 'science ' in flab:
        h["FF"]      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('0 0 0 0 0 0 0 0','Flatfield Screen 1:closed 0:open')                     
    elif 'arc     ' in flab:
        h["FF"]      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('1 1 1 1 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('1 1 1 1 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open') 
    h["MGDPOS"]  = ('C       ','MaNGA dither position (C,N,S,E)')             
    h["MGDRA"]   = (0.0,'MaNGA decenter in RA, redundant with MGDPOS') 
    h["MGDDEC"]  = (0.0,'MaNGA decenter in Dec, redundant with MGDPOS')  
    h["GUIDER1"] = ('proc-gimg-0500.fits.gz','The first guider image')              
    h["SLITID1"] = (16,'Normalized slithead ID. sp1&2 should match.') 
    h["SLITID2"] = (16,'Normalized slithead ID. sp1&2 should match.')  
    h["GUIDERN"] = ('proc-gimg-0529.fits.gz','The last guider image')  
    h["COLLA"]   = (1173,'The position of the A collimator motor')         
    h["COLLB"]   = (164,'The position of the B collimator motor')       
    h["COLLC"]   = (577,'The position of the C collimator motor')       
    h["HARTMANN"]= ('Out     ','Hartmanns: Left,Right,Out')   
    if '2' in typ:    
        h["MC2HUMHT"]= (32.3,'sp2 mech Hartmann humidity, %')                  
        h["MC2HUMCO"]= (25.6,'sp2 mech Central optics humidity, %')           
        h["MC2TEMDN"]= (7.2,'sp2 mech Median temp, C')                
        h["MC2THT"]  = (7.5,'sp2 mech Hartmann Top Temp, C')          
        h["MC2TRCB"] = (7.4,'sp2 mech Red Cam Bottom Temp, C')                
        h["MC2TRCT"] = (6.9,'sp2 mech Red Cam Top Temp, C')              
        h["MC2TBCB"] = (7.1,'sp2 mech Blue Cam Bottom Temp, C')               
        h["MC2TBCT"] = (7.2,'sp2 mech Blue Cam Top Temp, C')   
    else:
        h["MC2HUMHT"]= (32.3,'sp1 mech Hartmann humidity, %')                  
        h["MC2HUMCO"]= (25.6,'sp1 mech Central optics humidity, %')           
        h["MC2TEMDN"]= (7.2,'sp1 mech Median temp, C')                
        h["MC2THT"]  = (7.5,'sp1 mech Hartmann Top Temp, C')          
        h["MC2TRCB"] = (7.4,'sp1 mech Red Cam Bottom Temp, C')                
        h["MC2TRCT"] = (6.9,'sp1 mech Red Cam Top Temp, C')              
        h["MC2TBCB"] = (7.1,'sp1 mech Blue Cam Bottom Temp, C')               
        h["MC2TBCT"] = (7.2,'sp1 mech Blue Cam Top Temp, C')                 
    h["REQTIME"] = (expt,'requested exposure time')             
    h["EXPTIME"] = (expt+0.14,'measured exposure time, s')                      
    h["SHOPETIM"]= (0.6899999999999999,'open shutter transit time, s')                   
    h["SHCLOTIM"]= (0.63,'close shutter transit time, s')                  
    h["DARKTIME"]= (expt+9.4519929885864,'time between flush end and readout start')       
    h["LN2TEMP"] = 81.64100000000001                                                  
    h["CCDTEMP"] = -133.984                                                  
    h["IONPUMP"] = -6.17                                                  
    h["BSCALE"]  = 1                                                  
    h["BZERO"]   = 32768                                                  
    h["CHECKSUM"]= ('DrANEo1KDo8KDo8K','HDU checksum updated 2016-05-10T06:58:02')   
    h["DATASUM"] = ('516485492','data unit checksum updated 2016-05-10T06:58:02')                                                                          
    return h

    
def row_data_header(h,plate,mjd,exp,typ,flb='s',ra=0.0,dec=0.0,azim=180.0,alt=90.0,expt=900.0,expof=0.0):  
    if flb == 's':
        flab='science '
    elif flb == 'a':
        flab='arc     '
        expt=4.0
    elif flb == 'f':
        flab='flat    '  
        expt=25.0                                        
    h["TELESCOP"]= 'SDSS 2-5m'                                                           
    h["FILENAME"]= 'sdR-'+typ+'-'+id_str(exp,n_z=8)+'.fit'                                                 
    h["CAMERAS"] = typ+'      '                                                            
    h["EXPOSURE"]=  exp                                                  
    h["V_BOSS"]  = ('v4_0    ' ,'Active version of the BOSS ICC')              
    h["CAMDAQ"]  = '1.5.0:37'                                                            
    h["SUBFRAME"]= ('' ,'the subframe readout command')                                    
    h["ERRCNT"]  = 'NONE    '                                                            
    h["SYNCERR"] = 'NONE    '                                                            
    h["SLINES"]  = 'NONE    '                                                            
    h["PIXERR"]  = 'NONE    '                                                            
    h["PLINES"]  = 'NONE    '                                                            
    h["PFERR"]   = 'NONE    '                                                            
    h["DIDFLUSH"]= (True ,'CCD was flushed before integration')          
    h["FLAVOR"]  = (flab,'exposure type, SDSS spectro style')            
    h["MJD"]     = (np.int(mjd) ,'APO fMJD day at start of exposure')          
    h["TAI-BEG"] = ((np.float(mjd)+0.25)*24.0*3600.0+expof,'MJD(TAI) seconds at start of integration')  
    h["DATE-OBS"]= ('2012-03-20T06:00:00','TAI date at start of integration')           
    h["V_GUIDER"]= ('v3_4    ','version of the current guiderActor')            
    h["V_SOP"]   = ('v3_8_1  ','version of the current sopActor')           
    h["NAME"]    = (plate+'-'+mjd+'-01','The name of the currently loaded plate')     
    h["PLATEID"] = (np.int(plate),'The currently loaded plate')              
    h["CARTID"]  = (16,'The currently loaded cartridge')                 
    h["MAPID"]   = (1,'The mapping version of the loaded plate') 
    h["POINTING"]= ('A       ','The currently specified pointing')               
    h["PLATETYP"]= ('BOSS    ','Type of plate (e.g. BOSS, MANGA, APOGEE, APOGEE')
    h["SRVYMODE"]= ('None    ','Survey leading this observation and its mode')
    h["OBJSYS"]  = ('ICRS    ','The TCC objSys')             
    h["RA"]      = (ra,'RA of telescope boresight (deg)')             
    h["DEC"]     = (dec,'Dec of telescope boresight (deg)')              
    h["RADEG"]   = (ra+0.704,'RA of telescope pointing(deg)')                 
    h["DECDEG"]  = (dec+0.083,'Dec of telescope pointing (deg)')                
    h["SPA"]     = (-158.0698343797722,'TCC SpiderInstAng')                              
    h["ROTTYPE"] = ('Obj     ','Rotator request type')                           
    h["ROTPOS"]  = (0.0,'Rotator request position (deg)')            
    h["BOREOFFX"]= (0.0,'TCC Boresight offset, deg')               
    h["BOREOFFY"]= (0.0,'TCC Boresight offset, deg')                    
    h["ARCOFFX"] = (-8.8999999999999E-05,'TCC ObjArcOff, deg')                    
    h["ARCOFFY"] = (-0.000807,'TCC ObjArcOff, deg')                           
    h["CALOFFX"] = (0.0,'TCC CalibOff, deg')                           
    h["CALOFFY"] = (0.0,'TCC CalibOff, deg')                            
    h["CALOFFR"] = (0.0,'TCC CalibOff, deg')                            
    h["GUIDOFFX"]= (0.0,'TCC GuideOff, deg')                            
    h["GUIDOFFY"]= (0.0,'TCC GuideOff, deg')                            
    h["GUIDOFFR"]= (0.052684,'TCC GuideOff, deg')                            
    h["AZ"]      = (azim,'Azimuth axis pos. (approx, deg)')           
    h["ALT"]     = (alt,'Altitude axis pos. (approx, deg)')              
    h["IPA"]     = (21.60392,'Rotator axis pos. (approx, deg)')             
    h["FOCUS"]   = (10.7512,'User-specified focus offset (um)')              
    h["M2PISTON"]= (357.36,'TCC SecOrient')             
    h["M2XTILT"] = (7.19,'TCC SecOrient')                                
    h["M2YTILT"] = (-18.2,'TCC SecOrient')                                
    h["M2XTRAN"] = (24.89,'TCC SecOrient')                                
    h["M2YTRAN"] = (-110.34,'TCC SecOrient')                                
    h["M2ZROT"]  = (-19.77,'TCC SecOrient')                                
    h["M1PISTON"]= (-949.28,'TCC PrimOrient')                                
    h["M1XTILT"] = (-24.31,'TCC PrimOrient')                               
    h["M1YTILT"] = (6.14,'TCC PrimOrient')                               
    h["M1XTRAN"] = (356.01,'TCC PrimOrient')                               
    h["M1YTRAN"] = (1322.6,'TCC PrimOrient')                               
    h["M1ZROT"]  = (0.03,'TCC PrimOrient')                  
    h["SCALE"]   = (1.000096,'User-specified scale factor')              
    h["V_APO"]   = ('trunk+svn158476M','version of the current apoActor')              
    h["PRESSURE"]=21.413                                                  
    h["WINDD"]   =286.0                                                  
    h["WINDS"]   =18.6                                                  
    h["GUSTD"]   =295.6                                                  
    h["GUSTS"]   =25.1                                                  
    h["AIRTEMP2"]=8.1                                                  
    h["DEWPOINT"]=-4.2                                                  
    h["TRUSTEMP"]=7.79                                                  
    h["HUMIDITY"]=39.9                                                  
    h["DUSTA"]   =16084.0                                                  
    h["DUSTB"]   =1020.0                                                  
    h["WINDD25M"]=318.3                                                  
    h["WINDS25M"]=1.4    
    if 'flat    ' in flab:                                              
        h["FF"]      = ('1 1 1 1 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open')    
    elif 'science ' in flab:
        h["FF"]      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('0 0 0 0 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('0 0 0 0 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('0 0 0 0 0 0 0 0','Flatfield Screen 1:closed 0:open')                     
    elif 'arc     ' in flab:
        h["FF"]      = ('0 0 0 0 ','FF lamps 1:on 0:0ff')                       
        h["NE"]      = ('1 1 1 1 ','NE lamps 1:on 0:0ff')                          
        h["HGCD"]    = ('1 1 1 1 ','HGCD lamps 1:on 0:0ff')                          
        h["FFS"]     = ('1 1 1 1 1 1 1 1','Flatfield Screen 1:closed 0:open') 
    h["MGDPOS"]  = ('C       ','MaNGA dither position (C,N,S,E)')             
    h["MGDRA"]   = (0.0,'MaNGA decenter in RA, redundant with MGDPOS') 
    h["MGDDEC"]  = (0.0,'MaNGA decenter in Dec, redundant with MGDPOS')  
    h["GUIDER1"] = ('proc-gimg-0500.fits.gz','The first guider image')              
    h["SLITID1"] = (16,'Normalized slithead ID. sp1&2 should match.') 
    h["SLITID2"] = (16,'Normalized slithead ID. sp1&2 should match.')  
    h["GUIDERN"] = ('proc-gimg-0529.fits.gz','The last guider image')  
    h["COLLA"]   = (1173,'The position of the A collimator motor')         
    h["COLLB"]   = (164,'The position of the B collimator motor')       
    h["COLLC"]   = (577,'The position of the C collimator motor')       
    h["HARTMANN"]= ('Out     ','Hartmanns: Left,Right,Out')   
    if '2' in typ:    
        h["MC2HUMHT"]= (32.3,'sp2 mech Hartmann humidity, %')                  
        h["MC2HUMCO"]= (25.6,'sp2 mech Central optics humidity, %')           
        h["MC2TEMDN"]= (7.2,'sp2 mech Median temp, C')                
        h["MC2THT"]  = (7.5,'sp2 mech Hartmann Top Temp, C')          
        h["MC2TRCB"] = (7.4,'sp2 mech Red Cam Bottom Temp, C')                
        h["MC2TRCT"] = (6.9,'sp2 mech Red Cam Top Temp, C')              
        h["MC2TBCB"] = (7.1,'sp2 mech Blue Cam Bottom Temp, C')               
        h["MC2TBCT"] = (7.2,'sp2 mech Blue Cam Top Temp, C')   
    else:
        h["MC2HUMHT"]= (32.3,'sp1 mech Hartmann humidity, %')                  
        h["MC2HUMCO"]= (25.6,'sp1 mech Central optics humidity, %')           
        h["MC2TEMDN"]= (7.2,'sp1 mech Median temp, C')                
        h["MC2THT"]  = (7.5,'sp1 mech Hartmann Top Temp, C')          
        h["MC2TRCB"] = (7.4,'sp1 mech Red Cam Bottom Temp, C')                
        h["MC2TRCT"] = (6.9,'sp1 mech Red Cam Top Temp, C')              
        h["MC2TBCB"] = (7.1,'sp1 mech Blue Cam Bottom Temp, C')               
        h["MC2TBCT"] = (7.2,'sp1 mech Blue Cam Top Temp, C')                 
    h["REQTIME"] = (expt,'requested exposure time')             
    h["EXPTIME"] = (expt+0.14,'measured exposure time, s')                      
    h["SHOPETIM"]= (0.6899999999999999,'open shutter transit time, s')                   
    h["SHCLOTIM"]= (0.63,'close shutter transit time, s')                  
    h["DARKTIME"]= (expt+9.4519929885864,'time between flush end and readout start')       
    h["LN2TEMP"] = 81.64100000000001                                                  
    h["CCDTEMP"] = -133.984                                                  
    h["IONPUMP"] = -6.17                                                  
    h["BSCALE"]  = 1                                                  
    h["BZERO"]   = 32768                                                  
    h["CHECKSUM"]= ('DrANEo1KDo8KDo8K','HDU checksum updated 2016-05-10T06:58:02')   
    h["DATASUM"] = ('516485492','data unit checksum updated 2016-05-10T06:58:02')                                                                          
    return h

def star_pool(str='none'):
    start_pool=np.array(['F0II (25291)','F0III (89025)','F0IV (81937)','F0Ib (36673)','F0V (90277)','F2III (89254)','F2V (33256)','F3/F5V (30743)','F5Ib... (17463)','F6II (61295)','F6III (61064)','F6Iab: (187929)','F6V (16673)','F8Ibvar (45412)','F8V (30562)','F8V (G_243-63)','F9IV (136064)'])
    if not 'none' in str:
        for i in range(0, len(start_pool)):
            if str in start_pool[i]: 
                star_t=start_pool[i]
    else:
        star_t=start_pool[np.int(ran.rand(1)[0]*len(start_pool))]       
    return star_t
    

def mock_arc(dirtemp="./",):
    spec_b,spec_r=fib_arc_raw(dirtemp=dirtemp)
    return spec_b,spec_r

def mock_flat(dirtemp="./",):
    spec_b,spec_r=fib_flat_raw(dirtemp=dirtemp)
    return spec_b,spec_r

def mock_sky(idh,xx,yy,fib_id,idn='spectras',dirtemp="./",sp_res=2000.0,basename='star-',template2="libs/sky_model.txt",dir1='',SN=15.0,Fluxm=20.0,plots=1,ifutype="BOSS",expt=900.0):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    cubef=basename+id+'.sky'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+idn+'/'
    dir1=dir1+idn+'/'
    sycall('mkdir -p '+dir1t)
    dir1t=dir1t+'spectras/'
    dir1=dir1+'spectras/'   
    sycall('mkdir -p '+dir1t)    
    spec_b,spec_r=fib_sky_raw(cubef,xx,yy,fib_id,dirtemp=dirtemp,template4=template2,sp_res=sp_res,SNi=SN,Flux_m=Fluxm,dir_o=dir1,ifutype=ifutype,pdf=plots,expt=expt)
    return spec_b,spec_r

def mock_star(idh,xx,yy,fib_id,dsep=1.0,lshf=5070.0,mjd='56000',idn='spectras',star_t='F0II (25291)',ra_0=180.0,expof=0.0,Mag_s=17.5,dirtemp="./",Av_g=0.00,sp_res=2000.0,basename='star-',template2="libs/sky_model.txt",template="libs/spEigenStar-55734.fits",template_wd="libs/da012500_800.dat",wd=False,dir1='',psf=0,SN=15.0,Fluxm=20.0,plots=1,ifutype="BOSS",expt=900.0,outs=1,expot=3600.0,toffs=0.0,t_exp=0,apot=0,alp=0.0,bet=0.0,dec_0=0):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    cubef=basename+id+'.spec'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+idn+'/'
    dir1=dir1+idn+'/'
    sycall('mkdir -p '+dir1t)
    dir1t=dir1t+'spectras/'
    dir1=dir1+'spectras/'   
    sycall('mkdir -p '+dir1t)
    ha=lst_c(tai=((np.float(mjd)+.25)*3600.0*24.0+expof+expt/2.0))-ra_0/180.0*12.0 
    if t_exp == 0:
        ha1=ha
    elif t_exp == 1:
        ha1=lst_c(tai=((np.float(mjd)+0.25)*3600.0*24.0+toffs+expot/2.0))-ra_0/180.0*12.0
    else:
        ha1=0
    vel=0.0#+ran.randn(1)*0.1
    #Mag_s=Mag_s+ran.randn(1)[0]*0.1   
    dir1n=dir1.replace(" ","\ ")
    if  outs == 0:
        spec_b,spec_r,mag_v=fib_star_raw(cubef,Av_g,Mag_s,xx,yy,fib_id,dsep=dsep,lshf=lshf,dirtemp=dirtemp,alp=alp,bet=bet,ha=ha,ha1=ha1,star_t=star_t,template4=template2,template3=template,template6=template_wd,wd=wd,sp_res=sp_res,vel_0=vel,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,ifutype=ifutype,pdf=plots,expt=expt,outs=outs,apot=apot,dec_0=dec_0)
    else:
        if ptt.exists(dir1+cubef+'.fits.gz') == False:
            spec_b,spec_r,mag_v=fib_star_raw(cubef,Av_g,Mag_s,xx,yy,fib_id,dsep=dsep,lshf=lshf,dirtemp=dirtemp,alp=alp,bet=bet,ha=ha,ha1=ha1,star_t=star_t,template4=template2,template3=template,template6=template_wd,wd=wd,sp_res=sp_res,vel_0=vel,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,ifutype=ifutype,pdf=plots,expt=expt,apot=apot,dec_0=dec_0)
            sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        else:
            sycall('rm '+dir1+cubef+'*.fits.gz')
            spec_b,spec_r,mag_v=fib_star_raw(cubef,Av_g,Mag_s,xx,yy,fib_id,dsep=dsep,lshf=lshf,dirtemp=dirtemp,alp=alp,bet=bet,ha=ha,ha1=ha1,star_t=star_t,template4=template2,template3=template,template6=template_wd,wd=wd,sp_res=sp_res,vel_0=vel,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,ifutype=ifutype,pdf=plots,expt=expt,apot=apot,dec_0=dec_0)
            sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
    return spec_b,spec_r,mag_v

    
def mock_gali(idh,idhp,idt,xx,yy,fib_id,dir1='',mjd='56000',lshf=5070.0,dirtemp="libs/",template4="libs/sky_model.txt",outs=1,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",SN=15.0,Fluxm=20.0,psf=0,ho=0.704,Lam=0.7274,Om=0.2726,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="MaNGA",basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1',ra_0=180.0,expof=0.0,pdf=1,expt=900.0,Av_gal=0.0,expot=3600.0,toffs=0.0,t_exp=0,dec_0=0):
    ha=lst_c(tai=((np.float(mjd)+0.25)*3600.0*24.0+expof+expt/2.0))-ra_0/180.0*12.0
    mjdt=((np.float(mjd)+0.25)*3600.0*24.0+expof+expt/2.0)/(3600.0*24.0)
    if t_exp == 0:
        ha1=ha
    elif t_exp == 1:
        ha1=lst_c(tai=((np.float(mjd)+0.25)*3600.0*24.0+toffs+expot/2.0))-ra_0/180.0*12.0
    else:
        ha1=0.0
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idhp)
    idt=str(idt)
    #outf='image_'+id
    #cubef='ilust-'+str(idhp)+'-'+idt+'.spec'
    cubef='sdR-'+idt+'.spec'
    #dirf=dir1t+'spec_lin/'
    #sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    dir1t=dir1t+'spectras/'
    dir1=dir1+'spectras/'   
    sycall('mkdir -p '+dir1t)
    stars = il.snapshot.loadSubhalo(basePath,135,idh,'stars')
    gas = il.snapshot.loadSubhalo(basePath,135,idh,'gas')
    bhs = il.snapshot.loadSubhalo(basePath,135,idh,'blackholes')
    dx_g = gas['Coordinates'][:,0]#
    dy_g = gas['Coordinates'][:,1]#
    dz_g = gas['Coordinates'][:,2]#
    vx_g = gas['Velocities'][:,0]#
    vy_g = gas['Velocities'][:,1]#
    vz_g = gas['Velocities'][:,2]#
    volm = gas['Volume'][:]
    dens = gas['Density'][:]#
    sfri = gas['StarFormationRate'][:]#
    meta_g=gas['GFM_Metallicity'][:]#     
    dx = stars['Coordinates'][:,0]#
    dy = stars['Coordinates'][:,1]#
    dz = stars['Coordinates'][:,2]#
    phot=stars['GFM_StellarPhotometrics'][:,5]
    mass=stars['Masses'][:]#
    mas0=stars['GFM_InitialMass'][:]#
    meta=stars['GFM_Metallicity'][:]#
    denS=stars['SubfindDensity']
    time=stars['GFM_StellarFormationTime'][:]#
    vx = stars['Velocities'][:,0]#
    vy = stars['Velocities'][:,1]#
    vz = stars['Velocities'][:,2]#
    mass_g=gas['Masses'][:]
    TE =   gas['InternalEnergy'][:]
    if bhs['count'] > 0:
        bhmass = bhs['BH_Mass'][:]#
        bhacre =bhs['BH_Mdot'][:]#
    else:
        bhmass=0.0
        bhacre=0.0
    #ages=stars['GFM_StellarAge'][:]
    sfri=sfri*100.0#Check this part
    temp_g=TE*((1e5)**2.0)*(0.938*1.7827e-24)*(4.0/(8.0-5.0*0.245))/(1.3807e-16)
    nt=np.where(time > 0)[0]
    #print len(nt)
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    phot=phot[nt]
    mass_g=mass_g/ho*1e10
    mass=mass[nt]/ho*1e10
    #print np.log10(np.sum(mass))
    mas0=mas0[nt]/ho*1e10
    bhmass=bhmass*1e10/ho
    bhacre=bhacre*1e10/0.978/1e9
    meta=meta[nt]#/0.0127
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=volm/ho**3.0
    dens=dens*ho**2.0
    volm=float_((volm/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    dens=dens*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    Av_g=meta_g*(3.0*1.67262e-24*np.pi*dens*volm)/(4.0*np.log(10.)*3.0*5494e-8)
   # print meta_g,(10.0**(-0.59)*0.0127)
    nt1=np.where(meta_g > (10.0**(-0.59)*0.0127))[0]
    nt2=np.where(meta_g <= (10.0**(-0.59)*0.0127))[0]
    Av_g[nt1]=Av_g[nt1]/(10.0**(2.21-1.0))
    if len(nt2) > 0 :
        Av_g[nt2]=Av_g[nt2]/(10.0**(2.21-1.0)/(meta_g[nt2]/0.0127)**(3.1-1.0))
    #print np.amax(meta),np.amin(meta)
    zf=1/time[nt]-1
    xo=modes(dx,nit=7)
    yo=modes(dy,nit=7)
    zo=modes(dz,nit=7)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    #print Rc/1e3
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    #print np.dot(np.dot(R2,R1),np.array([[x0],[y0],[z0]]))
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    #print red_0
    #print len(zf)
    age_s=cd.lookback_time(zf, **cosmo)/31557600./1e9
    #[age_F,ind]=ages_definition(age_s,n_ages=55)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:    
        spec_b,spec_r,mag_v=fib_conv_raw(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,xx,yy,fib_id,bhmass,bhacre,lshf=lshf,dirtemp=dirtemp,template4=template4,outs=outs,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,ifutype=ifutype,expt=expt,pdf=pdf,Av_gal=Av_gal,ha=ha,ha1=ha1,dec_0=dec_0,ra_0=ra_0,mjd=mjdt)
        #fib_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=fib_n,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        #sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        #band_cube(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
    else:
        if outs == 1:
            sycall('rm '+dir1+cubef+'.fits.gz')
        spec_b,spec_r,mag_v=fib_conv_raw(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,xx,yy,fib_id,bhmass,bhacre,lshf=lshf,dirtemp=dirtemp,template4=template4,outs=outs,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,ifutype=ifutype,expt=expt,pdf=pdf,Av_gal=Av_gal,ha=ha,ha1=ha1,dec_0=dec_0,ra_0=ra_0,mjd=mjdt)
        #band_cube(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
        #sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)   
    return spec_b,spec_r,mag_v 
    
def mock_sill(fov_p=3.0,rw=1,scp_s=60.4,ra_o=180.0,dec_o=0.0,ntar=800,dir1='',basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'): 
    #obs=[-106534,-106534,-106534]
    import random as rat
    obs=[ 52355.0, 43067.0, 43718.0]
    obs=[-106534.0,-106534.0,-106534.0]
    obs=[50000.0,50000.0,-50000.0]
    obs=[0.0,0.0,200000.0]
    plots=1
    f=h5py.File(il.snapshot.snapPath(basePath,135),'r')
    header = dict( f['Header'].attrs.items() )
    f.close()
    Om=header['Omega0']
    Lam=header['OmegaLambda']
    ho=0.74
    HaloID = il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloGrNr'])
    SubHaloMas= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloMassInRadType'])
    Len= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloLenType'])
    Rad=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloHalfmassRad'])#SubhaloStellarPhotometricsRad'])
    XYZ=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloCM'])
    MAGs=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloStellarPhotometrics'])
    Lent=Len[:,0]
    x=XYZ[:,0]
    y=XYZ[:,1]
    z=XYZ[:,2]
    mag=MAGs[:,5]
    #Rad=Rad_T[:,4]
    Fin_mass=np.log10(SubHaloMas[:,4]/ho*1e10+1)
    n_sub=len(HaloID)
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    
    print np.amax(y),np.amin(y)
    xt=x-(np.amax(x)+np.amin(x))/2.0
    yt=y-(np.amax(y)+np.amin(y))/2.0
    zt=z-(np.amax(z)+np.amin(z))/2.0
    print np.amax(xt),np.amin(xt)
    x0=obs[0]#-xt/ho
    y0=obs[1]#-yt/ho
    z0=obs[2]#-zt/ho
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    cam=Rc
    rad=np.sqrt(xt**2.+yt**2.+(cam-zt)**2.)
    reds=reds_cos(rad/1e3)
    print red_0,Rc,'H'
    radA=rad/(1+reds)
    phi=np.arcsin(xt/radA)/2.5
    the=np.arcsin(yt/(radA*np.cos(phi)))/2.5
    the=the*180/np.pi+ra_o#*3600
    phi=phi*180/np.pi+dec_o#*3600
    phi_ma=np.amax(phi)
    phi_mi=np.amin(phi)
    the_ma=np.amax(the)
    the_mi=np.amin(the)
    the_o=ra_o#1.2+180.0
    phi_o=dec_o#1.0
    phi_o=(phi_ma+phi_mi)/2.0
    the_o=(the_ma+the_mi)/2.0
    print the_o,phi_o,phi_ma-phi_o,phi_mi-phi_o
    xp=(the-the_o)*3600.0*scp_s/1e3
    yp=(phi-phi_o)*3600.0*scp_s/1e3
    
    mag_a=mag+5.0*np.log10(rad*1e3)-5.0
    #print mag_a
    
    dm=0.5
    Ms=14.5
    Mi=8.0
    nr=500
    nx=int((Ms-Mi)/dm)
    indf=[]
    if ptt.exists("Selection.dat") == False or rw == 1:
        f=open("Selection.dat","w")
        for i in range(0, nx):
            nt=np.where((Fin_mass > (Mi+dm*i)) & (Fin_mass <= (Mi+dm*(i+1))))[0]
            temp=nt#[rat.sample(range(0, len(nt)), nr)]
            for jj in range(0, len(nt)):#nr):
                if mag_a[temp[jj]] < 22.7 and Lent[temp[jj]] > 100 and mag_a[temp[jj]] > 11.6 :
                    #print mag_a[temp[jj]]
                    indf.extend([temp[jj]])
                    f.write(str(temp[jj])+'\n')
    else:
        f=open("Selection.dat","r")
        for line in f:
            indf.extend([int(line.replace('\n',''))])
    indf=np.array(indf)
    f.close()
    if ptt.exists('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt') == False or rw == 1:
        nt=np.where(np.sqrt((the[indf]-the_o)**2.0+(phi[indf]-phi_o)**2.0) < fov_p/2.0)[0]
        print len(nt),"H"
        
        nt0=range(0, len(nt))
        nt1=[]
        ntt1=[]
        nt2=[]
        ntt2=[]
        ntt=rat.sample(nt0, ntar)
        if len(nt)-ntar > 0:
            nt1=np.array([i for i in nt0 if i not in ntt])
            if len(nt)-ntar > ntar:
                ntt1=rat.sample(nt1, ntar)
                if len(nt)-ntar -ntar > 0:
                    nt2=np.array([i for i in nt1 if i not in ntt1])
        #print len(nt0)
        #print len(nt2)
        #print len(nt1)
                
        f1=open('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','w')
        for i in range(0, len(ntt)):
            f1.write(str(indf[nt][ntt[i]])+','+str(xp[indf][nt][ntt[i]])+','+str(yp[indf][nt][ntt[i]])+','+str(the[indf][nt][ntt[i]])+','+str(phi[indf][nt][ntt[i]])+','+str(mag_a[indf][nt][ntt[i]])+'\n')
        f1.close()     
        if len(ntt1) > 0:
            f1=open('plate2_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','w')
            for i in range(0, len(ntt1)):
                f1.write(str(indf[nt][ntt1[i]])+','+str(xp[indf][nt][ntt1[i]])+','+str(yp[indf][nt][ntt1[i]])+','+str(the[indf][nt][ntt1[i]])+','+str(phi[indf][nt][ntt1[i]])+'\n')
            f1.close()     
        if len(nt2) > 0:
            f1=open('plate3_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','w')
            for i in range(0, len(ntt2)):
                f1.write(str(indf[nt][ntt2[i]])+','+str(xp[indf][nt][ntt2[i]])+','+str(yp[indf][nt][ntt2[i]])+','+str(the[indf][nt][ntt2[i]])+','+str(phi[indf][nt][ntt2[i]])+'\n')
            f1.close() 
        ind0=indf[nt][ntt]
        ra=the[indf][nt][ntt]
        dec=phi[indf][nt][ntt]
        xf=xp[indf][nt][ntt]
        yf=yp[indf][nt][ntt]
        magf=mag_a[indf][nt][ntt]
        #import matplotlib.pyplot as plt
        #n, bins, patches = plt.hist(magf, 50, density=True, facecolor='g', alpha=0.75)
        #plt.show()
    else:
        f1=open('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','r')
        ind0=[]
        xf=[]
        yf=[]
        ra=[]
        dec=[]
        for line in f1:
            data=line.replace("\n","").split(',')
            data=filter(None,data)
            ind0.extend([np.int(data[0])])
            xf.extend([np.float(data[1])])
            yf.extend([np.float(data[2])])
            ra.extend([np.float(data[3])])
            dec.extend([np.float(data[4])])
        f1.close()
        ind0=np.array(ind0)
        xf=np.array(xf)
        yf=np.array(yf)
        ra=np.array(ra)
        dec=np.array(dec)
    dec_0=phi_o
    ra_0=the_o
        #print red_0[indf][nt]
    print len(ra)
    import matplotlib.pyplot as plt
    plt.plot(ra,dec,'o',markersize=0.5)
    #plt.plot(the[indf][nt],phi[indf][nt],'o',markersize=0.1)
    plt.savefig('test2.pdf',dpi = 1000)
    #plt.show()
    #print np.amin(mag_a[indf][nt][nt0])
    #npt=np.argmin(mag_a[indf][nt][nt0])#,"tes"
    #print indf[nt][nt0][npt]
    plt.close()
    #print indf
    
    #from matplotlib import cm
    #import matplotlib.pyplot as plt
    #fig = plt.figure(figsize=(6,5.5))
    #plt.scatter(obs[1]-z[indf]/ho,obs[2]-y[indf]/ho,c=Fin_mass[indf],vmin=8.5,vmax=11.5,s=10.5)
    #plt.plot([0],[0],'o')
    #plt.plot(obs[1]-z[indf[0]]/ho,obs[2]-y[indf[0]]/ho,'o',markersize=10,color='green') # 
    #plt.plot(obs[0]-66111.5123788,obs[2]-10145.8658918,'o',markersize=5,color='red')
    #plt.savefig('test.pdf',dpi = 1000)
    #plt.show()
    #plt.close()
    
    #sys.exit()
    #id=i#indf[i]
    #id_p='7400-23000' #HaloID[indf[i]]
    #mock_shalo_ill(id,id_p,5,dirtemp=dirtemp,template3=template_1,template5=template_3,template2=template_2,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,plots=plots,observer=obs,ifutype='BOSS',basePath=basePath,expt=900.0)
    return obs,ra_0,dec_0,ra,dec,xf,yf,ind0,ho,Om,Lam
        
def mock_sill_2(fov_p=3.0,rw=1,scp_s=60.4,ra_o=180.0,dec_o=0.0,the_o=0.0,phi_o=0.0,dir1='',basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'): 
    #obs=[-106534,-106534,-106534]
    import random as rat
    obs=[ 52355.0, 43067.0, 43718.0]
    obs=[-106534.0,-106534.0,-106534.0]
    obs=[50000.0,50000.0,-50000.0]
    obs=[0.0,0.0,150000.0]
    plots=1
    f=h5py.File(il.snapshot.snapPath(basePath,135),'r')
    header = dict( f['Header'].attrs.items() )
    f.close()
    Om=header['Omega0']
    Lam=header['OmegaLambda']
    ho=0.74
    HaloID = il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloGrNr'])
    SubHaloMas= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloMassInRadType'])
    Len= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloLenType'])
    Rad=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloHalfmassRad'])#SubhaloStellarPhotometricsRad'])
    XYZ=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloCM'])
    MAGs=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloStellarPhotometrics'])
    bhmass=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloBHMass'])
    Lent=Len[:,0]
    x=XYZ[:,0]
    y=XYZ[:,1]
    z=XYZ[:,2]
    mag=MAGs[:,5]
    #Rad=Rad_T[:,4]
    Fin_mass=np.log10(SubHaloMas[:,4]/ho*1e10+1)
    bhmassT=np.log10(bhmass/ho*1e10+1)
    n_sub=len(HaloID)
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    
    print np.amax(y),np.amin(y)
    xt=x-(np.amax(x)+np.amin(x))/2.0
    yt=y-(np.amax(y)+np.amin(y))/2.0
    zt=z-(np.amax(z)+np.amin(z))/2.0
    print np.amax(xt),np.amin(xt)
    x0=obs[0]#-xt/ho
    y0=obs[1]#-yt/ho
    z0=obs[2]#-zt/ho
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    cam=Rc
    rad=np.sqrt(xt**2.+yt**2.+(cam-zt)**2.)
    reds=reds_cos(rad/1e3)
    print red_0,Rc,'H'
    radA=rad/(1+reds)
    phi=np.arcsin(xt/radA)/10.0#2.5
    the=np.arcsin(yt/(radA*np.cos(phi)))/10.0#2.5
    the=the*180/np.pi+ra_o#*3600
    phi=phi*180/np.pi+dec_o#*3600
    phi_ma=np.amax(phi)
    phi_mi=np.amin(phi)
    the_ma=np.amax(the)
    the_mi=np.amin(the)
    #the_o=ra_o#1.2+180.0
    #phi_o=dec_o#1.0
    if phi_o == 0.0 or the_o == 0.0:
        phi_o=(phi_ma+phi_mi)/2.0
        the_o=(the_ma+the_mi)/2.0
    print the_o,phi_o,phi_ma-phi_o,phi_mi-phi_o
    xp=(the-the_o)*3600.0*scp_s/1e3
    yp=(phi-phi_o)*3600.0*scp_s/1e3
    
    mag_a=mag+5.0*np.log10(rad*1e3)-5.0
    #print mag_a
    
    dm=0.5
    Ms=15.5
    Mi=8.0
    nr=500
    nx=int((Ms-Mi)/dm)
    indf=[]
    if ptt.exists("Selection.dat") == False or rw == 1:
        f=open("Selection.dat","w")
        for i in range(nx-1,-1,-1):#range(0, nx):
            nt=np.where((Fin_mass > (Mi+dm*i)) & (Fin_mass <= (Mi+dm*(i+1))))[0]
            temp=nt#[rat.sample(range(0, len(nt)), nr)]
            for jj in range(0, len(nt)):#nr):
                if mag_a[temp[jj]] < 18.5 and Lent[temp[jj]] > 100 and bhmassT[temp[jj]] > 6.0 and mag_a[temp[jj]] > 11.6 :
                    #print mag_a[temp[jj]]
                    indf.extend([temp[jj]])
                    f.write(str(temp[jj])+'\n')
    else:
        f=open("Selection.dat","r")
        for line in f:
            indf.extend([int(line.replace('\n',''))])
    indf=np.array(indf)
    f.close()
    if ptt.exists('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt') == False or rw == 1:
        nt=np.where(np.sqrt((the[indf]-the_o)**2.0+(phi[indf]-phi_o)**2.0) < fov_p/2.0)[0]
        print len(nt),"H"
            
        f1=open('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','w')
        for i in range(0, len(nt)):
            f1.write(str(indf[nt[i]])+','+str(xp[indf][nt[i]])+','+str(yp[indf][nt[i]])+','+str(the[indf][nt[i]])+','+str(phi[indf][nt[i]])+','+str(mag_a[indf][nt[i]])+'\n')
        f1.close()     
        ind0=indf[nt]
        ra=the[indf][nt]
        dec=phi[indf][nt]
        xf=xp[indf][nt]
        yf=yp[indf][nt]
        magf=mag_a[indf][nt]
        #import matplotlib.pyplot as plt
        #n, bins, patches = plt.hist(magf, 50, density=True, facecolor='g', alpha=0.75)
        #plt.show()
    else:
        f1=open('plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'.txt','r')
        ind0=[]
        xf=[]
        yf=[]
        ra=[]
        dec=[]
        for line in f1:
            data=line.replace("\n","").split(',')
            data=filter(None,data)
            ind0.extend([np.int(data[0])])
            xf.extend([np.float(data[1])])
            yf.extend([np.float(data[2])])
            ra.extend([np.float(data[3])])
            dec.extend([np.float(data[4])])
        f1.close()
        ind0=np.array(ind0)
        xf=np.array(xf)
        yf=np.array(yf)
        ra=np.array(ra)
        dec=np.array(dec)
    dec_0=phi_o
    ra_0=the_o
        #print red_0[indf][nt]
    print len(ra)
    import matplotlib.pyplot as plt
    plt.plot(ra,dec,'o',markersize=0.5)
    #plt.plot(the[indf][nt],phi[indf][nt],'o',markersize=0.1)
    plt.savefig('test2.pdf',dpi = 1000)
    #plt.show()
    #print np.amin(mag_a[indf][nt][nt0])
    #npt=np.argmin(mag_a[indf][nt][nt0])#,"tes"
    #print indf[nt][nt0][npt]
    plt.close()
    return obs,ra_0,dec_0,ra,dec,xf,yf,ind0,ho,Om,Lam

def mock_sill_3(obs=[0.0,0.0,150000.0],fov_p=3.0,rw=1,scp_s=60.4,ra_o=180.0,dec_o=0.0,the_o=0.0,phi_o=0.0,indx=0,dir1='',basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'): 
    import random as rat
    plots=1
    f=h5py.File(il.snapshot.snapPath(basePath,135),'r')
    header = dict( f['Header'].attrs.items() )
    f.close()
    Om=header['Omega0']
    Lam=header['OmegaLambda']
    ho=0.74
    HaloID = il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloGrNr'])
    SubHaloMas= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloMassInRadType'])
    Len= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloLenType'])
    Rad=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloHalfmassRad'])#SubhaloStellarPhotometricsRad'])
    XYZ=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloCM'])
    MAGs=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloStellarPhotometrics'])
    bhmass=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloBHMass'])
    Lent=Len[:,0]
    x=XYZ[:,0]
    y=XYZ[:,1]
    z=XYZ[:,2]
    mag=MAGs[:,5]
    #Rad=Rad_T[:,4]
    Fin_mass=np.log10(SubHaloMas[:,4]/ho*1e10+1)
    bhmassT=np.log10(bhmass/ho*1e10+1)
    n_sub=len(HaloID)
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    #print obs
    #print np.amax(y),np.amin(y),"limits"
    xt=x-(np.amax(x)+np.amin(x))/2.0
    yt=y-(np.amax(y)+np.amin(y))/2.0
    zt=z-(np.amax(z)+np.amin(z))/2.0
    #print np.amax(xt),np.amin(xt)
    #print np.amax(yt),np.amin(yt)
    #print np.amax(zt),np.amin(zt)
    x0=obs[0]#-xt/ho
    y0=obs[1]#-yt/ho
    z0=obs[2]#-zt/ho
    #print np.amax(xt-x0),np.amin(xt-x0),'x_lim'
    #print np.amax(yt-y0),np.amin(yt-y0),'y_lim'
    #print np.amax(zt-z0),np.amin(zt-z0),'z_lim'
    rad=np.sqrt((xt-x0)**2.+(yt-y0)**2.+(zt-z0)**2.)
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    phi=np.arcsin((xt-x0)/radA)/3.5#/10.0#2.5
    the=np.arcsin((yt-y0)/(radA*np.cos(phi)))/3.5#/10.0#2.5
    the=the*180/np.pi+ra_o#*3600
    phi=phi*180/np.pi+dec_o#*3600
    phi_ma=np.amax(phi)
    phi_mi=np.amin(phi)
    the_ma=np.amax(the)
    the_mi=np.amin(the)
    if phi_o == 0.0 or the_o == 0.0:
        phi_o=(phi_ma+phi_mi)/2.0
        the_o=(the_ma+the_mi)/2.0
    #print phi_o,phi_ma-phi_o,phi_mi-phi_o,np.nanmax(reds),'angle limits_A'
    #print the_o,the_ma-the_o,the_mi-the_o,np.nanmax(reds),'angle limits_B'
    xp=(the-the_o)*3600.0*scp_s/1e3
    yp=(phi-phi_o)*3600.0*scp_s/1e3
    
    mag_a=mag+5.0*np.log10(rad*1e3)-5.0

    
    dm=0.5
    Ms=15.5
    Mi=8.0
    nr=500
    nx=int((Ms-Mi)/dm)
    indf=[]
    sycall("mkdir -p Selection")
    if ptt.exists("Selection/Selection_"+str(indx)+".dat") == False or rw == 1:
        f=open("Selection/Selection_"+str(indx)+".dat","w")
        for i in range(nx-1,-1,-1):
            nt=np.where((Fin_mass > (Mi+dm*i)) & (Fin_mass <= (Mi+dm*(i+1))))[0]
            temp=nt
            for jj in range(0, len(nt)):
                if mag_a[temp[jj]] < 18.5 and Lent[temp[jj]] > 100 and bhmassT[temp[jj]] > 0.0 and mag_a[temp[jj]] > 11.6 :#bhmassT[temp[jj]] > 6.0
                    indf.extend([temp[jj]])
                    f.write(str(temp[jj])+' , '+str(HaloID[temp[jj]])+'\n')
                    #f.write(str(temp[jj])+'\n')
    else:
        f=open("Selection/Selection_"+str(indx)+".dat","r")
        for line in f:
            indf.extend([int(line.replace('\n',''))])
    indf=np.array(indf)
    f.close()
    if ptt.exists('Selection/plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'_'+str(indx)+'.txt') == False or rw == 1:
        nt=np.where(np.sqrt((the[indf]-the_o)**2.0+(phi[indf]-phi_o)**2.0) < fov_p/2.0)[0]
       # print len(nt),"H"
            
        f1=open('Selection/plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'_'+str(indx)+'.txt','w')
        for i in range(0, len(nt)):
            f1.write(str(indf[nt[i]])+','+str(xp[indf][nt[i]])+','+str(yp[indf][nt[i]])+','+str(the[indf][nt[i]])+','+str(phi[indf][nt[i]])+','+str(mag_a[indf][nt[i]])+'\n')
        f1.close()     
        ind0=indf[nt]
        ra=the[indf][nt]
        dec=phi[indf][nt]
        xf=xp[indf][nt]
        yf=yp[indf][nt]
        magf=mag_a[indf][nt]
    else:
        f1=open('Selection/plate1_'+str(np.round(the_o,1)).replace('.','')+'_'+str(np.round(phi_o,1)).replace('.','')+'_'+str(indx)+'.txt','r')
        ind0=[]
        xf=[]
        yf=[]
        ra=[]
        dec=[]
        magf=[]
        for line in f1:
            data=line.replace("\n","").split(',')
            data=filter(None,data)
            ind0.extend([np.int(data[0])])
            xf.extend([np.float(data[1])])
            yf.extend([np.float(data[2])])
            ra.extend([np.float(data[3])])
            dec.extend([np.float(data[4])])
            magf.extend([np.float(data[5])])
        f1.close()
        ind0=np.array(ind0)
        xf=np.array(xf)
        yf=np.array(yf)
        ra=np.array(ra)
        dec=np.array(dec)
        magf=np.array(magf)
    dec_0=phi_o
    ra_0=the_o
    obs_f=np.zeros([3,len(ra)])
    for i in range(0, len(ra)):
        obs_f[:,i]=obs
    #print len(ra)
    import matplotlib.pyplot as plt
    plt.plot(ra,dec,'o',markersize=0.5)
    plt.savefig('Selection/test_'+str(indx)+'.pdf',dpi = 1000)
    plt.close()
    return obs_f,ra_0,dec_0,ra,dec,xf,yf,ind0,ho,Om,Lam

def cone_mock_ill(fov_p=3.0,rw=1,scp_s=60.4,ra_o=180.0,dec_o=0.0,the_o=0.0,phi_o=0.0,dir1='',basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'):
    lt=75000.0
    do=lt/2.0#+100.0  
    nid=5
    cont=0
    for i in range(0, nid):
        for j in [0,1,-1]:
            for k in [0,1,-1]:
                obs_i=[0.0+lt*j,0.0+lt*k,0.0+lt*i+do]          
                obs1,ra_01,dec_01,raf1,decf1,xf1,yf1,ind01,ho1,Om1,Lam1=mock_sill_3(obs=obs_i,fov_p=fov_p,rw=rw,dir1=dir1,basePath=basePath,scp_s=scp_s,ra_o=ra_o,dec_o=dec_o,the_o=the_o,phi_o=phi_o,indx=cont)
                if i == 0 and j == 0 and k == 0:
                    #obs=np.zeros([3,1])
                    #obs[:,0]=obs1
                    obs=obs1
                    ra_0=ra_01
                    dec_0=dec_01
                    raf=raf1
                    decf=decf1
                    xf=xf1
                    yf=yf1
                    ind0=ind01
                    ho=ho1
                    Om=Om1
                    Lam=Lam1
                else:
                    if len(raf1) >= 1:
                        #obst=np.zeros([3,1])
                        #obst[:,0]=obs1
                        obs=np.concatenate((obs,obs1),axis=1)
                        raf=np.concatenate((raf,raf1))
                        decf=np.concatenate((decf,decf1))
                        xf=np.concatenate((xf,xf1))
                        yf=np.concatenate((yf,yf1))
                        ind0=np.concatenate((ind0,ind01))
                        print obs.shape
                #if len(raf1) >= 1:        
                print len(raf1),i,j,k,cont,'slide'
                cont=cont+1
    print len(raf)
    #print obs.shape
    fco=open('Selection_sc.csv','w')
    fco.write('RA,DEC\n')
    for i in range(0, len(raf)):
        fco.write(str(raf[i])+','+str(decf[i])+' \n')
    fco.close()
    import matplotlib.pyplot as plt
    plt.plot(raf,decf,'o',markersize=0.5)
    plt.savefig('test_cone_fin.pdf',dpi = 1000)
    plt.close()
    return obs,ra_0,dec_0,raf,decf,xf,yf,ind0,ho,Om,Lam           
 
    
=======
    the=the*180/np.pi*3600#+ran.randn(len(rad))*2.0
    phi=phi*180/np.pi*3600#+ran.randn(len(rad))*2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600#+ran.randn(len(rad_g))*2.0
    phi_g=phi_g*180/np.pi*3600#+ran.randn(len(rad_g))*2.0
    ns=3*nl*(nl-1)+1
    Dfib=fibA*scalep
    Rifu=Dfib*((2.0*nl-1.0)/2.0-0.5)
    xfib0=-Rifu
    yfib0=0
    dxf=1.0
    dyf=np.sin(60.*np.pi/180.)
    xifu=np.zeros(ns)
    yifu=np.zeros(ns)
    ini=0
    for i in range(0, nl):
        nt=nl*2-1-i
        yfib=yfib0+i*dyf*Dfib
        for j in range(0, nt):
            xfib=xfib0+(j*dxf+0.5*i)*Dfib
            xifu[ini]=xfib
            yifu[ini]=yfib
            ini=ini+1
        if i > 0:
            for j in range(0, nt):
                xfib=xfib0+(j*dxf+0.5*i)*Dfib
                xifu[ini]=xfib
                yifu[ini]=-yfib
                ini=ini+1
    ndt=35
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.0]
    dit[1,:]=[+0.00,+1.0]
    dit[2,:]=[+0.00,-1.0]
    dit[3,:]=[+0.25,+0.5]
    dit[4,:]=[+0.25,-0.5]
    dit[5,:]=[-0.25,+0.5]
    dit[6,:]=[-0.25,-0.5]
    dit[7,:]=[+0.50,+0.5]
    dit[8,:]=[+0.50,-0.5]
    dit[9,:]=[+0.50,+1.5]
    dit[10,:]=[+0.50,-1.5]
    dit[11,:]=[-0.50,+0.5]
    dit[12,:]=[-0.50,-0.5]
    dit[13,:]=[-0.50,+1.5]
    dit[14,:]=[-0.50,-1.5]
    dit[15,:]=[+1.00,+0.0]
    dit[16,:]=[+1.00,+0.5]
    dit[17,:]=[+1.00,-0.5]
    dit[18,:]=[+1.00,+1.0]
    dit[19,:]=[+1.00,-1.0]    
    dit[20,:]=[-1.00,+1.0]
    dit[21,:]=[-1.00,-1.0]
    dit[22,:]=[+1.50,+0.0]
    dit[23,:]=[+1.50,-0.5]
    dit[24,:]=[+1.50,+0.5]
    dit[25,:]=[+0.50,+0.0]
    dit[26,:]=[-0.50,+0.0]
    dit[27,:]=[+0.35,-0.9]
    dit[28,:]=[-0.35,-0.9]
    dit[29,:]=[+0.85,-1.3]
    dit[30,:]=[-0.85,-1.3]
    dit[31,:]=[+0.78,-0.25]
    dit[32,:]=[+0.78,+0.25]
    dit[33,:]=[-0.35,+0.73]
    dit[34,:]=[+0.50,+0.70]
    dyf=1.0
    dyt=Dfib/2.0/np.cos(30.0*np.pi/180.0)
    dxt=Dfib/2.0
    ndt=3
    dit=np.zeros([ndt,2])
    dit[0,:]=[+0.00,+0.00]+ran.randn(2)*0.025
    dit[1,:]=[+0.00,+dyt/1.0]+ran.randn(2)*0.025
    dit[2,:]=[-dxt, +dyt/2.0]+ran.randn(2)*0.025
#    ndt=35
#    print Dfib/2.0/np.cos(30.0*np.pi/180.0)
#    import matplotlib.pyplot as plt
#    matplotlib.use('Agg')
    #plt.plot(phi,the,'o')
#    for i in range(0, ndt):
#        plt.plot(xifu+dit[i,0],yifu+dyf*dit[i,1],'o')
#    plt.show()
#    sys.exit()
#    Av=1.1
    
    #template="../../Base_bc03/templete_bc03.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template5)
    [ssp_template3,wave3,age_ssp3,met_ssp3,ml_ssp3,crval_w3,cdelt_w3,crpix_w3]=ssp_extract(template3)
#    [ssp_template5,wave5,age_ssp5,met_ssp5,ml_ssp5,crval_w5,cdelt_w5,crpix_w5]=ssp_extract(template5)
    ml_ssp=1.0/ml_ssp
    #ml_ssp5=1.0/ml_ssp5
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    #pht_g=asosiate_pho2(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_s,met_s,age_s,x_g,y_g,z_g,x,y,z,Rs)
    #print np.amax(pht_g[np.where(np.log10(pht_g) != 0)[0]]),"Q(H) max"
    #print np.amin(pht_g[np.where(np.log10(pht_g) != 0)[0]]),"Q(H) min"
    #print np.amax(met_g)/0.02,"Z_g max"
    #print np.amin(met_g)/0.02,"Z_g min"
    #print np.amax(temp_g),"T max"
    #print np.amin(temp_g),"T min"
    #print np.amax(nh),"nh max"
    #print np.amin(nh),"nh min"
    #print np.amax(met_s),"Z_s max"
    #print np.amin(met_s),"Z_s min"
    #print np.amax(met_ssp),"Z_ssp max"
    #print np.amin(met_ssp),"Z_ssp min"
    #print tem_gas[0]
    #print den_gas[0]
    #print met_gas[0]
    #print pht_gas[0]
    #sys.exit()
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=sp_samp
    crval_w=3622.0
    crpix_w=1
    wave_f=np.arange(crval_w,10352.0,cdelt_w)
    if Flux_m != 20.0:
        Flux_m=1.0/(10.0**(-0.4*Flux_m)*np.pi*(fibB*scalep/2.0)**2.0*(466.9e-11*cdelt_w)/1e-16*SNi)
    else:
        Flux_m=20.0
    #print Flux_m
    #sys.exit()
    s_nr=noise_sig(wave_f,SNi-1.0)#15.0#35#55....70,105#140#80#30#15#60.0
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    spec_ifu=np.zeros([nw,ndt*ns])
    spec_ifu_e=np.zeros([nw,ndt*ns])
    spec_val=np.zeros([35,ndt*ns])
    n_ages=num_ages(age_ssp3)
    ages_r=arg_ages(age_ssp3)
    sim_imag=np.zeros([n_ages,ndt*ns])
    x_ifu=np.zeros(ndt*ns)
    y_ifu=np.zeros(ndt*ns)
    facto=(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)#*np.pi
    t_noise=1.0/Flux_m#0.007#2.0
    con=0
    for i in range(0, ndt):
        #print i
        phie=phi+ran.randn(len(rad))*seeing# 1.43#/2.0
        thee=the+ran.randn(len(rad))*seeing# 1.43#/2.0
        phieg=phi_g+ran.randn(len(rad_g))*seeing #1.43#/2.0
        theeg=the_g+ran.randn(len(rad_g))*seeing #1.43#/2.0
        for j in range(0, ns):
            #print i, j
            xo=xifu[j]+dit[i,0]
            yo=yifu[j]+dyf*dit[i,1]    
            r=np.sqrt((xo-phie)**2.0+(yo-thee)**2.0)
            r_g=np.sqrt((xo-phieg)**2.0+(yo-theeg)**2.0)
            #print np.amin(phi),np.amax(phi),np.amin(the),np.amax(the),Rifu
#            print np.amin(r),fibB*scalep/2.0,fibB,scalep
            nt=np.where(r <= fibB*scalep/2.0)[0]
            nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            noise=t_noise*ran.randn(nw)/s_nr#*0.0
            mass_t=0
            ml_t=0
            vel_t=0
            sfr_t=0
            Av_s=0
            Av_sg=0
            sve_t=0
            Lt=0
            Ltg=0
            ml_ts=0
            met_ligt=0
            met_mas=0
            age_ligt=0
            age_mas=0
            Av_ligt=0
            Av_flux=0
            Ve_ligt=0
            Ve_flux=0
            Avg_ligt=0
            Avg_flux=0
            Veg_ligt=0
            Veg_flux=0
            Ft=0
            Ftg=0
            Sig_flux=0
            Sig_ligt=0
            Sig_flux_g=0
            Sig_ligt_g=0
            wl_t=[]
            wf_t=[]
            wl_tg=[]
            wf_tg=[]
            va_1=[]
            va_1g=[]
            Mft=0
            age_flux=0
            age_Mflux=0
            met_ligt_g=0
            met_flux_g=0
#            noise2=t_noise/s_nr
#            plt.plot(wave_f,np.abs(noise*16.0))
#            plt.plot(wave_f,np.abs(noise2*16.0))
#            plt.show()
#            sys.exit()
            #print i,j,len(nt),"STARS"
            sycall('echo '+str(i)+'  '+str(j)+'  '+str(len(nt))+'  STARS')
            if len(nt) > 0:
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                sve_t=np.std(v_rad[nt])
                mass_t_t=asosiate_ages(age_ssp3,age_s[nt],mass_s[nt])
                sim_imag[:,con]=mass_t_t*facto
                for k in range(0, len(nt)):
                    #r_ge=np.sqrt((phi[nt[k]]-phi_g)**2.0+(the[nt[k]]-the_g)**2.0)
                    #nt_e=np.where((r_ge <= 0.05) & (rad_g <= rad[nt[k]]))[0]
                    #nt_e=np.where((r_g <= fibB*scalep/2.0) & (rad_g <= rad[nt[k]]))[0]
                    nt_e=np.where((abs(phi[nt[k]]-phi_g) <= d_r) & (abs(the[nt[k]]-the_g) <= d_r) & (rad_g <= rad[nt[k]]))[0]#DECOMENTAR
                    if len(nt_e) > 0:#DECOMENTAR
                        Av=np.sum(Av_g[nt_e])#DECOMENTAR
                    else:#DECOMENTAR
                        Av=0#DECOMENTAR
                    #Av=0#COMENTAR
                    Av_s=10**(-0.4*Av)+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        #print in_ssp[nt[k]]
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust/1e-16#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]#*20.0
                            Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]#/20.0
                            Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                            met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                            met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                            age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                            age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas                            
                            #if Av > 0:
                            ft_w=mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)
                            lt_w=mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]
                            fm_w=mass_s[nt[k]]*10**(-0.4*Av*0.44)
                            Mft=Mft+fm_w
                            age_flux=np.log10(age_s[nt[k]])*ft_w+age_flux
                            age_Mflux=np.log10(age_s[nt[k]])*fm_w+age_Mflux
                            Ve_ligt=v_rad[nt[k]]*lt_w+Ve_ligt
                            Ve_flux=v_rad[nt[k]]*ft_w+Ve_flux
                            Av_ligt=10**(-0.4*Av)*lt_w+Av_ligt
                            Av_flux=10**(-0.4*Av)*ft_w+Av_flux
                            va_1.extend([v_rad[nt[k]]])
                            wf_t.extend([ft_w])
                            wl_t.extend([lt_w])
                            #plt.plot(wave,ssp_template[in_ssp[nt[k]],:])
                            #plt.show()
                           # plt.plot(wave,spect_t)
                            #plt.show()
                            #print radL[nt[k]],mass_s[nt[k]],age_s[nt[k]],met_s[nt[k]]
                            #sys.exit()
                ml_t=ml_t/len(nt)
                ml_ts=mass_t/Lt
                if Lt > 0:
                    #met_ligt=10.0**(met_ligt/Lt)
                    met_ligt=met_ligt/Lt
                    age_ligt=10.0**(age_ligt/Lt)
                    age_flux=10.0**(age_flux/Ft)
                    age_Mflux=10.0**(age_Mflux/Mft)
                    Av_ligt=(Av_ligt/Lt)
                    Av_flux=(Av_flux/Ft)
                    Ve_ligt=(Ve_ligt/Lt)
                    Ve_flux=(Ve_flux/Ft)
                    #met_mas=10.0**(met_mas/mass_t)
                    met_mas=met_mas/mass_t
                    age_mas=10.0**(age_mas/mass_t)
                    va_1=np.array(va_1)
                    wf_t=np.array(wf_t)
                    wl_t=np.array(wl_t)
                    Sig_flux=np.sqrt(np.nansum(np.abs(wf_t)*(Ve_flux-va_1)**2.0)/(np.nansum(np.abs(wf_t))-np.nansum(wf_t**2.0)/np.nansum(np.abs(wf_t))))
                    Sig_ligt=np.sqrt(np.nansum(np.abs(wl_t)*(Ve_ligt-va_1)**2.0)/(np.nansum(np.abs(wl_t))-np.nansum(wl_t**2.0)/np.nansum(np.abs(wl_t))))
                spect_t[np.isnan(spect_t)]=0
                spect_i=inst_disp(wave,spect_t,sigma_inst)
                spect_ii=spec_resol(wave,spect_i,sp_res)
                spect=interp1d(wave,spect_ii,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
                spec_val[0,con]=Av_s/len(nt)#*facto
#                if j == 5:
#                    plt.plot(wave_f,spect*facto)
#                    plt.plot(wave_f,noise*facto)
#                    plt.savefig(dir_o+outf+'test5.pdf')
#                    plt.show()
#                    sys.exit()
            #print i,j,len(nt_g),"GAS"
            sycall('echo '+str(i)+'  '+str(j)+'  '+str(len(nt_g))+'  GAS')
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                #Ha_f=0
                #SFH_a=np.sum(sfri[nt_g])
                #Tem_a=np.average(temp_g[nt_g])
                #phi_a=np.average(pht_g[nt_g])
                #print SFH_a/(7.9e-42)
                for k in range(0, len(nt_g)):
                    if band_g[nt_g[k]] > 0:
                        #r_ge=np.sqrt((phi_g[nt_g[k]]-phi_g)**2.0+(the_g[nt_g[k]]-the_g)**2.0)
                        #nt_e=np.where((r_ge <= 0.05) & (rad_g <= rad_g[nt_g[k]]))[0]
                        #nt_e=np.where((r_g <= fibB*scalep/2.0) & (rad_g <= rad_g[nt_g[k]]))[0]
                        nt_e=np.where((abs(phi_g[nt_g[k]]-phi_g) <= d_r) & (abs(the_g[nt_g[k]]-the_g) <= d_r) & (rad_g <= rad_g[nt_g[k]]))[0]#DECOMENTAR
                        if len(nt_e) > 0:#DECOMENTAR
                            Av=np.sum(Av_g[nt_e])#DECOMENTAR
                        else:#DECOMENTAR
                            Av=0#DECOMENTAR
                    #    Av=0#COMENTAR
                        Av_sg=Av+Av_sg
                        if np.isnan(in_gas[nt_g[k]]):
                            spect_gf=spect_gf
                        else:
                            if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:
                                dust=10**(-0.4*Av*dust_rat_gas)
                                spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#*0.01#*mass_g[nt_g[k]]
                                spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                                lt_wg=np.nansum(gas_template[in_gas[nt_g[k]],:])*10.0**(-3.18)
                                ft_wg=np.nansum(spect_sfg)
                                Ltg=Ltg+lt_wg
                                Ftg=Ftg+ft_wg
                                Veg_ligt=v_rad_g[nt_g[k]]*lt_wg+Veg_ligt
                                Veg_flux=v_rad_g[nt_g[k]]*ft_wg+Veg_flux
                                Avg_ligt=10**(-0.4*Av)*lt_wg+Avg_ligt
                                Avg_flux=10**(-0.4*Av)*ft_wg+Avg_flux
                                met_ligt_g=np.log10(met_g[nt_g[k]])*lt_wg+met_ligt_g   
                                met_flux_g=np.log10(met_g[nt_g[k]])*ft_wg+met_flux_g                    
                                spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                                spect_g=spect_sfg+spect_g    
                                va_1g.extend([v_rad_g[nt_g[k]]])
                                wf_tg.extend([ft_wg])
                                wl_tg.extend([lt_wg])                            
                                #plt.xlim(3400,4000)
                            #plt.plot(wave_g,dust)
                            #plt.show()
                            #plt.xlim(3400,4000)                            
                            #plt.show()
                            #plt.plot(wave_g,gas_template[in_gas[nt_g[k]],:])
                            #plt.xlim(3400,4000)
                            #plt.plot(wave_g,spect_sg)
                            #plt.show()
                            #plt.xlim(3400,4000)
                            #plt.plot(wave_g,spect_g)
                            #plt.show()
                            #sys.exit()
                            #temp=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*mass_g[nt_g[k]]*3.846e33
                            #nt=np.where((wave_g < (6562+170)) & (wave_g > (6562-170)))[0]
                            #nt1=np.where((wave_g < (6262+170)) & (wave_g > (6262-170)))[0]
                            #Ha_f=np.sum(temp[nt])-np.sum(temp[nt1])+Ha_f
                #print Ha_f,(SFH_a/(7.9e-42)/Ha_f)**(1/3.)
                #print Tem_a,phi_a
#                plt.xlim(3400,4000)
#                plt.plot(wave_g,spect_g)   
                if Ltg > 0:
                    Avg_ligt=(Avg_ligt/Ltg)
                    Avg_flux=(Avg_flux/Ftg)
                    Veg_ligt=(Veg_ligt/Ltg)
                    Veg_flux=(Veg_flux/Ftg)
                    #met_ligt_g=10.0**(met_ligt_g/Ltg)
                    #met_flux_g=10.0**(met_flux_g/Ftg)
                    met_ligt_g=met_ligt_g/Ltg
                    met_flux_g=met_flux_g/Ftg
                    va_1g=np.array(va_1g)
                    wf_tg=np.array(wf_tg)
                    wl_tg=np.array(wl_tg)
                    Sig_flux_g=np.sqrt(np.nansum(np.abs(wf_tg)*(Veg_flux-va_1g)**2.0)/(np.nansum(np.abs(wf_tg))-np.nansum(wf_tg**2.0)/np.nansum(np.abs(wf_tg))))
                    Sig_ligt_g=np.sqrt(np.nansum(np.abs(wl_tg)*(Veg_ligt-va_1g)**2.0)/(np.nansum(np.abs(wl_tg))-np.nansum(wl_tg**2.0)/np.nansum(np.abs(wl_tg))))    
                spec_val[6,con]=np.sum(Av_g[nt_g])         
                spec_val[4,con]=Av_sg/len(nt_g)
                spect_g[np.isnan(spect_g)]=0
                spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
                spect_g_ii=spec_resol(wave_g,spect_g_i,sp_res)
                spect_gf=interp1d(wave_g,spect_g_ii,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
                #plt.plot(wave_f,spect_gf+spect)
                #plt.show()
                #sys.exit()
#            spec_ifu[:,con]=(spect_gf+spect)+noise#np.median(spect_gf+spect)*noise
#            spec_ifu[:,con]=np.sqrt((spect_gf+spect)**2.0+noise**2.0)            
            #spec_ifu[:,con]=np.sqrt((facto**2.0)*(spect_gf+spect)**2.0+noise**2.0)
            spec_ifu[:,con]=facto*(spect_gf+spect+noise)
            spec_ifu_e[:,con]=noise*facto
            spec_val[1,con]=mass_t*facto
            spec_val[2,con]=vel_t#*facto
            spec_val[3,con]=sfr_t*facto
            spec_val[5,con]=spec_val[0,con]*0+spec_val[4,con]
            spec_val[7,con]=sve_t
            spec_val[8,con]=ml_t#*facto
            spec_val[10,con]=Lt*facto
            spec_val[9,con]=ml_ts
            spec_val[11,con]=met_ligt
            spec_val[12,con]=met_mas
            spec_val[13,con]=age_ligt
            spec_val[14,con]=age_mas
            spec_val[15,con]=Av_ligt
            spec_val[16,con]=Ft*facto
            spec_val[17,con]=Av_flux
            spec_val[18,con]=Ve_ligt
            spec_val[19,con]=Ve_flux
            spec_val[20,con]=Sig_ligt
            spec_val[21,con]=Sig_flux
            spec_val[22,con]=Avg_ligt
            spec_val[23,con]=Avg_flux
            spec_val[24,con]=Veg_ligt
            spec_val[25,con]=Veg_flux
            spec_val[26,con]=Sig_ligt_g
            spec_val[27,con]=Sig_flux_g
            spec_val[28,con]=Ftg*facto
            spec_val[29,con]=Ltg*facto
            spec_val[30,con]=Mft*facto
            spec_val[31,con]=age_flux
            spec_val[32,con]=age_Mflux
            spec_val[33,con]=met_ligt_g
            spec_val[34,con]=met_flux_g
            #+ran.randn(nw)*np.median(spect)*0.05
            x_ifu[con]=xo
            y_ifu[con]=yo
            con=con+1
    nl=int(round((np.amax([np.amax(x_ifu),-np.amin(x_ifu),np.amax(y_ifu),-np.amin(y_ifu)])+1)*2/pix_s))
    print nl
    ifu=np.zeros([nw,nl,nl])
    ifu_e=np.ones([nw,nl,nl])
    ifu_1=np.ones([nw,nl,nl])
    ifu_m=np.zeros([nw,nl,nl])
    ifu_v=np.zeros([35,nl,nl])
    ifu_a=np.zeros([n_ages,nl,nl])
    xo=-nl/2*pix_s
    yo=-nl/2*pix_s
    print xo,yo
    xi=xo
    xf=xo
#    import matplotlib.pyplot as plt
                #plt.show()
    int_spect=np.zeros(nw)
    for i in range(0, nl):
        xi=xf
        xf=xf+pix_s
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+pix_s
            spt_new=np.zeros(nw)
            spt_err=np.zeros(nw)
            spt_val=np.zeros(35)
            spt_mas=np.zeros(n_ages)
            Wgt=0
            for k in range(0, len(x_ifu)):
                V1=np.sqrt((x_ifu[k]-xi)**2.0+(y_ifu[k]-yf)**2.0)
                V2=np.sqrt((x_ifu[k]-xf)**2.0+(y_ifu[k]-yf)**2.0)
                V3=np.sqrt((x_ifu[k]-xi)**2.0+(y_ifu[k]-yi)**2.0)
                V4=np.sqrt((x_ifu[k]-xf)**2.0+(y_ifu[k]-yi)**2.0)
                Vt=np.array([V1,V2,V3,V4])
                Rsp=np.sqrt((x_ifu[k]-(xf+xi)/2.0)**2.0+(y_ifu[k]-(yf+yi)/2.0)**2.0)
                Vmin=np.amin(Vt)
                Vmax=np.amax(Vt)
                #if Vmin <= fibB*scalep/2.0:
                #    if Vmax <= fibB*scalep/2.0:
                #        Wg=(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)
                #    else:
                #        Wg=(1.0-(Vmax-fibB*scalep/2.0)/(np.sqrt(2.0)*pix_s))*(pix_s)**2.0/(np.pi*(fibB*scalep/2.0)**2.0)
                #    spt_new=spec_ifu[:,k]*Wg+spt_new
                #    spt_err=(spec_ifu_e[:,k]*Wg)**2.0+spt_err**2.0
                #    Wgt=Wgt+Wg
                if Rsp <= fibA*scalep*1.4/2.0:
                    Wg=np.exp(-(Rsp/pix_s)**2.0/2.0)
                    spt_new=spec_ifu[:,k]*Wg+spt_new
                    spt_err=(spec_ifu_e[:,k]*Wg)**2.0+spt_err**2.0
                    spt_val=spec_val[:,k]*Wg+spt_val
                    spt_mas=sim_imag[:,k]*Wg+spt_mas
                    Wgt=Wgt+Wg
            if Wgt == 0:
                Wgt=1
            ifu[:,j,i]=spt_new/Wgt
            ifu_v[:,j,i]=spt_val/Wgt
            ifu_a[:,j,i]=spt_mas/Wgt
            #ifu_imag[:,j,i]=spt_imag/Wgt
            if np.sum(np.sqrt(spt_err/Wgt**2.0)) == 0:
                ifu_e[:,j,i]=1.0
            else:
                ifu_e[:,j,i]=np.sqrt(spt_err/Wgt**2.0)
               # ifu_m[:,j,i]=1.0
            int_spect=int_spect+spt_new/Wgt
                #plt.plot(wave,int_spect)
    #plt.plot(wave_f,int_spect)
    ##plt.show() 
    ##fig.tight_layout()
    #plt.savefig(dir_o+outf+'_int.pdf')
    #plt.close()
             
    h1=pyf.PrimaryHDU(ifu)#.header
    h2=pyf.ImageHDU(ifu_e)
    h3=pyf.ImageHDU(ifu_1)
    h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mock "+ifutype+" IFU"
    h["CRVAL1"]=0
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['R']=(sp_res,'Spectral Resolution')
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['IFUCON']=(str(np.int(ns))+' ','NFibers')
    h['UNITS']='1E-16 erg/s/cm^2'
    hlist=pyf.HDUList([h1,h2,h3,h4])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip  '+out_fit1)
    ifu_v[0,:,:]=-2.5*np.log10(ifu_v[0,:,:]+0.0001)
    ifu_v[1,:,:]=np.log10(ifu_v[1,:,:]+1.0)
    ifu_v[30,:,:]=np.log10(ifu_v[30,:,:]+1.0)
    ifu_v[10,:,:]=np.log10(ifu_v[10,:,:]+1.0)
    ifu_v[29,:,:]=np.log10(ifu_v[29,:,:]+1.0)
    ifu_v[15,:,:]=-2.5*np.log10(ifu_v[15,:,:]+0.0001)
    ifu_v[17,:,:]=-2.5*np.log10(ifu_v[17,:,:]+0.0001)
    ifu_v[22,:,:]=-2.5*np.log10(ifu_v[22,:,:]+0.0001)
    ifu_v[23,:,:]=-2.5*np.log10(ifu_v[23,:,:]+0.0001)
    h1t=pyf.PrimaryHDU(ifu_v)
    h=h1t.header
    h["NAXIS"]=3
    h["NAXIS3"]=35
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Real Values "+ifutype+" IFU"
    h["CRVAL1"]=0
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['Type0']=('Av_T    ','Mag')
    h['Type1']=('MASS    ','log10(Msun)')
    h['Type2']=('VEL     ','km/s')
    h['Type3']=('SFR     ','Msun/yr')
    h['Type4']=('DUST_G  ','Av BETA')
    h['Type5']=('DUST_T  ','Av BETA')
    h['Type6']=('DUST_Av ','Av BETA')
    h['Type7']=('DISP    ','km/s')
    h['Type8']=('aML     ','Msun/Lsun BETA')
    h['Type9']= ('tML     ','Msun/Lsun BETA')
    h['Type10']=('LUM     ','log10(Lsun)')
    h['Type11']=('Z_lw    ','log10(Z/H) add 1.77 to convert to log10(Z/Z_sun)')
    h['Type12']=('Z_mw    ','log10(Z/H) add 1.77 to convert to log10(Z/Z_sun)')
    h['Type13']=('AGE_lw  ','Gyr')
    h['Type14']=('AGE_mw  ','Gyr')
    h['Type15']=('Av_lw   ','Mag')
    h['Type16']=('FLUX    ','1e-16 ergs/s/cm2')
    h['Type17']=('Av_fw   ','Mag')
    h['Type18']=('VEL_lw  ','km/s')
    h['Type19']=('VEL_fw  ','km/s')
    h['Type20']=('DIS_lw   ','km/s BETA')
    h['Type21']=('DIS_fw   ','km/s BETA')
    h['Type22']=('Av_lw_g  ','Mag')
    h['Type23']=('Av_fw_g  ','Mag')
    h['Type24']=('VEL_lw_g ','km/s')
    h['Type25']=('VEL_fw_g ','km/s')
    h['Type26']=('DIS_l_gas','km/s BETA')
    h['Type27']=('DIS_f_gas','km/s BETA')
    h['Type28']=('FLUX_gas','1e-16 ergs/s/cm2 Bolometric')
    h['Type29']=('LUM_gas ','log10(Lsun) Bolometric')
    h['Type30']=('MASS_fw ','log10(Msun) BETA')
    h['Type31']=('AGE_fw  ','Gyr')
    h['Type32']=('AGE_mfw ','Gyr BETA')
    h['Type33']=('Z_lw_gas ','log10(Z/H) add 10.46 to convert to 12+log10(O/H) or add 1.77 to convert to log10(Z/Z_sun)')
    h['Type34']=('Z_fw_gas ','log10(Z/H) add 10.46 to convert to 12+log10(O/H) or add 1.77 to convert to log10(Z/Z_sun)')#8.69 
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['R']=(sp_res,'Spectral Resolution')
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist1=pyf.HDUList([h1t])
    hlist1.update_extend()
    out_fit=dir_o+outf+'_val.fits'
    wfits_ext(out_fit,hlist1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val.fits'
    sycall('gzip -f '+out_fit1)
    h1tt=pyf.PrimaryHDU(ifu_a)
    h=h1tt.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Real Values "+ifutype+" IFU"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*pix_s/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*pix_s/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=seeing
    h['FOV']=Rifu*2.0
    h['CAMX']=0
    h['CAMY']=0
    h['CAMZ']=cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist1t=pyf.HDUList([h1tt])
    hlist1t.update_extend()
    out_fit=dir_o+outf+'_val_mass_t.fits'
    wfits_ext(out_fit,hlist1t)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_val_mass_t.fits'
    sycall('gzip -f '+out_fit1)
    
def photo_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    seeing=sig    
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    #radL=10.0*3.08567758e16*100.0
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600+ran.randn(len(rad))*seeing#/2.0
    phi=phi*180/np.pi*3600+ran.randn(len(rad))*seeing#/2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600+ran.randn(len(rad_g))*seeing#/2.0
    phi_g=phi_g*180/np.pi*3600+ran.randn(len(rad_g))*seeing#/2.0
#    print np.amax(phi),np.amin(np.abs(phi))
#    print np.amax(x),np.amin(np.abs(x))
#    sys.exit()
    Av=1.1
    #template="../home/sanchez/ppak/legacy/gsd61_156.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=30.0#5.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)#25000
    s_nr=noise_sig(wave_f,25)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    #spec_ifu=np.zeros([nw,ndt*ns])
    #spec_ifu_e=np.zeros([nw,ndt*ns])
    #x_ifu=np.zeros(ndt*ns)
    #y_ifu=np.zeros(ndt*ns)
    #t_noise=2.0
    #con=0
    #for i in range(0, ndt):
    #    for j in range(0, ns):
    #        xo=xifu[j]+dit[i,0]
    #        yo=yifu[j]+dyf*dit[i,1]    
    #        r=np.sqrt((xo-phi)**2.0+(yo-the)**2.0)
    #        r_g=np.sqrt((xo-phi_g)**2.0+(yo-the_g)**2.0)
    #        nt=np.where(r <= fibB*scalep/2.0)[0]
    #        nt_g=np.where(r_g <= fibB*scalep/2.0)[0]
    #        spect_t=np.zeros(nw_s)
    #        spect=np.zeros(nw)
    #        spect_g=np.zeros(nw_g)
    #        spect_gf=np.zeros(nw)
    #        noise=t_noise*ran.randn(nw)/s_nr
            #if len(nt) > 0:
                #for k in range(0, len(nt)):
                #    if np.isnan(in_ssp[nt[k]]):
                #        spect=spect
                #    else:
                #        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                #            #dust=10**(-0.4*Av*dust_rat_ssp)
                #            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*1e10*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16#*dust
                #            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                #            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                #            spect_t=spect_sf+spect_t
                #spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                #spect[np.isnan(spect)]=0
            #if len(nt_g) > 0:
            #    for k in range(0, len(nt_g)):
            #        if np.isnan(in_gas[nt_g[k]]):
            #            spect_gf=spect_gf
            #        else:
            #            if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 157:
            #                dust=10**(-0.4*Av*dust_rat_gas)
            #                spect_sg=gas_template[in_gas[nt_g[k]],:]*10**(ha_gas[in_gas[nt_g[k]]])*sup[nt_g[k]]*fact[nt_g[k]]*band_g[nt_g[k]]/1e1/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust/1e-16
            #                spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
            #                spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
            #                spect_g=spect_sfg+spect_g
            #    spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
            #    spect_gf[np.isnan(spect_gf)]=0
            #spec_ifu[:,con]=(spect_gf+spect)+noise
            #spec_ifu_e[:,con]=noise
            #x_ifu[con]=xo
            #y_ifu[con]=yo
            #con=con+1
    #import matplotlib.pyplot as plt
    #plt.plot(phi,the,'o')
    #plt.ylabel('some numbers')
    #plt.show()
    photo_imag=np.zeros([nw,nl,nl])
    ext_temp=np.zeros([3,nl,nl])
    #mas_temp=np.zeros([nl,nl])
    #vel_temp=np.zeros([nl,nl])
    #sfr_temp=np.zeros([nl,nl])
    #ra1=ran.randn(len(phi))*leng_s/dkpcs
    #de1=ran.randn(len(the))*leng_s/dkpcs
    #ra2=ran.randn(len(phi_g))*leng_s/dkpcs
    #de2=ran.randn(len(the_g))*leng_s/dkpcs
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            nt_g=np.where((phi_gt < xf) & (phi_gt >= xi) & (the_gt < yf) & (the_gt >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            #noise=t_noise*ran.randn(nw)/s_nr
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                for k in range(0, len(nt)):
                    nt_e=np.where(rad_g[nt_g] <= rad[nt[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
#                    Av=0
                    Av_s=Av+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            #print spect_sf
                            #print dlam[nt[k]]
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            #plt.plot(wave,ssp_template[in_ssp[nt[k]],:])
                            #plt.show()
                            #plt.plot(wave,spect_sf)
                            #plt.show()
                            #sys.exit()
                #ext_temp[0,j,i]=Av_s/len(nt)
                ext_temp[0,j,i]=mass_t
                ext_temp[1,j,i]=vel_t
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                for k in range(0, len(nt_g)):
#                    nt_e=np.where((phi_g < xf) & (phi_g >= xi) & (the_g < yf) & (the_g >= yi) & (rad_g <= rad_g[nt_g[k]]))[0]
                    nt_e=np.where(rad_g[nt_g] <= rad_g[nt_g[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
#                    Av=0
                    if np.isnan(in_gas[nt_g[k]]):
                        spect_gf=spect_gf
                    else:             
                        if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:   
                            dust=10**(-0.4*Av*dust_rat_gas)
                            spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust*10.0**(-2.18+2.18-3.18)#+0.3+0.6)#0.01#*mass_g[nt_g[k]]
                            spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                            spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                            spect_g=spect_sfg+spect_g  
                ext_temp[2,j,i]=sfr_t
                spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
            photo_imag[:,j,i]=spect+spect_gf
             
    h1=pyf.PrimaryHDU(photo_imag)#.header
    #h2=pyf.ImageHDU(ifu_e)
    #h3=pyf.ImageHDU(ifu_1)
    #h4=pyf.ImageHDU(ifu_m)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dus_m=1
    if dus_m==1:
        h2=pyf.PrimaryHDU(ext_temp)
        hf=h2.header
        hf["NAXIS"]=3
        hf["NAXIS1"]=nl
        hf["NAXIS2"]=nl
        hf["NAXIS3"]=3
        hf["CRVAL1"]=0#oap
        hf["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
        hf["CRPIX1"]=nl/2
        hf["CTYPE1"]='RA---TAN'
        hf["CRVAL2"]=0#oap
        hf["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
        hf["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CRPIX2"]=nl/2
        hf["CTYPE2"]='DEC--TAN'
        hf['CUNIT1']='deg     '                                           
        hf['CUNIT2']='deg     '
        hf['RADECSYS']='ICRS    '
        hf['SYSTEM']='FK5     '
        hf['EQUINOX']=2000.00
       # hf['Type1']=('DUST    ','Av')
        hf['Type1']=('MASS    ','Msun')
        hf['Type2']=('VEL     ','km/s')
        hf['Type3']=('SFR     ','Msun/yr')
        hf['PSF']=sig
        hf['FOV']=fov
        hf['KPCSEC']=dkpcs
        hf['CAMX']=observer[0]#0
        hf['CAMY']=observer[1]#0
        hf['CAMZ']=observer[2]#cam
        hf['REDSHIFT']=float(red_0)
        hf['H0']=ho
        hf['Lambda_0']=Lam
        hf['Omega_m']=Om
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        out_fit_d=dir_o+outf+'_val.fits'
        wfits_ext(out_fit_d,hlist)
        dir_od1=dir_o.replace(" ","\ ")
        out_fitd1=dir_od1+outf+'_val.fits'
        sycall('gzip -f '+out_fitd1)
    #band_photo_r(photo_imag, h, name=outf+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir_o)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)

def sim_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="../../Base_bc03/templete_bc03_2.fits",template="../home/sanchez/ppak/legacy/gsd61_156.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    #radL=10.0*3.08567758e16*100.0
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600#+ran.randn(len(rad))*1.43/2.0
    phi=phi*180/np.pi*3600#+ran.randn(len(rad))*1.43/2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600#+ran.randn(len(rad_g))*1.43/2.0
    phi_g=phi_g*180/np.pi*3600#+ran.randn(len(rad_g))*1.43/2.0
    
    #template="/home/hjibarram/FIT3D_py/Base_bc03/Base_bc17/bc17_salp_Agelin_Metlin_330.fits"
    
    #template2="templete_gas.fits"
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    [ssp_template1,wave1,age_ssp1,met_ssp1,ml_ssp1,crval_w1,cdelt_w1,crpix_w1]=ssp_extract(template2)
    ml_ssp1=1./ml_ssp1
    #for i in range(0, len(age_ssp1)):
    #    print age_ssp1[i],"  ",met_ssp1[i],"  ",ml_ssp1[i]
    #sys.exit()
    n_ages=num_ages(age_ssp)
    ages_r=arg_ages(age_ssp)
    sim_imag=np.zeros([n_ages,nl,nl])
    sim_imag2=np.zeros([n_ages,nl,nl])
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            if len(nt) > 0:
                #if i == 220 & j == 220:
                mass_t=asosiate_ages(age_ssp,age_s[nt],mass_s[nt])
                ligh_t=asosiate_light(age_ssp,age_ssp1,met_ssp1,ml_ssp1,age_s[nt],mass_s[nt],met_s[nt])
                sim_imag[:,j,i]=mass_t   
                sim_imag2[:,j,i]=ligh_t  
    #print np.sum(sim_imag2)   
    #print np.sum(sim_imag)
    #print np.sum(sim_imag)/np.sum(sim_imag2)
    #sys.exit()
    h1=pyf.PrimaryHDU(sim_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mook Ilustris IFS Mass"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Msun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
    
    h2=pyf.PrimaryHDU(sim_imag2)
    h=h2.header
    h["NAXIS"]=3
    h["NAXIS3"]=n_ages 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
    h["COMMENT"]="Mook Ilustris IFS Light"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    h['UNITS']='Lsun'
    for kk in range(0, n_ages):
        h['AGE'+str(kk)]=ages_r[kk]
    hlist=pyf.HDUList([h2])
    hlist.update_extend()
    out_fit=dir_o+outf+'_L.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'_L.fits'
    sycall('gzip -f '+out_fit1)
    
    
def photosim_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    dlam=(1+(v_rad/vel_light+reds*1.0))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    cdelt_w=30.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)
    nw=len(wave_f)
    nw_s=len(wave)
    #import matplotlib.pyplot as plt
    photosim_imag=np.zeros([nw,nl,nl])
    phit=phi
    the1=the
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                for k in range(0, len(nt)):
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)#*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_t=spect_sf+spect_t
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            photosim_imag[:,j,i]=spect
             
    h1=pyf.PrimaryHDU(photosim_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
    
def photosimextgas_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir_o='',red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=200,fov=0.2,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
    sfri=sfri+1e-6
    mass_gssp=sfri*100e6
    Rs=vol#float_((vol/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    sup=4.0*np.pi*Rs**2.0#4.0*np.pi*(3.0*vol/4.0/np.pi)**(2.0/3.0)*(3.08567758e19*100)**2.0
    vel_light=299792.458
    no_nan=0
    fibA=150.
    fibB=120.
    leng_s=0.5#kpc
    scalep=2.0/fibB
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    cam=cd.comoving_distance(red_0, **cosmo)*1e3
    dap=fov/nl
    oap=-fov/2.0
    xima=np.zeros(nl)
    yima=np.zeros(nl)
    rad=np.sqrt(x**2.+y**2.+(cam-z)**2.)
    dkpcs=cam*(1./3600.)*(np.pi/180.)
    v_rad=(vx*x+vy*y+vz*(z-cam))/rad   
    rad_g=np.sqrt(x_g**2.+y_g**2.+(cam-z_g)**2.)
    v_rad_g=(vx_g*x_g+vy_g*y_g+vz_g*(z_g-cam))/rad_g 
    reds=reds_cos(rad/1e3)
    radA=rad/(1+reds)
    radL=np.array(rad*(1+reds)*(3.08567758e19*100))
    reds_g=reds_cos(rad_g/1e3)
    dlam=(1+(v_rad/vel_light+reds*1.0))
    dlam_g=(1+(v_rad_g/vel_light+reds_g))
    radA_g=rad_g/(1+reds_g)
    radL_g=np.array(rad_g*(1+reds_g)*(3.08567758e19*100))
    phi=np.arcsin(x/radA)
    the=np.arcsin(y/(radA*np.cos(phi)))
    the=the*180/np.pi*3600
    phi=phi*180/np.pi*3600
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600
    phi_g=phi_g*180/np.pi*3600
    Av=1.1
    #template="../home/sanchez/ppak/legacy/gsd61_156.fits"
    
    
    [ssp_template,wave,age_ssp,met_ssp,ml_ssp,crval_w,cdelt_w,crpix_w]=ssp_extract(template)
    ml_ssp=1.0/ml_ssp
    [gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,crval_g,cdelt_g,crpix_g]=gas_extract(template2)
    in_ssp=asosiate_ssp(ssp_template,wave,age_ssp,met_ssp,ml_ssp,age_s,met_s)
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
    in_gas=asosiate_gas(gas_template,wave_g,pht_gas,met_gas,den_gas,tem_gas,ha_gas,pht_g,met_g,nh,temp_g)
    dust_rat_ssp=A_l(3.1,wave)
    dust_rat_gas=A_l(3.1,wave_g)
    cdelt_w=30.0#5.0
    crval_w=500
    crpix_w=1.0
    wave_f=np.arange(crval_w,25000.0,cdelt_w)#25000
    s_nr=noise_sig(wave_f,25)
    band_g=np.ones(len(met_g))
    band_g[np.where((pht_g == 0))[0]]=1.0 # &(in_gas == -100)
    nw=len(wave_f)
    nw_s=len(wave)
    nw_g=len(wave_g)
    #import matplotlib.pyplot as plt
    #plt.plot(phi,the,'o')
    #plt.ylabel('some numbers')
    #plt.show()
    photo_imag=np.zeros([nw,nl,nl])
    ext_temp=np.zeros([10,nl,nl])
    phit=phi#+ra1
    the1=the#+de1
    phi_gt=phi_g#+ra2
    the_gt=the_g#+de2
    xo=-nl/2*dap
    yo=-nl/2*dap
    xi=xo
    xf=xo
    for i in range(0, nl):
        #print i,nl
        sycall('echo '+str(i)+'  '+str(nl))
        xi=xf
        xf=xf+dap
        yi=yo
        yf=yo
        for j in range(0, nl):
            yi=yf
            yf=yf+dap
            nt=np.where((phit < xf) & (phit >= xi) & (the1 < yf) & (the1 >= yi))[0]
            nt_g=np.where((phi_gt < xf) & (phi_gt >= xi) & (the_gt < yf) & (the_gt >= yi))[0]
            spect_t=np.zeros(nw_s)
            spect=np.zeros(nw)
            spect_g=np.zeros(nw_g)
            spect_gf=np.zeros(nw)
            if len(nt) > 0:
                Av_s=0
                mass_t=np.sum(mass_s[nt])
                vel_t=np.average(v_rad[nt])
                Lt=0
                Ft=0
                met_ligt=0
                met_mas=0
                age_ligt=0
                age_mas=0
                Av_ligt=0
                Av_flux=0
                for k in range(0, len(nt)):
                    nt_e=np.where(rad_g[nt_g] <= rad[nt[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
                    Av_s=Av+Av_s
                    if np.isnan(in_ssp[nt[k]]):
                        spect=spect
                    else:
                        if in_ssp[nt[k]] > 0 and in_ssp[nt[k]] < 1346:
                            dust=10**(-0.4*Av*dust_rat_ssp*0.44)
                            spect_s=ssp_template[in_ssp[nt[k]],:]/ml_ssp[in_ssp[nt[k]]]*mass_s[nt[k]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)*dust#/20.0
                            spect_sf=shifts(spect_s,wave,dlam[nt[k]])
                            spect_sf=spect_sf+ran.randn(nw_s)*np.median(spect_sf)*0.01
                            spect_t=spect_sf+spect_t
                            #ml_t=ml_t+ml_ssp[in_ssp[nt[k]]]#*20.0
                            Lt=Lt+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]#/20.0
                            Ft=Ft+mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)#/20.0
                            met_ligt=np.log10(met_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+met_ligt
                            met_mas=np.log10(met_s[nt[k]])*mass_s[nt[k]]+met_mas
                            age_ligt=np.log10(age_s[nt[k]])*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+age_ligt
                            age_mas=np.log10(age_s[nt[k]])*mass_s[nt[k]]+age_mas
                            Av_ligt=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]+Av_ligt
                            Av_flux=10**(-0.4*Av)*mass_s[nt[k]]/ml_ssp[in_ssp[nt[k]]]*3.846e33/(4.0*np.pi*radL[nt[k]]**2.0)/1e-16*10**(-0.4*Av*0.44)+Av_flux
                if Lt > 0:       
                    ext_temp[0,j,i]=Av_ligt/Lt
                    ext_temp[9,j,i]=Av_flux/Ft
                    ext_temp[5,j,i]=10.0**(met_ligt/Lt)
                    ext_temp[7,j,i]=10.0**(age_ligt/Lt)
                if mass_t > 0:
                    ext_temp[6,j,i]=10.0**(met_mas/mass_t)
                    ext_temp[8,j,i]=10.0**(age_mas/mass_t)
                ext_temp[1,j,i]=mass_t
                ext_temp[2,j,i]=vel_t
                ext_temp[4,j,i]=Lt
                spect=interp1d(wave,spect_t,bounds_error=False,fill_value=0.)(wave_f)
                spect[np.isnan(spect)]=0
            if len(nt_g) > 0:
                sfr_t=np.sum(sfri[nt_g])
                for k in range(0, len(nt_g)):
                    nt_e=np.where(rad_g[nt_g] <= rad_g[nt_g[k]])[0]
                    if len(nt_e) > 0:
                        Av=np.sum(Av_g[nt_g[nt_e]])
                    else:
                        Av=0
                    if np.isnan(in_gas[nt_g[k]]):
                        spect_gf=spect_gf
                    else:             
                        if in_gas[nt_g[k]] > 0 and in_gas[nt_g[k]] < 525:   
                            dust=10**(-0.4*Av*dust_rat_gas)
                            spect_sg=gas_template[in_gas[nt_g[k]],:]/ha_gas[in_gas[nt_g[k]]]*3.846e33*band_g[nt_g[k]]/(4.0*np.pi*radL_g[nt_g[k]]**2.0)*dust*10.0**(-2.18+2.18-1.18)#-3.18)#+0.3+0.6)#0.01#*mass_g[nt_g[k]]
                            spect_sfg=shifts(spect_sg,wave_g,dlam_g[nt_g[k]])
                            spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                            spect_g=spect_sfg+spect_g  
                ext_temp[3,j,i]=sfr_t
                spect_gf=interp1d(wave_g,spect_g,bounds_error=False,fill_value=0.)(wave_f)
                spect_gf[np.isnan(spect_gf)]=0
            photo_imag[:,j,i]=spect+spect_gf
    ext_temp[1,:,:]=np.log10(ext_temp[1,:,:]+1.0)
    ext_temp[4,:,:]=np.log10(ext_temp[4,:,:]+1.0)
    ext_temp[0,:,:]=-2.5*np.log10(ext_temp[0,:,:]+0.0001)
    ext_temp[9,:,:]=-2.5*np.log10(ext_temp[9,:,:]+0.0001)
    h1=pyf.PrimaryHDU(photo_imag)
    h=h1.header
    h["NAXIS"]=3
    h["NAXIS3"]=nw 
    h["NAXIS1"]=nl
    h["NAXIS2"]=nl
   # h["COMMENT"]="Mook Ilustris IFS"
    h["CRVAL1"]=0#oap
    h["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
    h["CRPIX1"]=nl/2
    h["CTYPE1"]='RA---TAN'
    h["CRVAL2"]=0#oap
    h["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
    h["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
    h["CRPIX2"]=nl/2
    h["CTYPE2"]='DEC--TAN'
    h['CUNIT1']='deg     '                                           
    h['CUNIT2']='deg     '
    h['CDELT3']=cdelt_w
    h['CRPIX3']=crpix_w
    h['CRVAL3']=crval_w
    h['CUNIT3']='Wavelength [A]'
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=fov
    h['KPCSEC']=dkpcs
    h['CAMX']=observer[0]#0
    h['CAMY']=observer[1]#0
    h['CAMZ']=observer[2]#cam
    h['REDSHIFT']=float(red_0)
    h['H0']=ho
    h['Lambda_0']=Lam
    h['Omega_m']=Om
    hlist=pyf.HDUList([h1])
    hlist.update_extend()
    out_fit=dir_o+outf+'.fits'
    wfits_ext(out_fit,hlist)
    dus_m=1
    if dus_m==1:
        h2=pyf.PrimaryHDU(ext_temp)
        hf=h2.header
        hf["NAXIS"]=3
        hf["NAXIS1"]=nl
        hf["NAXIS2"]=nl
        hf["NAXIS3"]=10
        hf["CRVAL1"]=0#oap
        hf["CD1_1"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CD1_2"]=np.sin(thet*np.pi/180.)*dap/3600.
        hf["CRPIX1"]=nl/2
        hf["CTYPE1"]='RA---TAN'
        hf["CRVAL2"]=0#oap
        hf["CD2_1"]=-np.sin(thet*np.pi/180.)*dap/3600.
        hf["CD2_2"]=np.cos(thet*np.pi/180.)*dap/3600.
        hf["CRPIX2"]=nl/2
        hf["CTYPE2"]='DEC--TAN'
        hf['CUNIT1']='deg     '                                           
        hf['CUNIT2']='deg     '
        hf['RADECSYS']='ICRS    '
        hf['SYSTEM']='FK5     '
        hf['EQUINOX']=2000.00
        #hf['Type1']=('DUST    ','Av')
        #hf['Type2']=('MASS    ','Msun')
        #hf['Type3']=('VEL     ','km/s')
        #hf['Type4']=('SFR     ','Msun/yr')
        hf['Type0']=('Av  ','Mag')
        hf['Type1']=('MASS    ','log10(Msun)')
        hf['Type2']=('VEL     ','km/s')
        hf['Type3']=('SFR     ','Msun/yr')
        hf['Type4']=('LUM     ','log10(Lsun)')
        hf['Type5']=('Z_lw    ','Z/H')
        hf['Type6']=('Z_mw    ','Z/H')
        hf['Type7']=('AGE_lw  ','Gyr')
        hf['Type8']=('AGE_mw  ','Gyr') 
        hf['Type8']=('Av_flx  ','Mag')    
        #hf['PSF']=sig
        hf['FOV']=fov
        hf['KPCSEC']=dkpcs
        hf['CAMX']=observer[0]#0
        hf['CAMY']=observer[1]#0
        hf['CAMZ']=observer[2]#cam
        hf['REDSHIFT']=float(red_0)
        hf['H0']=ho
        hf['Lambda_0']=Lam
        hf['Omega_m']=Om
        hlist=pyf.HDUList([h2])
        hlist.update_extend()
        out_fit_d=dir_o+outf+'_val.fits'
        wfits_ext(out_fit_d,hlist)
        dir_od1=dir_o.replace(" ","\ ")
        out_fitd1=dir_od1+outf+'_val.fits'
        sycall('gzip -f '+out_fitd1)
    dir_o1=dir_o.replace(" ","\ ")
    out_fit1=dir_o1+outf+'.fits'
    sycall('gzip -f '+out_fit1)
    
def noise_sig(wave,sn):
    w1=3900
    w2=10100
    dw2=100
    dw1=100
    s=sn/(1.+np.exp(-(wave-w1)/dw1))/(1.+np.exp((wave-w2)/dw2))+1.0
    return s

def inst_disp(pdl_lamb,pdl_flux,sigma_inst):
    vel_light=299792.458
    sigma=sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    rsigma=sigma/dpix
    if sigma > 0.1:
        box=int(3.0*rsigma*2.0)
        if box < 3:
            box=3
        kernel=np.zeros(2*box+1)
        norm=0
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
            kernel[j]=gaus
            norm=norm+gaus
        kernel=kernel/norm
        pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
    else:
        pdl_flux_conv_inst=pdl_flux
    return pdl_flux_conv_inst

def spec_resol(pdl_lamb,pdl_flux,R):
    vel_light=299792.458
    #dl=5500.0/R
    #sigma=dl#sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    pdl_flux_conv_inst=np.zeros(len(pdl_lamb))
    for i in range(0, len(pdl_lamb)):
        sigma=pdl_lamb[i]/R
        rsigma=sigma/dpix
        if sigma > 0.1:
            box=int(3.0*rsigma*2.0*10.0)
            if box < 3:
                box=3
            nk=2*box+1
            kernel=np.zeros(nk)
            norm=0
            for j in range(0, 2*box+1):
                gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
                kernel[j]=gaus
                norm=norm+gaus
            kernel=kernel/norm
            if i < box:
                kernel=kernel[box+i:nk]
            if len(pdl_lamb)-1-i < box+1:
                kernel=kernel[0:len(pdl_lamb)-i+box]
            nlk=len(kernel)
            if i < box:
                pdl_flux_t=pdl_flux[0:nlk]
                pdl_lamb_t=pdl_lamb[0:nlk]
            elif len(pdl_lamb)-1-i < box+1:
                pdl_flux_t=pdl_flux[i-box:len(pdl_lamb)]
                pdl_lamb_t=pdl_lamb[i-box:len(pdl_lamb)]
            else:
                pdl_flux_t=pdl_flux[i-box:i-box+nlk]
                pdl_lamb_t=pdl_lamb[i-box:i-box+nlk]
            #print pdl_lamb_t.shape
            #print kernel.shape
            #print nlk
            #print pdl_lamb_t[nlk-1]
            pdl_flux_t=pdl_flux_t*kernel
            val=simpson_r(pdl_flux_t,pdl_lamb_t,0,nlk-2)
            #pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
        else:
            val=pdl_flux[i]
            #pdl_flux_conv_inst=pdl_flux
        pdl_flux_conv_inst[i]=val
    return pdl_flux_conv_inst

def spec_resol_old(pdl_lamb,pdl_flux,R):
    vel_light=299792.458
    dl=5500.0/R
    sigma=dl#sigma_inst/vel_light*5000.0
    dpix=pdl_lamb[1]-pdl_lamb[2]
    rsigma=sigma/dpix
    if sigma > 0.1:
        box=int(3.0*rsigma*2.0)
        if box < 3:
            box=3
        kernel=np.zeros(2*box+1)
        norm=0
        for j in range(0, 2*box+1):
            gaus=np.exp(-0.5*(((j-box)/rsigma)**2.0))    
            kernel[j]=gaus
            norm=norm+gaus
        kernel=kernel/norm
        pdl_flux_conv_inst = convolve1d(pdl_flux,kernel)#,mode='same')
    else:
        pdl_flux_conv_inst=pdl_flux
    return pdl_flux_conv_inst
     
def shifts(spect_s,wave,dlam):
    wave_f=wave*dlam
    spect_f=interp1d(wave_f, spect_s,bounds_error=False,fill_value=0.)(wave)
    nt=np.where(spect_f == 0)[0]
    #spect_f[nt]=np.nan
    return  spect_f 

def asosiate_gas(ssp_temp,wave,age_s,met_s,den_s,tem_s,ml_s,age,met,den,tem):
    met_s=float_(met_s)*0.02#127
    age_a=[]
    met_a=[]
    den_a=[]
    tem_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    ban2=0
    ban3=0
    age_t=sorted(age_s)
    met_t=sorted(met_s)
    den_t=sorted(den_s)
    tem_t=sorted(tem_s)
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    den_a.extend([den_t[0]])
    tem_a.extend([tem_t[0]])
    ind_ssp=np.zeros((len(met)),dtype=np.int)
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
        if den_t[i-1] > den_t[i]:
                ban2 =1
        if den_t[i-1] < den_t[i] and ban2 == 0:
            den_a.extend([den_t[i]])
        if tem_t[i-1] > tem_t[i]:
                ban3 =1
        if tem_t[i-1] < tem_t[i] and ban3 == 0:
            tem_a.extend([tem_t[i]])
    n_age=len(age_a)
    n_met=len(met_a)
    n_den=len(den_a)
    n_tem=len(tem_a)
    ind_age=ages_definition_l(age,age_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    ind_den=val_definition_l(den[ind_age[i][ind_met[j]]],den_a) 
                    for k in range(0, n_den):
                        if len(ind_den[k]) > 0:
                            ind_tem=val_definition_l(tem[ind_age[i][ind_met[j][ind_den[k]]]],tem_a)
                            for h in range(0, n_tem):
                                if len(ind_tem[h]) >  0:               
                                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]) & (den_s == den_a[k]) & (tem_s == tem_a[h]))[0]
                                    ind_ssp[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]=nt[0]
                                    #print nt[0]
                                    #print "QT=",age_s[nt[0]],";",age_a[i+1],age_a[i-1]
                                    #print age[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "Z=",met_s[nt[0]],";",met_a[j+1],met_a[j-1]
                                    #print met[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "nh=",den_s[nt[0]],";",den_a[k+1],den_a[k-1]
                                    #print den[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #print "Tem=",tem_s[nt[0]],";",tem_a[h+1],tem_a[h-1]
                                    #print tem[ind_age[i][ind_met[j][ind_den[k][ind_tem[h]]]]]
                                    #sys.exit()
    return ind_ssp

def asosiate_pho2(ssp_temp,wave,age_s,met_s,ml,massp_s,metp_s,agep_s,x_g,y_g,z_g,x_s,y_s,z_s,R_g):  
    mass=np.zeros(len(z_g))
    met=np.zeros(len(z_g))
    age=np.zeros(len(z_g))  
    for i in range(0, len(z_g)):
        r=np.sqrt((x_s-x_g[i])**2.0+(y_s-y_g[i])**2.0+(z_s-z_g[i])**2.0)
        nt=np.where(r <= R_g[i])[0]
        met[i]=np.nanmean(metp_s[nt])
        age[i]=np.nanmean(agep_s[nt])
        mass[i]=np.nansum(massp_s[nt])    
    vel_light=299792458.0
    h_p=6.62607004e-34
    #age=np.ones(len(met))*2.5e6/1e9
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    photo=np.zeros(len(age))
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
                    flux_0=ssp_temp[nt[0],:]/ml[nt[0]]/(h_p*vel_light/wave/1e-10/1e-7)*3.846e33
                    j1=0#int(0.47*n_c)
                    j2=int(0.63*len(wave))
                    norm=simpson_r(flux_0,wave,j1,j2)
                    photo[ind_age[i][ind_met[j]]]=norm*mass[ind_age[i][ind_met[j]]]+1#/(4.0*np.pi*Rs[ind_age[i][ind_met[j]]]**2.0*n_h[ind_age[i][ind_met[j]]])+1
    return np.log10(photo)


def asosiate_pho(ssp_temp,wave,age_s,met_s,ml,mass,met,Rs,n_h):
    vel_light=299792458.0
    h_p=6.62607004e-34
    age=np.ones(len(met))*2.5e6/1e9
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    photo=np.zeros(len(age))
    ind_ssp[:]=-100
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
                ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
                    flux_0=ssp_temp[nt[0],:]/ml[nt[0]]/(h_p*vel_light/wave/1e-10/1e-7)*3.846e33
                    j1=0#int(0.47*n_c)
                    j2=int(0.63*len(wave))
                    norm=simpson_r(flux_0,wave,j1,j2)
                    photo[ind_age[i][ind_met[j]]]=norm*mass[ind_age[i][ind_met[j]]]+1#/(4.0*np.pi*Rs[ind_age[i][ind_met[j]]]**2.0*n_h[ind_age[i][ind_met[j]]])+1
    return np.log10(photo)

def num_ages(age_s):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    n_age=len(age_a)
    return n_age

def reverse_numeric(x, y):
    return y - x

def arg_ages(age_s):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    return age_a[::-1]
  
def asosiate_ages(age_s,age,mass):
    age_a=[]
    nssp=len(age_s)
    ban=0
    age_t=sorted(age_s)
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    mass_f=np.zeros(n_age)
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            mass_f[i]=np.sum(mass[ind_age[i]])
    return mass_f[::-1]

def asosiate_light(age_s,age1_s,met1_s,ml1_s,age,mass,met):
    age_a=[]
    age1_a=[]
    met1_a=[]
    nssp=len(age_s)
    nssp1=len(age1_s)
    ban=0
    ban0=0
    ban1=0
    age_t=sorted(age_s)
    age1_t=sorted(age1_s)
    met1_t=met1_s#sorted(met_s)
    #ml_t=ml_t[np.argsort(age_s)]
    age_a.extend([age_t[0]])
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
    age1_a.extend([age1_t[0]])
    met1_a.extend([met1_t[0]])
    for i in range(1, nssp1):
        if age1_t[i-1] > age1_t[i]:
            ban0 =1
        if age1_t[i-1] < age1_t[i] and ban0 == 0:
            age1_a.extend([age1_t[i]])
        if met1_t[i-1] > met1_t[i]:
            ban1 =1
        if met1_t[i-1] < met1_t[i] and ban1 == 0:
            met1_a.extend([met1_t[i]])
    ind_age=ages_definition_l(age,age_a)
    ind_age1=ages_definition_l(age1_s,age_a)
    n_age=len(age_a)
    n_age1=len(age1_a)
    n_met1=len(met1_a)
    ligh_f=np.zeros(n_age)
    for i in range(0, n_age):
#        print "Age=",age_a[i]
        if len(ind_age[i]) > 0:
            temp=0.0
            ind_age2=ages_definition_l(age[ind_age[i]],age1_a)
            for j in range(0, n_age1):
                if len(ind_age2[j]) > 0:   
                    ind_met=met_definition_l(met[ind_age[i][ind_age2[j]]],met1_a)
                    for k in range(0, n_met1):
                        if len(ind_met[k]) > 0:
                            nt=np.where((age1_s == age1_a[j]) & (met1_s == met1_a[k]))[0]
                            #temp=np.sum(mass[ind_age[i][ind_age2[j][ind_met[k]]]])+temp
                            temp=np.sum(mass[ind_age[i][ind_age2[j][ind_met[k]]]]/ml1_s[nt[0]])+temp
                            #print "MLr=",ml1_s[nt[0]],age1_s[nt[0]],met1_s[nt[0]]
            ligh_f[i]=temp
 #   print "_________________________________________"
    return ligh_f[::-1]
        
def asosiate_ssp(ssp_temp,wave,age_s,met_s,ml_s,age,met):
    age_a=[]
    met_a=[]
    nssp=len(met_s)
    ban=0
    ban1=0
    age_t=sorted(age_s)
    met_t=met_s
    age_a.extend([age_t[0]])
    met_a.extend([met_t[0]])
    ind_ssp=np.zeros((len(age)),dtype=np.int)
    ind_ssp[:]=np.nan
    for i in range(1, nssp):
        if age_t[i-1] > age_t[i]:
            ban =1
        if age_t[i-1] < age_t[i] and ban == 0:
            age_a.extend([age_t[i]])
        if met_t[i-1] > met_t[i]:
            ban1 =1
        if met_t[i-1] < met_t[i] and ban1 == 0:
            met_a.extend([met_t[i]])
    ind_age=ages_definition_l(age,age_a)
    n_age=len(age_a)
    n_met=len(met_a)
    print n_age
    print n_met
    for i in range(0, n_age):
        if len(ind_age[i]) > 0:
            ind_met=met_definition_l(met[ind_age[i]],met_a)
            for j in range(0, n_met):
                if len(ind_met[j]) > 0:
                    nt=np.where((age_s == age_a[i]) & (met_s == met_a[j]))[0]
                    ind_ssp[ind_age[i][ind_met[j]]]=nt[0]
    return ind_ssp

def gas_extract(template):
    [pdl_flux_c_ini,hdr]=gdata(template, 0, header=True)
    [nf,n_c]=pdl_flux_c_ini.shape
    coeffs=np.zeros([nf,3])
    crpix=hdr["CRPIX1"]
    cdelt=hdr["CDELT1"]
    crval=hdr["CRVAL1"]
    tem_mod=[]
    pht_mod=[]
    den_mod=[]
    met_mod=[]
    Ha=[]
    name=[]
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_gas_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        TEM=data[3]
        PHT=data[2]
        DEN=data[1]
        MET=data[0]
        tem=float_(TEM.replace('t',''))
        pht=float_(PHT.replace('q',''))
        den=float_(DEN.replace('n',''))
        met=float_(MET.replace('z',''))
        tem_mod.extend([tem])    
        pht_mod.extend([pht])
        den_mod.extend([den])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr[header])
        Ha.extend([val_ml])
    wave_c=[]
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c.extend([(crval+cdelt*(j+1-crpix))])
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])
    wave_c=np.array(wave_c)
    return [pdl_flux_c_ini,wave_c,pht_mod,met_mod,den_mod,tem_mod,Ha,crval,cdelt,crpix]
    
def ssp_extract(template):
    [pdl_flux_c_ini,hdr]=gdata(template, 0, header=True)
    [nf,n_c]=pdl_flux_c_ini.shape
    coeffs=np.zeros([nf,3])
    crpix=hdr["CRPIX1"]
    cdelt=hdr["CDELT1"]
    crval=hdr["CRVAL1"]
    age_mod=[]
    met_mod=[]
    ml=[]
    name=[]
    for iii in range(0, nf):
        header="NAME"+str(iii)
        name.extend([hdr[header]]);
        name_min=name[iii]
        name_min=name_min.replace('spec_ssp_','')
        name_min=name_min.replace('.spec','')    
        name_min=name_min.replace('.dat','')
        data=name_min.split('_')
        AGE=data[0]
        MET=data[1]
        if 'Myr' in AGE:
            age=AGE.replace('Myr','')
            age=float_(age)/1000.
        else:
            age=AGE.replace('Gyr','')
            age=float_(age)
        met=float_(MET.replace('z','0.'))
        age_mod.extend([age])
        met_mod.extend([met])
        header="NORM"+str(iii)    
        val_ml=float_(hdr[header])
        if val_ml != 0:
            ml.extend([1/val_ml])
        else:
            ml.extend([1])
    wave_c=[]
    dpix_c_val=[]
    for j in range(0, n_c):
        wave_c.extend([(crval+cdelt*(j+1-crpix))])
        if j > 0:
            dpix_c_val.extend([wave_c[j]-wave_c[j-1]])
    wave_c=np.array(wave_c)
    ml=np.array(ml)
    age_mod=np.array(age_mod)
    met_mod=np.array(met_mod)

    return [pdl_flux_c_ini,wave_c,age_mod,met_mod,ml,crval,cdelt,crpix]

def simpson_r(f,x,i1,i2,typ=0):
    n=(i2-i1)*1.0
    if n % 2:
        n=n+1.0
        i2=i2+1
    b=x[i2]
    a=x[i1]
    h=(b-a)/n
    s= f[i1]+f[i2]
    n=int(n)
    dx=b-a
    for i in range(1, n, 2):
        s += 4 * f[i1+i]
    for i in range(2, n-1, 2):
        s += 2 * f[i1+i]
    if typ == 0:
        return s*h/3.0
    if typ == 1:
        return s*h/3.0/dx
             
def read_data(dir,file1,file2):
#    file1='all_stars_a1.001.out'
#    file1='all_stars_z0.0_r500kpc_centered.out'
    #file2='gas_component_5Rvir_a1.001.out'
    f1=open(dir+file1,'r')
    f2=open(dir+file2,'r')
    vel_str_x=[]
    vel_str_y=[]
    vel_str_z=[]
    pos_str_x=[]
    pos_str_y=[]
    pos_str_z=[]
    massf=[]
    mass0=[]
    age=[]
    metal1=[]
    metal2=[]
    scalf=[]
    for line in f1:
#    for j in range(0, 100000):
#        line=f1.readline()
        if not "#" in line: 
            line=line.replace('\n','')
            data=line.split(',')
            data=filter(None,data)
            pos_str_x.extend([float_(data[5])])
            pos_str_y.extend([float_(data[6])])
            pos_str_z.extend([float_(data[7])])
            vel_str_x.extend([float_(data[8])])
            vel_str_y.extend([float_(data[9])])
            vel_str_z.extend([float_(data[10])])
            massf.extend([float_(data[4])])
            mass0.extend([float_(data[3])])
            metal1.extend([float_(data[12])])
            metal2.extend([float_(data[13])])        
            age.extend([float_(data[1])])
            scalf.extend([float_(data[2])])
    f1.close()
    pos_str_x=np.array(pos_str_x)*1e3
    pos_str_y=np.array(pos_str_y)*1e3
    pos_str_z=np.array(pos_str_z)*1e3
    vel_str_x=np.array(vel_str_x)
    vel_str_y=np.array(vel_str_y)
    vel_str_z=np.array(vel_str_z)
    massf=np.array(massf)
    mass0=np.array(mass0)
    metal1=np.array(metal1)
    metal2=np.array(metal2)    
    #metal=1./((10.**metal*(1.0/0.0127-1.0))+1.)
    metal=(metal1+metal2)
    age=np.array(age)
    scalf=1.0/np.array(scalf)-1
    pos_str=np.zeros([len(pos_str_x),3])
    pos_str[:,0]=pos_str_x
    pos_str[:,1]=pos_str_y
    pos_str[:,2]=pos_str_z
    vel_str=np.zeros([len(pos_str_x),3])
    vel_str[:,0]=vel_str_x
    vel_str[:,1]=vel_str_y
    vel_str[:,2]=vel_str_z
    vel_gas_x=[]
    vel_gas_y=[]
    vel_gas_z=[]
    pos_gas_x=[]
    pos_gas_y=[]
    pos_gas_z=[]
    masses_g=[]
    dens_g=[]
    metal_g=[]
    sfr_g=[]
    Av_g=[]
    T_g=[]
    R_g=[]
    fac=[]
    for line in f2:
#    for j in range(0, 100000):
#        line=f2.readline()
        if not "*" in line:
            line=line.replace('\n','')
            data=line.split(',')
            data=filter(None,data)
            pos_gas_x.extend([float_(data[0])])
            pos_gas_y.extend([float_(data[1])])
            pos_gas_z.extend([float_(data[2])])
            vel_gas_x.extend([float_(data[3])])
            vel_gas_y.extend([float_(data[4])])
            vel_gas_z.extend([float_(data[5])])
            masses_g.extend([float_(data[9])])
            metal_g.extend([float_(data[10])])        
            dens_g.extend([float_(data[7])])
            sfr_g.extend([float_(data[11])])
            Av_g.extend([float_(data[6])])
            T_g.extend([float_(data[8])])
            R_g.extend([float_(data[12])])
            if float_(data[10]) > (10.0**(-0.59)*0.017):
                fac.extend([10.0**(2.21-1.0)])
            else:
                fac.extend([(10.0**(2.21-1.0))/(float_(data[10])/0.017)**(3.1-1.0)])
    f2.close()
    fac=np.array(fac)
    pos_gas_x=np.array(pos_gas_x)*1e3
    pos_gas_y=np.array(pos_gas_y)*1e3
    pos_gas_z=np.array(pos_gas_z)*1e3
    vel_gas_x=np.array(vel_gas_x)
    vel_gas_y=np.array(vel_gas_y)
    vel_gas_z=np.array(vel_gas_z)
    masses_g=np.array(masses_g)
    metal_g=np.array(metal_g)
    dens_g=np.array(dens_g)
    sfr_g=np.array(sfr_g)
    Av_g=np.array(Av_g)/0.017/fac/0.017#*0.5/0.017
    R_g=np.array(R_g)
    T_g=np.array(T_g)
    pos_gas=np.zeros([len(pos_gas_x),3])
    pos_gas[:,0]=pos_gas_x
    pos_gas[:,1]=pos_gas_y
    pos_gas[:,2]=pos_gas_z
    vel_gas=np.zeros([len(pos_gas_x),3])
    vel_gas[:,0]=vel_gas_x
    vel_gas[:,1]=vel_gas_y
    vel_gas[:,2]=vel_gas_z
    #f2.close()
    gas={'Coordinates':pos_gas,'Velocities':vel_gas,'Density':dens_g,'GFM_Metallicity':metal_g,'Radius':R_g,'StarFormationRate':sfr_g,'Extintion':Av_g,'Temperature':T_g,'GFM_Masses':masses_g}
    stars={'Coordinates':pos_str,'Velocities':vel_str,'Masses':massf,'GFM_InitialMass':mass0,'GFM_StellarFormationTime':scalf,'GFM_StellarAge':age,'GFM_Metallicity':metal}
    return [gas, stars]

def mock_halo_test(idh,basename='artsp8-',dir1='',fib_n=7,psf=0,nl=110,fov=30.0,thet=0.0,ifutype="MaNGA"):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
    cubef=basename+id+'.cube_test'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    dir1n=dir1.replace(" ","\ ")
    cube_conv_t(cubef,psfi=psf,dir_o=dir1,nl=fib_n,fov=fov,thet=thet,ifutype=ifutype)
    sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
                 
def mock_halo(idh,sp_res=0.0,sp_samp=1.25,basename='artsp8-',template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",file=file,file2='',dir1='',fib_n=7,psf=0,ho=0.704,Lam=0.7274,SN=15.0,Fluxm=20.0,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="MaNGA"):
    #if not "CALIFA" in ifutype or not "MaNGA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
#    cubef='cube_'+id
    cubef=basename+id+'.cube'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    [gas,stars]=read_data(dir0,file,file2)
    dx_g = gas['Coordinates'][:,0]
    dy_g = gas['Coordinates'][:,1]
    dz_g = gas['Coordinates'][:,2]
    vx_g = gas['Velocities'][:,0]
    vy_g = gas['Velocities'][:,1]
    vz_g = gas['Velocities'][:,2]
    ra_g = gas['Radius'][:]
    dens = gas['Density'][:]
    sfri = gas['StarFormationRate'][:]
    meta_g=gas['GFM_Metallicity'][:]
    mass_g=gas['GFM_Masses'][:]
    temp_g=gas['Temperature'][:]
    Av_g = gas['Extintion'][:]        
    dx = stars['Coordinates'][:,0]
    dy = stars['Coordinates'][:,1]
    dz = stars['Coordinates'][:,2]
    mass=stars['Masses'][:]
    mas0=stars['GFM_InitialMass'][:]
    meta=stars['GFM_Metallicity'][:]
    time=stars['GFM_StellarFormationTime'][:]
    vx = stars['Velocities'][:,0]
    vy = stars['Velocities'][:,1]
    vz = stars['Velocities'][:,2]
    ages=stars['GFM_StellarAge'][:]
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    mass_g=mass_g#/ho
    mass=mass[nt]#/ho
    mas0=mas0[nt]#/ho
    meta=meta[nt]#/0.0127
    ages=ages[nt]
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=ra_g/ho#**3.0
    dens=dens*ho**2.0
    zf=time[nt]
    xo=0#modes(dx,nit=7,n_s=0.8)
    yo=0#modes(dy,nit=7,n_s=0.8)
    zo=0#modes(dz,nit=7,n_s=0.8)
#    print xo,yo,zo,"HOLA",len(dx)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
#    print modes(x,nit=7,n_s=0.8),modes(y,nit=7,n_s=0.8),modes(z,nit=7,n_s=0.8)
#    sys.exit()
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    dir1n=dir1.replace(" ","\ ")
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        cube_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,sp_res=sp_res,sp_samp=sp_samp,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=fib_n,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
        #sys.exit()
        band_cube(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
    else:
        band_cube(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
    sycall('mv *'+basename+id+'* '+dir1n)
        

def mock_halo_s(idh,basename='artsp8-',template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",file=file,file2='',dir1='',psf=0,ho=0.704,Lam=0.7274,SN=15.0,Fluxm=20.0,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="SDSS"):
    #if not "CALIFA" in ifutype or not "MaNGA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
#    cubef='cube_'+id
    cubef=basename+id+'.spec'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    sycall('mkdir -p '+dir1t)
    [gas,stars]=read_data(dir0,file,file2)
    dx_g = gas['Coordinates'][:,0]
    dy_g = gas['Coordinates'][:,1]
    dz_g = gas['Coordinates'][:,2]
    vx_g = gas['Velocities'][:,0]
    vy_g = gas['Velocities'][:,1]
    vz_g = gas['Velocities'][:,2]
    ra_g = gas['Radius'][:]
    dens = gas['Density'][:]
    sfri = gas['StarFormationRate'][:]
    meta_g=gas['GFM_Metallicity'][:]
    mass_g=gas['GFM_Masses'][:]
    temp_g=gas['Temperature'][:]
    Av_g = gas['Extintion'][:]        
    dx = stars['Coordinates'][:,0]
    dy = stars['Coordinates'][:,1]
    dz = stars['Coordinates'][:,2]
    mass=stars['Masses'][:]
    mas0=stars['GFM_InitialMass'][:]
    meta=stars['GFM_Metallicity'][:]
    time=stars['GFM_StellarFormationTime'][:]
    vx = stars['Velocities'][:,0]
    vy = stars['Velocities'][:,1]
    vz = stars['Velocities'][:,2]
    ages=stars['GFM_StellarAge'][:]
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    mass_g=mass_g#/ho
    mass=mass[nt]#/ho
    mas0=mas0[nt]#/ho
    meta=meta[nt]#/0.0127
    ages=ages[nt]
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=ra_g/ho#**3.0
    dens=dens*ho**2.0
    zf=time[nt]
    xo=0#modes(dx,nit=7,n_s=0.8)
    yo=0#modes(dy,nit=7,n_s=0.8)
    zo=0#modes(dz,nit=7,n_s=0.8)
#    print xo,yo,zo,"HOLA",len(dx)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    dir1n=dir1.replace(" ","\ ")
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        fib_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val_mass_t.fits.gz '+dirf)

def mock_photo(id,basename='artsp8-',template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photo'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
        #[age_F,ind]=ages_definition(age_s,n_ages=55)
        photo_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
     #   band_photo(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
    #    [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".r_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
    #    sycall("mv "+cubef+".* "+dir1n)
        if ptt.exists(dir1+cubef+".rg_Lc_rad.fits") == True:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=0)
        else:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=[-100,-100,-100,-100,-100,-100]
        hdulist = fits.open(dir1+cubef+".rg_Lc.fits")
        hd=hdulist[0].header
        red_0=hd["REDSHIFT"]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    #sys.exit()
                
    if np.abs(not_fit) > 0:
        R50=-100
    #f4.write(cubef+','+str(R50)+','+str(e_rad)+','+str(red_0)+','+str(ab)+','+str(PA)+','+str(not_fit)+'\n')
    #return [PA,np.sqrt(1.0-ab**2.),e_rad,not_fit]

        
def mock_sim(id,basename='artsp8-',template2="../../Base_bc03/templete_bc03_2.fits",template="../home/sanchez/ppak/legacy/gsd61_156.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,1.5],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='simu-'+id
    cubef=basename+id+'.simu'
    dirf=dir1t+'photo_lin/'
    #dirf2=dir1t+'Ensamble/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
       # nt=np.where(np.abs(z) <= 2.0)[0]
       # x=x[nt]
        #y=y[nt]
        #z=z[nt]
        #vx=vx[nt]
        #vy=vy[nt]
        #vz=vz[nt]
        #ages=ages[nt]
        #meta=meta[nt]
        #mass=mass[nt]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3#cd.lookback_time(zf, **cosmo)/31557600./1e9
        #[age_F,ind]=ages_definition(age_s,n_ages=55)
        sim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
        #    ensamble_s(cubef+'.fits.gz',dir1=dir1)
        #    ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            #sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            #sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        #ensamble_int(cubef+'.fits.gz',dir1=dir1)
        #ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        #sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
        #    ensamble_s(cubef+'.fits.gz',dir1=dir1)
        #    ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            #sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            #sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        #ensamble_int(cubef+'.fits.gz',dir1=dir1)
        #ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        #sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        #sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
        
def mock_photosim(id,basename='artsp8-',template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photosim'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3
        photosim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #band_photo(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
        
        
def mock_photosimextgas(id,basename='artsp8-',template2="templete_gas.fits",template="/home/hjibarram/FIT3D_py/Base_bc03/templete_bc03_2.fits",dir1='',file=file,file2='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0]):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    outf='image_'+id
    #cubef='photo-'+id
    cubef=basename+id+'.photosimextgas'
    dirf=dir1t+'photo_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [gas,stars]=read_data(dir0,file,file2)
        dx_g = gas['Coordinates'][:,0]
        dy_g = gas['Coordinates'][:,1]
        dz_g = gas['Coordinates'][:,2]
        vx_g = gas['Velocities'][:,0]
        vy_g = gas['Velocities'][:,1]
        vz_g = gas['Velocities'][:,2]
        ra_g = gas['Radius'][:]
        dens = gas['Density'][:]
        sfri = gas['StarFormationRate'][:]
        meta_g=gas['GFM_Metallicity'][:]
        mass_g=gas['GFM_Masses'][:]
        temp_g=gas['Temperature'][:]
        Av_g = gas['Extintion'][:]           
        dx = stars['Coordinates'][:,0]
        dy = stars['Coordinates'][:,1]
        dz = stars['Coordinates'][:,2]
        mass=stars['Masses'][:]
        mas0=stars['GFM_InitialMass'][:]
        meta=stars['GFM_Metallicity'][:]
        time=stars['GFM_StellarFormationTime'][:]
        vx = stars['Velocities'][:,0]
        vy = stars['Velocities'][:,1]
        vz = stars['Velocities'][:,2]
        ages=stars['GFM_StellarAge'][:]
        nt=np.where(time > 0)[0]
        dx=dx[nt]/ho
        dy=dy[nt]/ho
        dz=dz[nt]/ho
        vx=vx[nt]
        vy=vy[nt]
        vz=vz[nt]
        mass_g=mass_g#/ho
        mass=mass[nt]#/ho
        mas0=mas0[nt]#/ho
        meta=meta[nt]#/0.0127
        ages=ages[nt]
        dx_g=dx_g/ho
        dy_g=dy_g/ho
        dz_g=dz_g/ho
        volm=ra_g/ho#**3.0
        dens=dens*ho**2.0
        zf=time[nt]
        xo=modes(dx,nit=7,n_s=0.8)
        yo=modes(dy,nit=7,n_s=0.8)
        zo=modes(dz,nit=7,n_s=0.8)
        x=dx-xo
        y=dy-yo
        z=dz-zo
        x_g=dx_g-xo
        y_g=dy_g-yo
        z_g=dz_g-zo
        x0=observer[0]-xo
        y0=observer[1]-yo
        z0=observer[2]-zo
        Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
        red_0=reds_cos(Rc/1e3)
        Ra=Rc#/(1+red_0)
        A1=np.arctan2(y0,z0)
        A2=np.arcsin(x0/Ra)
        R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
        R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
        R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
        Ve=np.array([x,y,z])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x=Vf[0]
        y=Vf[1]
        z=Vf[2]
        Ve=np.array([vx,vy,vz])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx=Vf[0]
        vy=Vf[1]
        vz=Vf[2]
        Ve=np.array([x_g,y_g,z_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        x_g=Vf[0]
        y_g=Vf[1]
        z_g=Vf[2]
        Ve=np.array([vx_g,vy_g,vz_g])
        Vf=np.dot(np.dot(R2,R1),Ve)
        vx_g=Vf[0]
        vy_g=Vf[1]
        vz_g=Vf[2]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        age_s=ages/1e3
        photosimextgas_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template2=template2,template=template,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)        
        if ptt.exists(dir1+cubef+".rg_Lc_rad.fits") == True:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=0)
        else:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
            
def ensamble_s(name, dir='', dir1='',rx=[0.0,0.5,1.0,1.5],lig=0):
    file2=dir1+name.replace('simu','photo').replace('_L','').replace(".fits.gz",".rg_Lc_rad.fits")
#    if ptt.exists(file2) == True:
    [pdl_rad, hdr2]=gdata(file2,0, header=True)
    ind=[]
    inda=[]
    nr=len(rx)
    rad=1.0
    for ii in range(0, nr-1):
        nt=np.where(pdl_rad< rx[ii+1]*rad)
        nta=np.where(pdl_rad[nt]>= rx[ii]*rad)
        ind.extend([nt])
        inda.extend([nta])
    vel_light=299792.458
    [mass_age_t,hdr]=gdata(dir1+name, 0, header=True)
    [n_age,nx,ny]=mass_age_t.shape
    ages=np.zeros(n_age)
    for i in range(0, n_age):
        ages[i]=np.log10(hdr["AGE"+str(i)])+9.0
    mass_age_t_2=mass_age_t
    mass_age_t_e=np.zeros([n_age,nx,ny])
    MassT=0
    massN=np.zeros(nr)
    massN_e=np.zeros(nr)
    mass=np.zeros([n_age,nr])
    mass_e=np.zeros([n_age,nr])
    sfrt=np.zeros([n_age,nr])
    sfdt=np.zeros([n_age])
    mass_temp_total=0
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
        if i == 0:
            temp1=temp
        else:
            temp1=temp1+temp
        MassT=MassT+np.sum(temp)
        for ii in range(0, nr-1):
            Dt_mass=np.sum(temp[ind[ii]][inda[ii]])
            Dt_mass_e=np.sum(temp_e[ind[ii]][inda[ii]])
            Dt_mass_2=np.sum(temp_2[ind[ii]][inda[ii]])
            massN[ii]=Dt_mass+massN[ii]
            massN_e[ii]=Dt_mass_e+massN_e[ii]
            mass[i,ii]=np.log10(massN[ii])
            mass_e[i,ii]=massN_e[ii]
            sfrt[i,ii]=Dt_mass_2/Dt_age
#    print mass[n_age-1,:]
        #print massN
        #print nr
#    sys.exit()
    mass_temp_total=np.log10(np.sum(10**mass[n_age-1,:]))
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    mass_n=10**(10**(mass-mass[n_age-1,:]))
    mass_n_e=np.sqrt((10**(mass-mass[n_age-1,:]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[n_age-1,:]/10**(2.0*mass[n_age-1,:]))))
    if lig == 1:
        f2=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_L.csv","w")
        f2.write("#  LOG_AGE  N_LIGHR_1  N_LIGHR_2  N_LIGHR_3  N_LIGHR_4  LOG_LIGHR_1  LOG_LIGHR_2  LOG_LIGHR_3  LOG_LIGHR_4 \n")
    else:
        f2=open(dir1+name.replace(".fits.gz","")+"_Ensemble.csv","w")
        f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    for i in range(0, n_age):
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
    f2.close()
    
    
    nrb=1+4*3+1
    ensamb=np.zeros([n_age,nrb])
    co=0
    if lig == 1:
        f=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_L.csv","r")
    else:
        f=open(dir1+name.replace(".fits.gz","")+"_Ensemble.csv","r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1

    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1:4])
    sfh=ensamb[:,7:10]#/10.0**ensamb[:,4:7]
    emgh=(ensamb[:,10:13])#*10
    #print ensamb[:,7:10]
    #print 10.0**ensamb[:,4:7]
    dth=ensamb[:,13]/1e3
    Ddth=np.zeros(len(dth))
    [nx,ny]=mgh.shape
    Dmgh=np.zeros([nx,ny])
    exy1=np.zeros([nx,ny])
    exy2=np.zeros([nx,ny])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i,:]-mgh[i+1,:])
            Ddth[i]=np.abs(dth[i]-dth[i+1])#/2.0
        elif i == len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i-1,:]-mgh[i,:])
            Ddth[i]=np.abs(dth[i-1]-dth[i])#/2.0 
    exy1[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)#np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]
    exy1[:,1]=np.amin([np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1],np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]],axis=0)#np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]
    exy1[:,2]=np.amin([np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2],np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]],axis=0)#np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]
    exy2[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)#np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0]
    exy2[:,1]=np.amin([np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1],np.abs(Dmgh[:,1]/(dth+Ddth)*dth)+emgh[:,1]],axis=0)#np.abs(Dmgh[:,1]/(dth-Ddth)*dth)+emgh[:,1]
    exy2[:,2]=np.amin([np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2],np.abs(Dmgh[:,2]/(dth+Ddth)*dth)+emgh[:,2]],axis=0)#np.abs(Dmgh[:,2]/(dth-Ddth)*dth)+emgh[:,2]

    for i in range(0, 3):
        val1=np.abs(Dmgh[:,i]/(dth+Ddth)*dth)+emgh[:,i]
        val2=np.abs(Dmgh[:,i]/(dth-Ddth)*dth)+emgh[:,i]
        for j in range(0, len(val1)):
            if val1[j] >= 0.15:
                val1[j]=0.1
            if val2[j] >= 0.15:
                val2[j]=0.1
        exy1[:,i]=val1
        exy2[:,i]=val2
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title("SIMULATION LGH INPUT")
    else:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title("SIMULATION MGH INPUT")
    plt.semilogx(age,mgh[:,0],'',color="b",label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=1.5)
    plt.semilogx(age,mgh[:,1],'--',color="g",label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=1.5)
    plt.semilogx(age,mgh[:,2],':',color="r",label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=1.5)
    ax.fill_between(age,mgh[:,0]+exy1[:,0],mgh[:,0]-exy2[:,0],alpha='0.20',color="b")
    ax.fill_between(age,mgh[:,1]+exy1[:,1],mgh[:,1]-exy2[:,1],alpha='0.20',color="g")
    ax.fill_between(age,mgh[:,2]+exy1[:,2],mgh[:,2]-exy2[:,2],alpha='0.20',color="r")

    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'LGH_'+name.replace('_L','')+'.pdf',dpi = 1000)        
    else:
        plt.savefig(dir1+'MGH_'+name+'.pdf',dpi = 1000)
    plt.close()
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)

    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$SFHL(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFHL INPUT")
    else:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFH INPUT")
    plt.semilogx(age,sfh[:,0],'',color="b",drawstyle='steps',label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=2)
    plt.semilogx(age,sfh[:,1],'--',color="g",drawstyle='steps',label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=2)
    plt.semilogx(age,sfh[:,2],':',color="r",drawstyle='steps',label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'SFHL_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'SFH_'+name+'.pdf',dpi = 1000)
    plt.close()
    
def ensamble_int(name, dir='', dir1='',lig=0):
#    file2=dir1+name.replace('simu','photo').replace(".fits.gz",".r_Lc_rad.fits")
#    if ptt.exists(file2) == True:
#    [pdl_rad, hdr2]=gdata(file2,0, header=True)
    ind=[]
    inda=[]
    rad=12.0
    vel_light=299792.458
    [mass_age_t,hdr]=gdata(dir1+name, 0, header=True)
    dpix=hdr['CD1_1']*3600.0
    [n_age,nx,ny]=mass_age_t.shape
    pdl_rad=np.zeros([nx,ny])
    for i in range(0, nx):
        for j in range(0, ny):
            pdl_rad[i,j]=np.sqrt((i-nx/2.)**2.0+(j-ny/2.0)**2.0)*dpix/(2.5*0.5)
    nt=np.where(pdl_rad< rad)
    ages=np.zeros(n_age)
    for i in range(0, n_age):
        ages[i]=np.log10(hdr["AGE"+str(i)])+9.0
    mass_age_t_2=mass_age_t
    mass_age_t_e=np.zeros([n_age,nx,ny])
    MassT=0
    massN=0
    massN_e=0
    mass=np.zeros([n_age])
    mass_e=np.zeros([n_age])
    sfrt=np.zeros([n_age])
    sfdt=np.zeros([n_age])
    mass_temp_total=0
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
        if i == 0:
            temp1=temp
        else:
            temp1=temp1+temp
        MassT=MassT+np.sum(temp)
        
        Dt_mass=np.sum(temp[nt])
        Dt_mass_e=np.sum(temp_e[nt])
        Dt_mass_2=np.sum(temp_2[nt])
        massN=Dt_mass+massN
        massN_e=Dt_mass_e+massN_e
        mass[i]=np.log10(massN)
        mass_e[i]=massN_e
        sfrt[i]=Dt_mass_2/Dt_age
#    print mass[n_age-1,:]
        #print massN
        #print nr
#    sys.exit()
    mass_temp_total=np.log10(10**mass[n_age-1])
    MassT=np.log10(MassT)
    MassT=mass_temp_total
    mass_n=10**(10**(mass-mass[n_age-1]))
    mass_n_e=np.sqrt((10**(mass-mass[n_age-1]))**2.0*((mass_e/10**(2.0*mass))+(mass_e[n_age-1]/10**(2.0*mass[n_age-1]))))
    if lig == 1:
        f2=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_int_L.csv","w")
        f2.write("#  LOG_AGE  N_LIGHR_1  N_LIGHR_2  N_LIGHR_3  N_LIGHR_4  LOG_LIGHR_1  LOG_LIGHR_2  LOG_LIGHR_3  LOG_LIGHR_4 \n")
    else:
        f2=open(dir1+name.replace(".fits.gz","")+"_Ensemble_int.csv","w")
        f2.write("#  LOG_AGE  N_MASSR_1  N_MASSR_2  N_MASSR_3  N_MASSR_4  LOG_MASSR_1  LOG_MASSR_2  LOG_MASSR_3  LOG_MASSR_4 \n")
    for i in range(0, n_age):
        line=''
        line=line+str(ages[i])
        line=line+';'+str(mass_n[i])
        line=line+';'+str(mass[i])
        line=line+';'+str(sfrt[i])
        line=line+';'+str(mass_n_e[i])
        line=line+';'+str(sfdt[i])
        line=line+' \n'
        f2.write(line)
    f2.close()
    
    
    nrb=6
    ensamb=np.zeros([n_age,nrb])
    co=0
    if lig == 1:
        f=open(dir1+name.replace("_L.fits.gz","")+"_Ensemble_int_L.csv","r")
    else:
        f=open(dir1+name.replace(".fits.gz","")+"_Ensemble_int.csv","r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1

    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1])
    sfh=ensamb[:,3]#/10.0**ensamb[:,4:7]
    emgh=(ensamb[:,4])#*10
    #print ensamb[:,7:10]
    #print 10.0**ensamb[:,4:7]
    dth=ensamb[:,5]/1e3
    Ddth=np.zeros(len(dth))
    nx=len(mgh)
    Dmgh=np.zeros([nx])
    exy1=np.zeros([nx])
    exy2=np.zeros([nx])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i]-mgh[i+1])
            Ddth[i]=np.abs(dth[i]-dth[i+1])#/2.0
        elif i == len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i-1]-mgh[i])
            Ddth[i]=np.abs(dth[i-1]-dth[i])#/2.0 
    exy1=np.amin([np.abs(Dmgh/(dth-Ddth)*dth)+emgh,np.abs(Dmgh/(dth+Ddth)*dth)+emgh],axis=0)#np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]
    exy2=np.amin([np.abs(Dmgh/(dth-Ddth)*dth)+emgh,np.abs(Dmgh/(dth+Ddth)*dth)+emgh],axis=0)#np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0]

    val1=np.abs(Dmgh/(dth+Ddth)*dth)+emgh
    val2=np.abs(Dmgh/(dth-Ddth)*dth)+emgh
    for j in range(0, len(val1)):
        if val1[j] >= 0.15:
            val1[j]=0.1
        if val2[j] >= 0.15:
            val2[j]=0.1
    exy1=val1
    exy2=val2
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    #plt.title(name)
    if lig == 1:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title("SIMULATION LGH INT INPUT")
    else:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title("SIMULATION MGH INT INPUT")
    plt.semilogx(age,mgh,'',color="b",label='$R=15arcsec$',lw=1.5)
    ax.fill_between(age,mgh+exy1,mgh-exy2,alpha='0.20',color="b")
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'LGH_int_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'MGH_int_'+name+'.pdf',dpi = 1000)
    plt.close()
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)

    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 1:
        plt.ylabel(r'$SFHL(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFHL INT INPUT")
    else:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title("SIMULATION SFH INT INPUT")
    plt.semilogx(age,sfh,'',color="b",drawstyle='steps',label='$R=15arcsec$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 1:
        plt.savefig(dir1+'SFHL_int_'+name.replace('_L','')+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir1+'SFH_int_'+name+'.pdf',dpi = 1000)
    plt.close()
    
def band_photo_r(pdl_flux, hdr, name='test.fits.gz', dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt']#,'U_Johnson.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt','J_2MASS.txt','H_2MASS.txt','K_2MASS.txt']
    band=['u','g','rg','ig','z']#,'U','B','V','I','R','J','H','K']
    zerop=[3631.0,3730.0,4490.0,4760.0,4810.0]#,1810.0,4260.0,3640.0,2550.0,3080.0,1600.0,1080.0,670.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                flux_t=spec1*trans
                #photo_a[i,j]=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2)/simpson_r(trans,wave,0,len(wave)-2)
                photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                photo_b[i,j]=f_fin*zerop[k]*jans
                photo_c[i,j]=f_fi2
#                if photo_a[i,j] > 0:
#                    print photo_a[i,j]
#                    import matplotlib.pyplot as pl
#                    pl.plot(wave,trans)
#                    pl.show()
#                    pl.plot(wave,flux_t)
#                    pl.show() 
#                    pl.plot(wave,spec1)
#                    pl.plot(wave_s,spec)
#                    pl.show()
#                    print simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)
#                    print simpson_r(trans,wave,0,len(wave)-2,typ=1)
#                    print simpson_r((trans/(wave**2.0)*vel_light/ang),wave,0,len(wave)-2)
#                    print simpson_r(trans/((wave*ang)**2.0)*vel_light,wave*ang,0,len(wave)-2)
#                    print simpson_r(trans,vel_light/(wave*ang),0,len(wave)-2)
#                    sys.exit()
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')
        imag_F2=photo_b#convolve(photo_b, PSF)
        imag_F3=photo_c#convolve(photo_c, PSF)
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr2)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    
def band_photo(name, dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [pdl_flux,hdr]=gdata(dir1+name, 0, header=True)
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','NUV_GALEX.txt','FUV_GALEX.txt','U_Johnson.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt','J_2MASS.txt','H_2MASS.txt','K_2MASS.txt']
    band=['ha','u','g','rg','ig','z','NUV','FUV','U','B','V','I','R','J','H','K']
#    zerop=[3631.0,3631.0,3730.0,4490.0,4760.0,4810.0,3801.4,3619.9,1810.0,4260.0,3640.0,2550.0,3080.0,1600.0,1080.0,670.0]
    zerop=[3631.0,3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                flux_t=spec1*trans*d_wave
                f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                photo_b[i,j]=f_fin*zerop[k]*jans
                photo_c[i,j]=f_fi2
#                if photo_a[i,j] > 0:
#                    print photo_a[i,j]
#                    import matplotlib.pyplot as pl
#                    pl.plot(wave,trans)
#                    pl.show()
#                    pl.plot(wave,flux_t)
#                    pl.show() 
#                    pl.plot(wave,spec1)
#                    pl.plot(wave_s,spec)
#                    pl.show()
#                    print simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)
#                    print simpson_r(trans,wave,0,len(wave)-2,typ=1)
#                    print simpson_r((trans/(wave**2.0)*vel_light/ang),wave,0,len(wave)-2)
#                    print simpson_r(trans/((wave*ang)**2.0)*vel_light,wave*ang,0,len(wave)-2)
#                    print simpson_r(trans,vel_light/(wave*ang),0,len(wave)-2)
#                    sys.exit()
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
            del hdr["CUNIT3"]
            del hdr['RADECSYS']
            hdr['RADECSYSa']='ICRS    '
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2/A'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2=photo_b#convolve(photo_b, PSF)#photo_b#
        imag_F3=photo_c#convolve(photo_c, PSF)#photo_c#
        imag_F1C=convolve(photo_a, PSF,  boundary='symm')#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2C=convolve(photo_b, PSF)#photo_b#
        imag_F3C=convolve(photo_c, PSF)#photo_c#
        imag_F1C[np.where(imag_F1C == 0)]=35.0
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr2)
        wfits(dir1+name_f+band[k]+'_c.fits',imag_F1C,hdr)
        wfits(dir1+name_f+band[k]+'_Fc.fits',imag_F2C,hdr2)
        wfits(dir1+name_f+band[k]+'_Lc.fits',imag_F3C,hdr2)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#5e-15)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)#7e-23)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    from PIL import Image
    max=20.0
    min=35.0
    file00u=dir1+name_f+"u_c.fits"
    file00g=dir1+name_f+"g_c.fits"
    file00r=dir1+name_f+"rg_c.fits"
    file00i=dir1+name_f+"ig_c.fits"
    file00z=dir1+name_f+"z_c.fits"
    file00U=dir1+name_f+"U_c.fits"
    file00B=dir1+name_f+"B_c.fits"
    file00V=dir1+name_f+"V_c.fits"
    file00R=dir1+name_f+"R_c.fits"
    file00I=dir1+name_f+"I_c.fits"    
    [pdl_00u,hdrt]=gdata(file00u, 0, header=True)
    [pdl_00g,hdrt]=gdata(file00g, 0, header=True)
    [pdl_00r,hdrt]=gdata(file00r, 0, header=True)
    [pdl_00i,hdrt]=gdata(file00i, 0, header=True)
    [pdl_00z,hdrt]=gdata(file00z, 0, header=True)
    [pdl_00U,hdrt]=gdata(file00U, 0, header=True)
    [pdl_00B,hdrt]=gdata(file00B, 0, header=True)
    [pdl_00V,hdrt]=gdata(file00V, 0, header=True)
    [pdl_00R,hdrt]=gdata(file00R, 0, header=True)
    [pdl_00I,hdrt]=gdata(file00I, 0, header=True)
    nx,ny=pdl_00g.shape
    pdl_00u=(np.flipud(pdl_00u)-min)/(max-min)*256
    pdl_00g=(np.flipud(pdl_00g)-min)/(max-min)*256
    pdl_00r=(np.flipud(pdl_00r)-min)/(max-min)*256
    pdl_00i=(np.flipud(pdl_00i)-min)/(max-min)*256
    pdl_00z=(np.flipud(pdl_00z)-min)/(max-min)*256
    pdl_00U=(np.flipud(pdl_00U)-min)/(max-min)*256
    pdl_00B=(np.flipud(pdl_00B)-min)/(max-min)*256
    pdl_00V=(np.flipud(pdl_00V)-min)/(max-min)*256
    pdl_00R=(np.flipud(pdl_00R)-min)/(max-min)*256
    pdl_00I=(np.flipud(pdl_00I)-min)/(max-min)*256
    pdl_00u[np.where(pdl_00u < 0)]=0
    pdl_00g[np.where(pdl_00g < 0)]=0
    pdl_00r[np.where(pdl_00r < 0)]=0
    pdl_00i[np.where(pdl_00i < 0)]=0
    pdl_00z[np.where(pdl_00z < 0)]=0
    pdl_00U[np.where(pdl_00U < 0)]=0
    pdl_00B[np.where(pdl_00B < 0)]=0
    pdl_00V[np.where(pdl_00V < 0)]=0
    pdl_00R[np.where(pdl_00R < 0)]=0
    pdl_00I[np.where(pdl_00I < 0)]=0
    pdl_00u[np.where(pdl_00u > 255)]=255
    pdl_00g[np.where(pdl_00g > 255)]=255
    pdl_00r[np.where(pdl_00r > 255)]=255
    pdl_00i[np.where(pdl_00i > 255)]=255
    pdl_00z[np.where(pdl_00z > 255)]=255
    pdl_00U[np.where(pdl_00U > 255)]=255
    pdl_00B[np.where(pdl_00B > 255)]=255
    pdl_00V[np.where(pdl_00V > 255)]=255
    pdl_00R[np.where(pdl_00R > 255)]=255
    pdl_00I[np.where(pdl_00I > 255)]=255
    pdl_00=np.zeros([nx,ny,3],dtype="uint8")
    pdl_00[:,:,0]=pdl_00i
    pdl_00[:,:,1]=pdl_00r
    pdl_00[:,:,2]=pdl_00g
    im = Image.fromarray(pdl_00)
    im.save(dir1+name_f+"gri.jpeg")
    pdl_00[:,:,0]=pdl_00R
    pdl_00[:,:,1]=pdl_00V
    pdl_00[:,:,2]=pdl_00B
    im1 = Image.fromarray(pdl_00)
    im1.save(dir1+name_f+"BVR.jpeg")
    
def band_cube(name, dir='', dir1=''):
    vel_light=299792458.0
    ang=1e-10
    jans=1e-23
    [pdl_flux,hdr]=gdata(dir1+name, 0, header=True)
    #pdl_noise=gdata(dir1+name, 1)
    pdl_flux=pdl_flux*1e-16
    #pdl_noise=pdl_noise*1e-16
    [nw,nx,ny]=pdl_flux.shape
    crpix=hdr["CRPIX3"]
    cdelt=hdr["CDELT3"]
    crval=hdr["CRVAL3"]
    int_spec1=np.zeros(nw)
    int_spec2=np.zeros(nw)
    wave_s=np.zeros(nw)
    for j in range(0, nw):
        wave_s[j]=crval+cdelt*(j+1-crpix)
        int_spec1[j]=np.sum(pdl_flux[j,:,:])#/cdelt
        int_spec2[j]=np.sum(pdl_flux[j,:,:])*wave_s[j]**2.0/cdelt/vel_light*ang
    file=['ha_filter_sh.txt','SDSS_u.txt','SDSS_g.txt','SDSS_r.txt','SDSS_i.txt','SDSS_z.txt','B_Johnson.txt','V_Johnson.txt','I_Cousins.txt','R_Cousins.txt']
    band=['ha','u','g','rg','ig','z','B','V','I','R']
#    zerop=[3631.0,3631.0,3730.0,4490.0,4760.0,4810.0,4260.0,3640.0,2550.0,3080.0]
    zerop=[3631.0,3631.0,3730.0,3730.0,3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]#3631.0,3631.0,3631.0,3631.0,3631.0,3631.0]
    for k in range(0, len(band)):
        photo_a=np.zeros([nx,ny])
        photo_b=np.zeros([nx,ny])
        photo_c=np.zeros([nx,ny])
        photo_ae=np.zeros([nx,ny])
        photo_be=np.zeros([nx,ny])
        photo_ce=np.zeros([nx,ny])
        f=open(dir+file[k],'r')
        wave=[]
        trans=[]
        for line in f:
            line=line.replace('\n','')
            data=line.split(' ')
            data=filter(None,data)
            if len(data) > 1:
                wave.extend([float_(data[1])])
                trans.extend([float_(data[2])])
        f.close()
        d_wave=np.zeros(len(wave))
        for kk in range(1,len(wave)):
            d_wave[kk]=wave[kk]-wave[kk-1]
        d_wave[0]=d_wave[1]
        trans=np.array(trans)
        wave=np.array(wave)
        for i in range(0, nx):
            for j in range(0, ny):
                spec=pdl_flux[:,i,j]
                #nois=pdl_noise[:,i,j]
                if np.sum(spec) > 0:
                    spec1=interp1d(wave_s, spec,bounds_error=False,fill_value=0.)(wave)
                    #nois1=interp1d(wave_s, nois,bounds_error=False,fill_value=0.)(wave)
                    flux_t=spec1*trans*d_wave
                    #nois_t=nois1*trans
                    f_fin=simpson_r(flux_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                    f_fi2=simpson_r(flux_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                    #f_fin_e=simpson_r(nois_t*wave**2.0/d_wave/vel_light*ang,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)/jans/zerop[k]
                    #f_fi2_e=simpson_r(nois_t/d_wave,wave,0,len(wave)-2,typ=1)/simpson_r(trans,wave,0,len(wave)-2,typ=1)
                    #photo_a[i,j]=-2.5*np.log10(f_fin+1e-14)
                    if f_fin <= 0:
                        photo_a[i,j]=-2.5*np.log10(1e-10)#14
                    else:
                        photo_a[i,j]=-2.5*np.log10(f_fin+1e-10)
                    photo_b[i,j]=f_fin*zerop[k]*jans
                    photo_c[i,j]=f_fi2
                    #photo_ae[i,j]=-2.5*np.log10(f_fin_e+1e-12)
                    #photo_be[i,j]=f_fin_e*zerop[k]*jans
                    #photo_ce[i,j]=f_fi2_e
        if k == 0:
            hdr["NAXIS"]=2
            hdr["UNITS"]=('Magnitudes', '-2.5log(F/F0) with F0 as ZPOINT')
            del hdr["NAXIS3"]
            del hdr["CRPIX3"]
            del hdr["CDELT3"]
            del hdr["CRVAL3"]
        hdr["PBAND"]=file[k].replace('.txt','')
        hdr["ZPOINT"]=(zerop[k], 'Zero Point in Jy')
        hdr2=hdr
        hdr3=hdr
        if k == 0:
            hdr2["UNITS"]='ergs/s/cm^2/Hz'
            hdr3["UNITS"]='ergs/s/cm^2/A'
        dv=2
        PSF=Gaussian2DKernel(stddev=dv)
        imag_F1=photo_a#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        imag_F2=photo_b#convolve(photo_b, PSF)#photo_b#
        imag_F3=photo_c#convolve(photo_c, PSF)#photo_c#
        #imag_E1=photo_ae#convolve(photo_a, PSF)#, mode='full')#,  boundary='symm')#photo_a#
        #imag_E2=photo_be#convolve(photo_b, PSF)#photo_b#
        #imag_E3=photo_ce#convolve(photo_c, PSF)#photo_c#
        name_f=name.replace('.fits.gz','.')
        wfits(dir1+name_f+band[k]+'.fits',imag_F1,hdr)
        wfits(dir1+name_f+band[k]+'_F.fits',imag_F2,hdr2)
        wfits(dir1+name_f+band[k]+'_L.fits',imag_F3,hdr3)
        #wfits(dir1+name_f+band[k]+'e.fits',imag_E1,hdr)
        #wfits(dir1+name_f+band[k]+'_Fe.fits',imag_E2,hdr2)
        #wfits(dir1+name_f+band[k]+'_Le.fits',imag_E3,hdr3)
    name_f=name.replace('.fits.gz','.')
    nt_s=np.where((wave_s > 4500) & (wave_s < 10000))[0]
    max_val1=np.amax(int_spec1[nt_s])*1.25
    max_val2=np.amax(int_spec2[nt_s])*1.25
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val1)#0.75e-14)#
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/A]$",fontsize=14)
    plt.plot(wave_s,int_spec1)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec.pdf')
    plt.close()
    fig, ax = plt.subplots(figsize=(8,5.5))
    ax.set_ylim(0,max_val2)
    ax.set_xlabel("$Wavelength [A]$",fontsize=14)
    ax.set_ylabel("Flux $[erg/s/cm^2/Hz]$",fontsize=14)
    plt.plot(wave_s,int_spec2)
    fig.tight_layout()
    plt.savefig(dir1+name_f+'spec_Hz.pdf')
    plt.close()
    from PIL import Image
    file00u=dir1+name_f+"u.fits"
    file00g=dir1+name_f+"g.fits"
    file00r=dir1+name_f+"rg.fits"
    file00i=dir1+name_f+"ig.fits"
    file00z=dir1+name_f+"z.fits"
    file00B=dir1+name_f+"B.fits"
    file00V=dir1+name_f+"V.fits"
    file00R=dir1+name_f+"R.fits"
    file00I=dir1+name_f+"I.fits"    
    [pdl_00u,hdrt]=gdata(file00u, 0, header=True)
    [pdl_00g,hdrt]=gdata(file00g, 0, header=True)
    [pdl_00r,hdrt]=gdata(file00r, 0, header=True)
    [pdl_00i,hdrt]=gdata(file00i, 0, header=True)
    [pdl_00z,hdrt]=gdata(file00z, 0, header=True)
    [pdl_00B,hdrt]=gdata(file00B, 0, header=True)
    [pdl_00V,hdrt]=gdata(file00V, 0, header=True)
    [pdl_00R,hdrt]=gdata(file00R, 0, header=True)
    [pdl_00I,hdrt]=gdata(file00I, 0, header=True)
    max=np.amin(pdl_00r[np.where(pdl_00r > 0)])-0.5#24.0#17.5
    min=25.0
    nx,ny=pdl_00g.shape
    pdl_00u=(np.flipud(pdl_00u)-min)/(max-min)*256
    pdl_00g=(np.flipud(pdl_00g)-min)/(max-min)*256
    pdl_00r=(np.flipud(pdl_00r)-min)/(max-min)*256
    pdl_00i=(np.flipud(pdl_00i)-min)/(max-min)*256
    pdl_00z=(np.flipud(pdl_00z)-min)/(max-min)*256
    pdl_00B=(np.flipud(pdl_00B)-min)/(max-min)*256
    pdl_00V=(np.flipud(pdl_00V)-min)/(max-min)*256
    pdl_00R=(np.flipud(pdl_00R)-min)/(max-min)*256
    pdl_00I=(np.flipud(pdl_00I)-min)/(max-min)*256
    pdl_00u[np.where(pdl_00u < 0)]=0
    pdl_00g[np.where(pdl_00g < 0)]=0
    pdl_00r[np.where(pdl_00r < 0)]=0
    pdl_00i[np.where(pdl_00i < 0)]=0
    pdl_00z[np.where(pdl_00z < 0)]=0
    pdl_00B[np.where(pdl_00B < 0)]=0
    pdl_00V[np.where(pdl_00V < 0)]=0
    pdl_00R[np.where(pdl_00R < 0)]=0
    pdl_00I[np.where(pdl_00I < 0)]=0
    pdl_00u[np.where(pdl_00u > 255)]=255
    pdl_00g[np.where(pdl_00g > 255)]=255
    pdl_00r[np.where(pdl_00r > 255)]=255
    pdl_00i[np.where(pdl_00i > 255)]=255
    pdl_00z[np.where(pdl_00z > 255)]=255
    pdl_00B[np.where(pdl_00B > 255)]=255
    pdl_00V[np.where(pdl_00V > 255)]=255
    pdl_00R[np.where(pdl_00R > 255)]=255
    pdl_00I[np.where(pdl_00I > 255)]=255
    pdl_00=np.zeros([nx,ny,3],dtype="uint8")
    pdl_00[:,:,0]=pdl_00i
    pdl_00[:,:,1]=pdl_00r
    pdl_00[:,:,2]=pdl_00g
    im = Image.fromarray(pdl_00)
    im.save(dir1+name_f+"gri.jpeg",quality=100)
    pdl_00[:,:,0]=pdl_00R
    pdl_00[:,:,1]=pdl_00V
    pdl_00[:,:,2]=pdl_00B
    im1 = Image.fromarray(pdl_00)
    im1.save(dir1+name_f+"BVR.jpeg",quality=100)
    
def id_str(id):
    id=int(np.float_(id))
    if id < 10:
        idt='0'+str(id)
        #idt=str(id)
    else:
        idt=str(id)
    return idt


def mock_test(fib_n,typef1="MaNGA",id1='A2-0',psf=0,dir1=''):
    sig=2.5
    thet=0.0
    plots=1
    nl=110
    n_pix=440
    cam=(142135.5*2.0)*7.0/np.float_(fib_n)#/2.34 # 2.34 sp6#/3.0#5.5#/2.2#IF CALIFA 5.5
    fov_p=np.round(142135.5/np.abs(cam)*91.0)
    fov=np.round(142135.5/np.abs(cam)*62.0)#30
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=300.0#150.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name='test-'+id1
    print name
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    Om=0.2726
    Lam=0.7274
    ho=0.704
    mock_halo_test(id1,psf=psf,dir1=dir1,fib_n=fib_n,nl=nl,fov=fov,thet=thet,ifutype=typef1)
    


def mock_sp(fib_n,ang,sp_res=2000.0,sp_samp=1.25,modt=0,template_1="libs/gsd61_156.fits",template_2="libs/templete_gas.fits",template_3="libs/templete_bc03_5.fits",template="libs/templete_bc03_2.fits",n_pix=440,fov_p=0,fov=0,rx=142135.5,Om=0.2726,Lam=0.7274,ho=0.704,cam=0,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',typef1="MaNGA",id1='A2-0',psf=0,redo=0,SN=15.0,Fluxm=20.0,dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat",file_out='mock_mass_ill_0.out',file_out_f='mock_mass_ill.out'):
    thet=0.0
    plots=1
    nl=110
    if cam == 0:
        if "CALIFA" in typef1:
            fib_n=14
        cam=(rx*2.0)*7.0/np.float_(fib_n)#/2.34 # 2.34 sp6#/3.0#5.5#/2.2#IF CALIFA 5.5
    cam=-cam
    if fov_p == 0:
        fov_p=np.round(rx/np.abs(cam)*91.0)
    if fov == 0:
        fov=np.round(rx/np.abs(cam)*62.0)#30
    if "MaNGA" in typef1:
        if psf <= 0:
            sig=1.43
        else:
            sig=psf
    elif "CALIFA" in typef1:
        if psf <= 0:
            sig=0.7
        else:
            sig=psf
    elif "MUSE" in typef1:
        if psf <= 0:
            sig=0.6
        else:
            sig=psf
    else:
        if psf == 0:
            sig=1.43
        else:
            sig=psf  
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=300.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name=base_name+id1
    #print name
    sycall('echo '+name)
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    R3=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]])
    ex=np.array([1,0,0])
    Ev=np.transpose(np.array([vx,vy,vz]))
    obs_ang1=np.dot(np.dot(Ev,R3),ex)*cam 
    if modt == 0:
        mock_photosim(id1,basename=base_name,file=file_red,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_photosimextgas(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_photo(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_sim(id1,basename=base_name,file=file_red,template2=template,template=template_1,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    if modt == 1:
        mock_photo(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    if modt == 2:
        mock_halo(id1,sp_res=sp_res,sp_samp=sp_samp,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
    if modt == 3:
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    #fit3d(name,fo,fot_par,dir2=dir1,reens=1)#,redo=1)
    #fit3d_cen(name,fa,reens=0)

    
def mock_cen_re(id1='A2-12700',dir1='',file_out='mock_mass_ill_0.out',n_mt=500,n_mt0=0):
    if ptt.exists(dir1+id1) == False:
        print "The File does not exist"
        sys.exit()
    name='artsp8-'+id1
    print name
    fa=open('FIT3D_c_std'+file_out,'w')                     
    fit3d_cen_stat(name,fa,n_m=n_mt,n_mo=n_mt0)
    fa.close()
    
    
def mock_re(fib_n,ang,template_2="libs/templete_gas.fits",template="libs/templete_bc03_2.fits",ns=0,fov_p=0,fov=0,rx=142135.5,n_pix=440,cam=0,Om=0.2726,Lam=0.7274,ho=0.704,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',typef1="MaNGA",id1='A2-0',redo=0,SN=15.0,dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat"):
    sig=2.5
    thet=0.0
    plots=1
    nl=110
    if cam == 0:
        cam=(rx*2.0)*7.0/np.float_(fib_n)#/2.2
    cam=-cam
    if fov_p == 0:
        fov_p=np.round(rx/np.abs(cam)*91.0)
    if fov == 0:
        fov=np.round(rx/np.abs(cam)*30.0)
    if "CALIFA" in typef1:
        fib_n=11
    elif "MUSE" in typef1:
        scp_s=150.0
        fibA=150.0
        fib_n=int(fov*scp_s/fibA/2)+1
    if ns == 0:
        ns=str(int(3*fib_n*(fib_n-1)+1))
    id1t=id1.split('-')
    if len(id1t) == 1:
        id=0
        while True:
            if ptt.exists(dir1+id1t[0]+'-'+ns+id_str(id)) == False:
                break
            else:
                id=id+1
        id1=id1t[0]+'-'+ns+id_str(id)
    else:
        id1=id1t[0]+'-'+ns+id_str(id1t[1])
    name=base_name+id1
    print name
    fov1=0.06
    rads=[0,0.5,1.0,1.5]
    #Om=0.2726
    #Lam=0.7274
    #ho=0.704
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    #f4=open('Ntemp','w')
    R3=np.array([[np.cos(ang*np.pi/180.0),-np.sin(ang*np.pi/180.0),0],[np.sin(ang*np.pi/180.0),np.cos(ang*np.pi/180.0),0],[0,0,1]])
    ex=np.array([1,0,0])
    Ev=np.array([vx,vy,vz])
    obs_ang1=np.dot(np.dot(Ev,R3),ex)*cam  
    #obs_ang1=[(+0.8972639*np.cos(ang*np.pi/180.0)+0.4414944*np.sin(ang*np.pi/180.0))*cam,(-0.07978683*np.cos(ang*np.pi/180.0)+0.1621534*np.sin(ang*np.pi/180.0))*cam,(-0.4342251*np.cos(ang*np.pi/180.0)+0.8824902*np.sin(ang*np.pi/180.0))*cam] 
    #obs_ang1=[(-0.7156755*np.cos(ang*np.pi/180.0)+0.6984330*np.sin(ang*np.pi/180.0))*cam,(-0.5130859*np.cos(ang*np.pi/180.0)-0.5257526*np.sin(ang*np.pi/180.0))*cam,(0.4738687*np.cos(ang*np.pi/180.0)+0.4855672*np.sin(ang*np.pi/180.0))*cam]
    mock_photosim(id1,file=file_red,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
    mock_photosimextgas(id1,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
    #f4.close()
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
    
def mock_halo_ill(idh,idhp,dir1='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="MaNGA",basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
    cubef='ilust-'+str(idhp)+'-'+id+'.cube'
    dirf=dir1t+'spec_lin/'
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    stars = il.snapshot.loadSubhalo(basePath,135,idh,'stars')
    gas = il.snapshot.loadSubhalo(basePath,135,idh,'gas')
    dx_g = gas['Coordinates'][:,0]#
    dy_g = gas['Coordinates'][:,1]#
    dz_g = gas['Coordinates'][:,2]#
    vx_g = gas['Velocities'][:,0]#
    vy_g = gas['Velocities'][:,1]#
    vz_g = gas['Velocities'][:,2]#
    volm = gas['Volume'][:]
    dens = gas['Density'][:]#
    sfri = gas['StarFormationRate'][:]#
    meta_g=gas['GFM_Metallicity'][:]#     
    dx = stars['Coordinates'][:,0]#
    dy = stars['Coordinates'][:,1]#
    dz = stars['Coordinates'][:,2]#
    phot=stars['GFM_StellarPhotometrics'][:,5]
    mass=stars['Masses'][:]#
    mas0=stars['GFM_InitialMass'][:]#
    meta=stars['GFM_Metallicity'][:]#
    denS=stars['SubfindDensity']
    time=stars['GFM_StellarFormationTime'][:]#
    vx = stars['Velocities'][:,0]#
    vy = stars['Velocities'][:,1]#
    vz = stars['Velocities'][:,2]#
    mass_g=gas['Masses'][:]
    TE =   gas['InternalEnergy'][:]     
    #ages=stars['GFM_StellarAge'][:]
    temp_g=TE*((1e5)**2.0)*(0.938*1.7827e-24)*(4.0/(8.0-5.0*0.245))/(1.3807e-16)
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    phot=phot[nt]
    mass_g=mass_g/ho*1e10
    mass=mass[nt]/ho*1e10
    mas0=mas0[nt]/ho*1e10
    meta=meta[nt]#/0.0127
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=volm/ho**3.0
    dens=dens*ho**2.0
    volm=float_((volm/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    dens=dens*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    Av_g=meta_g*(3.0*1.67262e-24*np.pi*dens*volm)/(4.0*np.log(10.)*3.0*5494e-8)
   # print meta_g,(10.0**(-0.59)*0.0127)
    nt1=np.where(meta_g > (10.0**(-0.59)*0.0127))[0]
    nt2=np.where(meta_g <= (10.0**(-0.59)*0.0127))[0]
    Av_g[nt1]=Av_g[nt1]/(10.0**(2.21-1.0))
    if len(nt2) > 0 :
        Av_g[nt2]=Av_g[nt2]/(10.0**(2.21-1.0)/(meta_g[nt2]/0.0127)**(3.1-1.0))
    #print np.amax(meta),np.amin(meta)
    zf=1/time[nt]-1
    xo=modes(dx,nit=7)
    yo=modes(dy,nit=7)
    zo=modes(dz,nit=7)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    #print np.dot(np.dot(R2,R1),np.array([[x0],[y0],[z0]]))
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    print red_0
    print len(zf)
    age_s=cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        cube_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=fib_n,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
        band_cube(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
    else:
        #band_cube(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        sycall('cp '+dir1n+cubef+'_val.fits.gz '+dirf)
  
def mock_photo_ill(idh,idhp,f4,dir1='',ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
    cubef='ilust-'+str(idhp)+'-'+id+'.photo'
    dirf=dir1t+'photo_lin/'   
    sycall('mkdir -p '+dirf)
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    stars = il.snapshot.loadSubhalo(basePath,135,idh,'stars')
    gas = il.snapshot.loadSubhalo(basePath,135,idh,'gas')
    dx_g = gas['Coordinates'][:,0]#
    dy_g = gas['Coordinates'][:,1]#
    dz_g = gas['Coordinates'][:,2]#
    vx_g = gas['Velocities'][:,0]#
    vy_g = gas['Velocities'][:,1]#
    vz_g = gas['Velocities'][:,2]#
    volm = gas['Volume'][:]
    dens = gas['Density'][:]#
    sfri = gas['StarFormationRate'][:]#
    meta_g=gas['GFM_Metallicity'][:]#     
    dx = stars['Coordinates'][:,0]#
    dy = stars['Coordinates'][:,1]#
    dz = stars['Coordinates'][:,2]#
    phot=stars['GFM_StellarPhotometrics'][:,5]
    mass=stars['Masses'][:]#
    mas0=stars['GFM_InitialMass'][:]#
    meta=stars['GFM_Metallicity'][:]#
    denS=stars['SubfindDensity']
    time=stars['GFM_StellarFormationTime'][:]#
    vx = stars['Velocities'][:,0]#
    vy = stars['Velocities'][:,1]#
    vz = stars['Velocities'][:,2]#
    mass_g=gas['Masses'][:]
    TE =   gas['InternalEnergy'][:]     
    #ages=stars['GFM_StellarAge'][:]
    temp_g=TE*((1e5)**2.0)*(0.938*1.7827e-24)*(4.0/(8.0-5.0*0.245))/(1.3807e-16)
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    phot=phot[nt]
    mass_g=mass_g/ho*1e10
    mass=mass[nt]/ho*1e10
    mas0=mas0[nt]/ho*1e10
    meta=meta[nt]#/0.0127
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=volm/ho**3.0
    dens=dens*ho**2.0
    volm=float_((volm/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    dens=dens*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    Av_g=meta_g*(3.0*1.67262e-24*np.pi*dens*volm)/(4.0*np.log(10.)*3.0*5494e-8)
   # print meta_g,(10.0**(-0.59)*0.0127)
    nt1=np.where(meta_g > (10.0**(-0.59)*0.0127))[0]
    nt2=np.where(meta_g <= (10.0**(-0.59)*0.0127))[0]
    Av_g[nt1]=Av_g[nt1]/(10.0**(2.21-1.0))
    if len(nt2) > 0 :
        Av_g[nt2]=Av_g[nt2]/(10.0**(2.21-1.0)/(meta_g[nt2]/0.0127)**(3.1-1.0))
    #print np.amax(meta),np.amin(meta)
    zf=1/time[nt]-1
    xo=modes(dx,nit=7)
    yo=modes(dy,nit=7)
    zo=modes(dz,nit=7)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
#    print Ra,"HOLA 2",x0,y0,z0,red_0
#    sys.exit()
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
#    print modes(x,nit=7,n_s=0.8),modes(y,nit=7,n_s=0.8),modes(z,nit=7,n_s=0.8)
#    sys.exit()
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    #print np.dot(np.dot(R2,R1),np.array([[x0],[y0],[z0]]))
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    print red_0
    print len(zf)
    age_s=cd.lookback_time(zf, **cosmo)/31557600./1e9
    if ptt.exists(dir1+cubef+'.fits.gz') == False:
        [age_F,ind]=ages_definition(age_s,n_ages=55)
        photo_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        band_photo(cubef+'.fits.gz', dir='legacy/', dir1=dir1)
        [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=1)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    else:
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        #band_photo(cubef+'.fits.gz', dir='/home/hjibarram/FIT3D_py/soft_f/legacy/', dir1=dir1)
        if ptt.exists(dir1+cubef+".r_Lc_rad.fits") == True:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=galfit_param(cubef,band=".rg_Lc",psf=0.2,dx=nl/2,dy=nl/2,dir=dir1,dir2=dir1,repro=0)
        else:
            [e_rad,e_rad_pix,n_ser,ab,PA,not_fit]=[-100,-100,-100,-100,-100,-100]
        #sys.exit()
        hdulist = fits.open(dir1+cubef+".rg_Lc.fits")
        hd=hdulist[0].header
        red_0=hd["REDSHIFT"]
        cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
        cosmo = cd.set_omega_k_0(cosmo)
        cam=cd.comoving_distance(red_0, **cosmo)*1e3
        R50=e_rad/3600.*cam/(1+red_0)*np.pi/180.
        sycall("mv "+cubef+".* "+dir1n)
        sycall("mv fit.log "+dir1n)
        print R50, e_rad, not_fit
    if not_fit == 0:
        ifl=(e_rad*2.1/2.5)
        fib_a=np.array([6.5,5.5,4.5,3.5,2.5])
        nt=np.argmin(np.abs(fib_a/ifl-1.))
        nfib=np.int(fib_a[nt]+0.5)
        #print nfib
    else:
        nfib=7
        R50=-100
    f4.write(cubef+','+str(nfib)+','+str(R50)+','+str(e_rad)+','+str(red_0)+','+str(ab)+','+str(PA)+','+str(not_fit)+'\n')
    return nfib

def mock_sim_ill(idh,idhp,dir1='',fib_n=7,ho=0.704,Lam=0.7274,Om=0.2726,nl=110,fov=0.2,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'):
    dir1t=dir1.replace(' ','\ ')
    dirs=dir1t.split('/')
    dirs=filter(None,dirs)
    dirf=''
    dirf2=''
    for i in range(0, len(dirs)):
        dirf=dirf+'/'+dirs[i]
        dirf2=dirf2+'/'+dirs[i]
        if i > 2:
            sycall('mkdir -p '+dirf)
            #sycall('mkdir -p '+dirf)
    id=str(idh)
    outf='image_'+id
    cubef='ilust-'+str(idhp)+'-'+id+'.simu'
    dirf=dir1t+'photo_lin/'
    dirf2=dir1t+'Ensamble/'
    sycall('mkdir -p '+dirf)
    sycall('mkdir -p '+dirf2)
    sycall('mkdir -p '+dirf2+'Plots/')
    dir0=""
    dir1t=dir1t+id+'/'
    dir1=dir1+id+'/'
    dir1n=dir1.replace(" ","\ ")
    sycall('mkdir -p '+dir1t)
    stars = il.snapshot.loadSubhalo(basePath,135,idh,'stars')
    gas = il.snapshot.loadSubhalo(basePath,135,idh,'gas')
    dx_g = gas['Coordinates'][:,0]#
    dy_g = gas['Coordinates'][:,1]#
    dz_g = gas['Coordinates'][:,2]#
    vx_g = gas['Velocities'][:,0]#
    vy_g = gas['Velocities'][:,1]#
    vz_g = gas['Velocities'][:,2]#
    volm = gas['Volume'][:]
    dens = gas['Density'][:]#
    sfri = gas['StarFormationRate'][:]#
    meta_g=gas['GFM_Metallicity'][:]#     
    dx = stars['Coordinates'][:,0]#
    dy = stars['Coordinates'][:,1]#
    dz = stars['Coordinates'][:,2]#
    phot=stars['GFM_StellarPhotometrics'][:,5]
    mass=stars['Masses'][:]#
    mas0=stars['GFM_InitialMass'][:]#
    meta=stars['GFM_Metallicity'][:]#
    denS=stars['SubfindDensity']
    time=stars['GFM_StellarFormationTime'][:]#
    vx = stars['Velocities'][:,0]#
    vy = stars['Velocities'][:,1]#
    vz = stars['Velocities'][:,2]#
    mass_g=gas['Masses'][:]
    TE =   gas['InternalEnergy'][:]     
    #ages=stars['GFM_StellarAge'][:]
    temp_g=TE*((1e5)**2.0)*(0.938*1.7827e-24)*(4.0/(8.0-5.0*0.245))/(1.3807e-16)
    nt=np.where(time > 0)[0]
    dx=dx[nt]/ho
    dy=dy[nt]/ho
    dz=dz[nt]/ho
    vx=vx[nt]
    vy=vy[nt]
    vz=vz[nt]
    phot=phot[nt]
    mass_g=mass_g/ho*1e10
    mass=mass[nt]/ho*1e10
    mas0=mas0[nt]/ho*1e10
    meta=meta[nt]#/0.0127
    dx_g=dx_g/ho
    dy_g=dy_g/ho
    dz_g=dz_g/ho
    volm=volm/ho**3.0
    dens=dens*ho**2.0
    volm=float_((volm/(4.0*np.pi/3.0))**(1./3.0)*(3.08567758e19*100))
    dens=dens*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    Av_g=meta_g*(3.0*1.67262e-24*np.pi*dens*volm)/(4.0*np.log(10.)*3.0*5494e-8)
   # print meta_g,(10.0**(-0.59)*0.0127)
    nt1=np.where(meta_g > (10.0**(-0.59)*0.0127))[0]
    nt2=np.where(meta_g <= (10.0**(-0.59)*0.0127))[0]
    Av_g[nt1]=Av_g[nt1]/(10.0**(2.21-1.0))
    if len(nt2) > 0 :
        Av_g[nt2]=Av_g[nt2]/(10.0**(2.21-1.0)/(meta_g[nt2]/0.0127)**(3.1-1.0))
    #print np.amax(meta),np.amin(meta)
    zf=1/time[nt]-1
    xo=modes(dx,nit=7)
    yo=modes(dy,nit=7)
    zo=modes(dz,nit=7)
    x=dx-xo
    y=dy-yo
    z=dz-zo
    x_g=dx_g-xo
    y_g=dy_g-yo
    z_g=dz_g-zo
    x0=observer[0]-xo
    y0=observer[1]-yo
    z0=observer[2]-zo
    Rc=np.sqrt(x0**2.+y0**2.+z0**2.)
    red_0=reds_cos(Rc/1e3)
    Ra=Rc#/(1+red_0)
    A1=np.arctan2(y0,z0)
    A2=np.arcsin(x0/Ra)
    R1=np.array([[1,0,0],[0,np.cos(A1),-np.sin(A1)],[0,np.sin(A1),np.cos(A1)]])
    R2=np.array([[np.cos(A2),0,-np.sin(A2)],[0,1,0],[np.sin(A2),0,np.cos(A2)]])
    R3=np.array([[np.cos(A2),-np.sin(A2),0],[np.sin(A2),np.cos(A2),0],[0,0,1]])
    Ve=np.array([x,y,z])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x=Vf[0]
    y=Vf[1]
    z=Vf[2]
#    print modes(x,nit=7,n_s=0.8),modes(y,nit=7,n_s=0.8),modes(z,nit=7,n_s=0.8)
#    sys.exit()
    Ve=np.array([vx,vy,vz])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx=Vf[0]
    vy=Vf[1]
    vz=Vf[2]
    Ve=np.array([x_g,y_g,z_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    x_g=Vf[0]
    y_g=Vf[1]
    z_g=Vf[2]
    Ve=np.array([vx_g,vy_g,vz_g])
    Vf=np.dot(np.dot(R2,R1),Ve)
    vx_g=Vf[0]
    vy_g=Vf[1]
    vz_g=Vf[2]
    #print np.dot(np.dot(R2,R1),np.array([[x0],[y0],[z0]]))
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    print red_0
    print len(zf)
    age_s=cd.lookback_time(zf, **cosmo)/31557600./1e9
    [age_F,ind]=ages_definition(age_s,n_ages=55)
    if ptt.exists(dir1+cubef+'.fits.gz') == False or ptt.exists(dir1+cubef+'_L.fits.gz') == False:
        sim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
            ensamble_s(cubef+'.fits.gz',dir1=dir1)
            ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        ensamble_int(cubef+'.fits.gz',dir1=dir1)
        ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
    else:
       # sim_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,sig=sig,thet=thet,observer=observer)
        sycall('cp '+dir1n+cubef+'.fits.gz '+dirf)
        if ptt.exists(dir1+cubef.replace('simu','photo')+".r_Lc_rad.fits") == True:
            ensamble_s(cubef+'.fits.gz',dir1=dir1)
            ensamble_s(cubef+'_L.fits.gz',dir1=dir1,lig=1)
            sycall('cp '+dir1n+cubef+'_Ensemble.csv '+dirf2)
            sycall('cp '+dir1n+cubef+'_Ensemble_L.csv '+dirf2)
        ensamble_int(cubef+'.fits.gz',dir1=dir1)
        ensamble_int(cubef+'_L.fits.gz',dir1=dir1,lig=1)
        sycall('cp '+dir1n+cubef+'_Ensemble_int.csv '+dirf2)
        sycall('cp '+dir1n+cubef+'_Ensemble_int_L.csv '+dirf2)
        sycall('cp '+dir1n+'*.simu.fits.gz.pdf '+dirf2+'Plots/')
        
def mock_ill(fov=30.0,j=0,n_tot=1100,typef1="MaNGA",dir1='',file_out='mock_mass_ill_0.out',tempf="Temp.dat",basePath='/media/hjibarram/ADATA NH03/ILLUSTRIS/Illustris-1'): 
    #obs=[-106534,-106534,-106534]
    import random as rat
    obs=[ 52355.0, 43067.0, 43718.0]
    obs=[-106534.0,-106534.0,-106534.0]
    nl=110
    n_pix=440
    fov1=0.06
    sig=2.5
    thet=0.0
    plots=1
    rads=[0,0.5,1.0,2.0]
    f=h5py.File(il.snapshot.snapPath(basePath,135),'r')
    header = dict( f['Header'].attrs.items() )
    f.close()
    Om=header['Omega0']
    Lam=header['OmegaLambda']
    ho=0.704
    HaloID = il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloGrNr'])
    SubHaloMas= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloMassInRadType'])
    Len= il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloLenType'])
    Rad=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloHalfmassRad'])#SubhaloStellarPhotometricsRad'])
    XYZ=il.groupcat.loadSubhalos(basePath,135,fields=['SubhaloCM'])
    Lent=Len[:,0]
    x=XYZ[:,0]
    y=XYZ[:,1]
    z=XYZ[:,2]
    #Rad=Rad_T[:,4]
    Fin_mass=np.log10(SubHaloMas[:,4]/ho*1e10+1)
    n_sub=len(HaloID)
    cosmo = {'omega_M_0' : Om, 'omega_lambda_0' : Lam, 'h' : ho}
    cosmo = cd.set_omega_k_0(cosmo)
    dm=0.5
    Ms=11.5
    Mi=9.5
    nr=500
    nx=int((Ms-Mi)/dm)
    indf=[]
    if ptt.exists("Selection.dat") == False:
        f=open("Selection.dat","w")
        for i in range(0, nx):
            nt=np.where((Fin_mass > (Mi+dm*i)) & (Fin_mass <= (Mi+dm*(i+1))))[0]
            temp=nt[rat.sample(range(0, len(nt)), nr)]
            for j in range(0, nr):
                indf.extend([temp[j]])
                f.write(str(temp[j])+'\n')
    else:
        f=open("Selection.dat","r")
        for line in f:
            indf.extend([int(line.replace('\n',''))])
    indf=np.array(indf)
    f.close()
    #print indf
    
    #from matplotlib import cm
    #import matplotlib.pyplot as plt
    #fig = plt.figure(figsize=(6,5.5))
    #plt.scatter(obs[0]-x[indf]/ho,obs[1]-y[indf]/ho,c=Fin_mass[indf],vmin=8.5,vmax=11.5,s=10.5)
    #plt.plot([0],[0],'o')
    #plt.plot(obs[0]-x[indf[0]]/ho,obs[2]-z[indf[0]]/ho,'o',markersize=10,color='green') # 
    #plt.plot(obs[0]-66111.5123788,obs[2]-10145.8658918,'o',markersize=5,color='red')
    #plt.savefig(dev,dpi = 1000)
    #plt.show()
    #plt.close()
    
    #sys.exit()
    cont=1
    cont2=0
    #f3=open('Nonly_'+file_out,'w')
    if j == 0:
        fo=open(tempf,'w')
        f4=open('N'+file_out,'w')
        for i in range(0, 500):#len(indf)):
        #for i in range(0, n_sub):
            #if Fin_mass[i] > 9.5 and Fin_mass[i] < 11.5 and cont2 < 1100 and
            if  Lent[indf[i]] > 0:
                print indf[i], HaloID[indf[i]],Fin_mass[indf[i]],cont,cont2,Lent[indf[i]]
                name='ilust-'+str(HaloID[indf[i]])+'-'+str(indf[i])
                id=indf[i]
                id_p=HaloID[indf[i]]
                Rx=np.sqrt((obs[0]-x[id]/ho)**2.0+(obs[1]-y[id]/ho)**2.0+(obs[0]-z[id]/ho)**2.0)
                fov_p=Rad[id]*5./ho/Rx*180./np.pi*3600.
                fib_n=mock_photo_ill(id,id_p,f4,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,basePath=basePath)
                mock_sim_ill(id,id_p,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,basePath=basePath)
                mock_halo_ill(id,id_p,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,ifutype=typef1,basePath=basePath)
                fit3d_cen(name,fo,reens=1)
                cont=cont+1
                cont2=cont2+1
#                if i < n_sub-1:
#                    if HaloID[indf[i]] < HaloID[indf[i]+1]:
#                        cont=0
    else:
        fo=open(tempf+'temp','w')
        f4=open('N'+file_out+'temp','w')
        i=j
        if Fin_mass[i] > 9.5 and Fin_mass[i] < 11.5 and Lent[i] > 0: #cont2 < n_tot
            print i, HaloID[i],Fin_mass[i],cont,cont2
            name='ilust-'+str(HaloID[i])+'-'+str(i)
            id=i
            id_p=HaloID[i]
            Rx=np.sqrt((obs[0]-x[id]/ho)**2.0+(obs[1]-y[id]/ho)**2.0+(obs[0]-z[id]/ho)**2.0)
            fov_p=Rad[id]*5./ho/Rx*180./np.pi*3600.
            fib_n=mock_photo_ill(id,id_p,f4,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,basePath=basePath)
            mock_sim_ill(id,id_p,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,basePath=basePath)
            mock_halo_ill(id,id_p,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs,ifutype=typef1,basePath=basePath)
            fit3d_cen(name,fo,reens=0)
        else:
            print "NOT enough particles"
            print i, HaloID[i],Fin_mass[i],cont,cont2
    #f3.close()
    f4.close()
    fo.close()

<<<<<<< HEAD
def fit3d_only(name,dir_root='/home/hjibarram/FIT3D_py/new/',cat_root='spec',redo=0):
    dir1=dir_root+cat_root+'_lin_ana'
    n1='/'+name.split('-')[1]+'/'+name.split('-')[2]+'/coeffs_auto_ssp.'+name+'.int.out'
    if (ptt.exists(dir1+n1) == False) or (redo==1):
        call="./my_script.sh "+name
        sycall(call)
=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8

def fit3d(name,fo,fot_par,dir2='',dir_root='/home/hjibarram/FIT3D_py/new/',cat_root='spec',dir4='/home/hjibarram/FIT3D_py',redo=0,reens=0):
    import imp
    dir0=dir_root+cat_root+'_lin'
    dir1=dir_root+cat_root+'_lin_ana'
    dir3=dir_root+'Ensamble'
    dir_proc='/home/hjibarram/FIT3D_py/new/Pipe3D/mocks/proc_elines'
    #sifu.sycall("mkdir -p "+dir3)
    ens = imp.load_source('ensemble', '/home/hjibarram/FIT3D_py/soft_f/Mocks/sin_ensamble.py')
    n1='/'+name.split('-')[1]+'/'+name.split('-')[2]+'/coeffs_auto_ssp.'+name+'.int.out'
    #print ptt.exists(dir1+n1) == False,redo
    #sys.exit()
    [pa,eli,rad,nfit]=fot_par
    if (ptt.exists(dir1+n1) == False) or (redo==1):
        call="./my_script.sh "+name
        sycall(call)
        if np.abs(nfit) == 0:
            call=dir_proc+"/proc_elines.pl "+str(name)+" "+str(90.0+pa)+" "+str(eli)+" "+str(rad)+" .ps/CPS"
            sycall(call)
            ens.ensamble(name,dir1,dir3,dir4,fo,pdf=1,fits_f=1)
            file_e=dir3+'/'+name+'_Ensemble.csv'
            if ptt.exists(file_e) == True:
                age1,mgh1,sfh1,exy11,exy21 = ensam_mgh(file_e)
                plots_gp(age1,mgh1,sfh1,exy11,exy21,tit1=name)
                sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
            file_e=dir3+'/'+name+'_Ensemble_L.csv'
            if ptt.exists(file_e) == True:
                age1,mgh1,sfh1,exy11,exy21 = ensam_mgh(file_e)
                plots_gp(age1,mgh1,sfh1,exy11,exy21,tit1=name,lig=1)
                sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
    if reens == 1:
        if np.abs(nfit) == 0:
            call=dir_proc+"/proc_elines.pl "+str(name)+" "+str(90.0+pa)+" "+str(eli)+" "+str(rad)+" .ps/CPS"
            sycall(call)
            ens.ensamble(name,dir1,dir3,dir4,fo,pdf=1,fits_f=1)
            file_e=dir3+'/'+name+'_Ensemble.csv'
            if ptt.exists(file_e) == True:
                age1,mgh1,sfh1,exy11,exy21 = ensam_mgh(file_e)
                plots_gp(age1,mgh1,sfh1,exy11,exy21,tit1=name)
                sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
            file_e=dir3+'/'+name+'_Ensemble_L.csv'
            if ptt.exists(file_e) == True:
                age1,mgh1,sfh1,exy11,exy21 = ensam_mgh(file_e)
                plots_gp(age1,mgh1,sfh1,exy11,exy21,tit1=name,lig=1)
                sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')

def fit3d_cen(name,fo,dir_root='/home/hjibarram/FIT3D_py/new/',cat_root='spec',dir4='/home/hjibarram/FIT3D_py',redo=0,reens=0):
    import imp
    dir0=dir_root+cat_root+'_lin'
    dir1=dir_root+cat_root+'_lin_ana_cen'
    dir3=dir_root+'Ensamble'
    #sifu.sycall("mkdir -p "+dir3)
    ens = imp.load_source('ensemble', '/home/hjibarram/FIT3D_py/soft_f/Mocks/sin_ensamble.py')
    n1='/'+name.split('-')[1]+'/'+name.split('-')[2]+'/coeffs_auto_ssp.'+name+'.int.out'
    #print ptt.exists(dir1+n1) == False,redo
    #sys.exit()
    if (ptt.exists(dir1+n1) == False) or (redo==1):
        call="./my_script_cen.sh "+name
        sycall(call)
        ens.int_ensamble(name,dir1,dir3,dir4,fo,pdf=1,fits_f=0)
        file_e=dir3+'/'+name+'_Ensemble_int.csv'
        if ptt.exists(file_e) == True:
            age1,mgh1,sfh1,exy11,exy21 = ensam_mgh_int(file_e)
            plots_gp_int(age1,mgh1,sfh1,exy11,exy21,tit1=name)
            sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
        file_e=dir3+'/'+name+'_Ensemble_int_L.csv'
        if ptt.exists(file_e) == True:
            age1,mgh1,sfh1,exy11,exy21 = ensam_mgh_int(file_e)
            plots_gp_int(age1,mgh1,sfh1,exy11,exy21,tit1=name,lig=1)
            sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
    if reens == 1:
        ens.int_ensamble(name,dir1,dir3,dir4,fo,pdf=1,fits_f=0)
        file_e=dir3+'/'+name+'_Ensemble_int.csv'
        if ptt.exists(file_e) == True:
            age1,mgh1,sfh1,exy11,exy21 = ensam_mgh_int(file_e)
            plots_gp_int(age1,mgh1,sfh1,exy11,exy21,tit1=name)
            sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
        file_e=dir3+'/'+name+'_Ensemble_int_L.csv'
        if ptt.exists(file_e) == True:
            age1,mgh1,sfh1,exy11,exy21 = ensam_mgh_int(file_e)
            plots_gp_int(age1,mgh1,sfh1,exy11,exy21,tit1=name,lig=1)
            sycall('mv *'+name+'*.pdf '+dir3+'/Plots/')
            
def mask_rad(name):
    dir1=name.split("-")[1]
    dir0='/home/hjibarram/FIT3D_py/new/no_name/'
    names=name.split("-")[2]
    filer=dir0+dir1+"/"+names+"/"+name+".p_e.pdl_r.fits"
    filem=dir0+dir1+"/"+names+"/"+name+".p_e.pdl_r_mask.fits"
    if ptt.exists(filer) == False:
        print "The File does not exist"
        sys.exit()
    else:
        [pdl_rad, hdr]=gdata(filer,0, header=True)
        mask=pdl_rad
        mask[np.where(pdl_rad<=1.5)]=1.0
        mask[np.where(pdl_rad>1.5)]=0.0
        wfits(filem,mask,hdr)
    return filem
                    
def fit3d_cen_stat(name,fo,dir_root='/home/hjibarram/FIT3D_py/new/',cat_root='spec',dir4='/home/hjibarram/FIT3D_py',redo=0,n_m=500,n_mo=0):
    file_m=mask_rad(name)
    import imp
    dir0=dir_root+cat_root+'_lin'
    dir1=dir_root+cat_root+'_lin_ana_cen'
    dir3=dir_root+'Ensamble'
    ens = imp.load_source('ensemble', '/home/hjibarram/FIT3D_py/soft_f/Mocks/sin_ensamble.py')
    for i in range(n_mo, n_m):
        call="./my_script_cen_m.sh "+name+" "+str(i)+" "+file_m
        sycall(call)
        ens.int_ensamble(name,dir1,dir3,dir4,fo,pdf=1,fits_f=0,m_t='_'+str(i))
        
def ensam_mgh(file):
    n_age=66#39
    nrb=1+4*3+1
    ensamb=np.zeros([n_age,nrb])
    co=0
    f=open(file,"r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1
    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1:4])
    sfh=ensamb[:,7:10]
    emgh=(ensamb[:,10:13])
    dth=ensamb[:,13]/1e3
    Ddth=np.zeros(len(dth))
    [nx,ny]=mgh.shape
    Dmgh=np.zeros([nx,ny])
    exy1=np.zeros([nx,ny])
    exy2=np.zeros([nx,ny])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i,:]-mgh[i+1,:])
            Ddth[i]=np.abs(dth[i]-dth[i+1])
        elif i == len(Ddth)-1:
            Dmgh[i,:]=np.abs(mgh[i-1,:]-mgh[i,:])
            Ddth[i]=np.abs(dth[i-1]-dth[i])
    #exy1[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)

    for i in range(0, 3):
        val1=np.abs(Dmgh[:,i]/(dth+Ddth)*dth)+emgh[:,i]
        val2=np.abs(Dmgh[:,i]/(dth-Ddth)*dth)+emgh[:,i]
        for j in range(0, len(val1)):
            if val1[j] >= 0.15:
                val1[j]=0.1
            if val2[j] >= 0.15:
                val2[j]=0.1
        exy1[:,i]=val1
        exy2[:,i]=val2
    return age,mgh,sfh,exy1,exy2

def plots_gp(age1,mgh1,sfh1,exy11,exy21,dir3="",tit1="MOCK",labe="",lig=0):
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    tit2=""
    tit3=""
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    if lig == 0:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title(tit1+" MGH "+tit2)
    else:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title(tit1+" LGH "+tit2)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    plt.semilogx(age1,mgh1[:,0],'',color="b",label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=1.5)
    plt.semilogx(age1,mgh1[:,1],'--',color="g",label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=1.5)
    plt.semilogx(age1,mgh1[:,2],':',color="r",label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=1.5)
    ax.fill_between(age1,mgh1[:,0]+exy11[:,0],mgh1[:,0]-exy21[:,0],alpha='0.20',color="b")
    ax.fill_between(age1,mgh1[:,1]+exy11[:,1],mgh1[:,1]-exy21[:,1],alpha='0.20',color="g")
    ax.fill_between(age1,mgh1[:,2]+exy11[:,2],mgh1[:,2]-exy21[:,2],alpha='0.20',color="r")
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 0:
        plt.savefig(dir3+'MGH_'+tit1+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir3+'LGH_'+tit1+'.pdf',dpi = 1000)
    plt.close()

    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 0:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title(tit1+" SFH "+tit2)
    else:
        plt.ylabel(r'$SFH_L(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title(tit1+" SFHL "+tit2)
    plt.semilogx(age1,sfh1[:,0],'',color="b",drawstyle='steps',label='$'+('%6.1f' % 0.0)+'R_{50}<R<'+('%6.1f' % 0.5)+'R_{50}$',lw=2)
    plt.semilogx(age1,sfh1[:,1],'--',color="g",drawstyle='steps',label='$'+('%6.1f' % 0.5)+'R_{50}<R<'+('%6.1f' % 1.0)+'R_{50}$',lw=2)
    plt.semilogx(age1,sfh1[:,2],':',color="r",drawstyle='steps',label='$'+('%6.1f' % 1.0)+'R_{50}<R<'+('%6.1f' % 1.5)+'R_{50}$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 0:
        plt.savefig(dir3+'SFH_'+tit1+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir3+'SFHL_'+tit1+'.pdf',dpi = 1000)
    plt.close()
    
def ensam_mgh_int(file):
    n_age=39
    nrb=6
    ensamb=np.zeros([n_age,nrb])
    co=0
    f=open(file,"r")
    for line in f:
        if not "#" in line:
            data=line.split(";")
            data=filter(None,data)
            for k in range(0,nrb):
                ensamb[co,k]=np.abs(float_(data[k]))
            co=co+1 
    age=10**(ensamb[:,0]-9)
    mgh=np.log10(ensamb[:,1])
    sfh=ensamb[:,3]
    emgh=(ensamb[:,4])
    dth=ensamb[:,5]/1e3
    Ddth=np.zeros(len(dth))
    nx=len(mgh)
    Dmgh=np.zeros([nx])
    exy1=np.zeros([nx])
    exy2=np.zeros([nx])
    for i in range(0, len(Ddth)):
        if i < len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i]-mgh[i+1])
            Ddth[i]=np.abs(dth[i]-dth[i+1])
        elif i == len(Ddth)-1:
            Dmgh[i]=np.abs(mgh[i-1]-mgh[i])
            Ddth[i]=np.abs(dth[i-1]-dth[i])
    #exy1[:,0]=np.amin([np.abs(Dmgh[:,0]/(dth-Ddth)*dth)+emgh[:,0],np.abs(Dmgh[:,0]/(dth+Ddth)*dth)+emgh[:,0]],axis=0)

    val1=np.abs(Dmgh/(dth+Ddth)*dth)+emgh
    val2=np.abs(Dmgh/(dth-Ddth)*dth)+emgh
    for j in range(0, len(val1)):
        if val1[j] >= 0.15:
            val1[j]=0.1
        if val2[j] >= 0.15:
            val2[j]=0.1
        if mgh[j]+val1[j] > 1.0:
            val1[j]=1.0-mgh[j]
    exy1=val1
    exy2=val2
    return age,mgh,sfh,exy1,exy2

def plots_gp_int(age1,mgh1,sfh1,exy11,exy21,dir3="",tit1="MOCK",labe="",lig=0):
    matplotlib.use('agg')
    import matplotlib.pyplot as plt
    tit2=""
    tit3=""
    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.ylim(0.4,1.09)
    if lig == 0:
        plt.ylabel(r'$\mathcal{M}(t)/\mathcal{M}_{0}$',fontsize=16)
        plt.title(tit1+" MGH "+tit2)
    else:
        plt.ylabel(r'$\mathcal{L}(t)/\mathcal{L}_{0}$',fontsize=16)
        plt.title(tit1+" LGH "+tit2)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    plt.semilogx(age1,mgh1,'',color="b",label='$'+('%6.1f' % 0.0)+'R=15arcsec$',lw=1.5)
    ax.fill_between(age1,mgh1+exy11,mgh1-exy21,alpha='0.20',color="b")
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.90,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.70,'--',color='black')
    plt.semilogx(np.arange(-20,20,.1),np.ones(400)*0.50,'--',color='black')
    plt.legend(loc=3)
    fig.tight_layout()
    if lig == 0:
        plt.savefig(dir3+'MGH_int_'+tit1+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir3+'LGH_int_'+tit1+'.pdf',dpi = 1000)
    plt.close()

    fig, ax = plt.subplots(figsize=(6,5.5))
    plt.xlim(0.1,18)
    plt.xlabel(r'$Look Back Time\ [Gyr]$',fontsize=16)
    if lig == 0:
        plt.ylabel(r'$SFH(t)\ [\mathcal{M}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title(tit1+" SFH "+tit2)
    else:
        plt.ylabel(r'$SFH_L(t)\ [\mathcal{L}_{\odot}yr^{-1}]$',fontsize=16)
        plt.title(tit1+" SFHL "+tit2)
    plt.semilogx(age1,sfh1,'',color="b",drawstyle='steps',label='$'+('%6.1f' % 0.0)+'R=15arcsec$',lw=2)
    plt.legend(loc=2)
    fig.tight_layout()
    if lig == 0:
        plt.savefig(dir3+'SFH_int_'+tit1+'.pdf',dpi = 1000)
    else:
        plt.savefig(dir3+'SFHL_int_'+tit1+'.pdf',dpi = 1000)
    plt.close()
<<<<<<< HEAD
    
def copy_data(mjd,obs,field):
    call='mkdir -p $BOSS_SPECTRO_DATA/'+mjd+'/'
    sycall(call)
    call='mkdir -p $SDSSCORE/'+mjd+'/'
    sycall(call)
    call='mkdir -p $SDSSCORE/photofield/'+field
    sycall(call)
    call='cp '+obs+'/bhm/'+mjd+'/raw_mock/*.gz $BOSS_SPECTRO_DATA/'+mjd+'/'
    sycall(call)
    call='cp '+obs+'/bhm/'+mjd+'/raw_mock/*.par $SDSSCORE/'+mjd+'/'
    sycall(call)
    call='cp '+obs+'/bhm/'+mjd+'/raw_mock/*'+field+'.fits $SDSSCORE/photofield/'+field+'/'
    sycall(call)
    
def run_spec2d(mjd='56000',field=0,pbs=0):
    field_id=id_str(field,n_z=5)
    if pbs == 1:
        call='./scritp.sh '+mjd
    else:
        call='./scritp_2.sh '+mjd+' '+field_id
    sycall(call)
    
    
=======
>>>>>>> 043939233f85f60c203ef4e87794d197c76e7dc8
