import matplotlib
import sys
import numpy as np
import cosmolopy.distance as cd
from numpy import float_
from numpy import absolute as abs
from numpy import random as ran
import time
import os.path as ptt
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
                Rt=np.sqrt( ((xi+xf)/2.0-x)**2+((yi+yf)/2.0-y)**2. )
                if Rt <= dr/2.0 :
                    Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0+1*dt**2.0
            else:
                Ft=Ft+np.exp(-((((xi+xf)/2.0-xo)**2+((yi+yf)/2.0-yo)**2.)/sig**2.)/2.0)/(2*np.pi*sig**2.0)*dt**2.0+1*dt**2.0
    print x,y
    return Ft
        
        
                

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
    
def fib_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="SDSS"):
    nh=dens
    fact=nh/10.0
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
        if psfi <= 0:
            seeing=1.43
        else:
            seeing=psfi
    elif "CALIFA" in ifutype:
        pix_s=1.0#arcsec
        scp_s=56.02#microns per arcsec
        fibA=197.4
        fibB=150.0
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
    pht_g =asosiate_pho(ssp_template,wave,age_ssp,met_ssp,ml_ssp,mass_gssp,met_g,Rs,nh)
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
        spect=interp1d(wave,spect_i,bounds_error=False,fill_value=0.)(wave_f)
        spect[np.isnan(spect)]=0
        spec_val[0]=Av_s/len(nt)
    sycall('echo '+str(len(nt))+'  GAS')
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
        spect_gf=interp1d(wave_g,spect_g_i,bounds_error=False,fill_value=0.)(wave_f)
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
    
   


def cube_conv(outf,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,met_s,mass_s,met_g,vol,dens,sfri,temp_g,Av_g,mass_g,template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",dir_o='',Flux_m=20.0,psfi=0,SNi=15.0,red_0=0.01,ho=0.704,Lam=0.7274,Om=0.2726,nl=7,fov=30.0,sig=2.5,thet=0.0,pdf=2,rx=[0,0.5,1.0,2.0],ifutype="MaNGA"):
    #if not "MaNGA" in ifutype or not "CALIFA" in ifutype or not "MUSE" in ifutype:
    #    ifutype="MaNGA"
    nh=dens#*1e10/(3.08567758e19*100)**3.0*1.9891e30/1.67262178e-27
    fact=nh/10.0
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
    cdelt_w=1.25
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
    spec_val=np.zeros([30,ndt*ns])
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
                    Av_s=Av+Av_s
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
                    met_ligt=10.0**(met_ligt/Lt)
                    age_ligt=10.0**(age_ligt/Lt)
                    Av_ligt=(Av_ligt/Lt)
                    Av_flux=(Av_flux/Ft)
                    Ve_ligt=(Ve_ligt/Lt)
                    Ve_flux=(Ve_flux/Ft)
                    met_mas=10.0**(met_mas/mass_t)
                    age_mas=10.0**(age_mas/mass_t)
                    va_1=np.array(va_1)
                    wf_t=np.array(wf_t)
                    wl_t=np.array(wl_t)
                    Sig_flux=np.sqrt(np.nansum(np.abs(wf_t)*(Ve_flux-va_1)**2.0)/(np.nansum(np.abs(wf_t))-np.nansum(wf_t**2.0)/np.nansum(np.abs(wf_t))))
                    Sig_ligt=np.sqrt(np.nansum(np.abs(wl_t)*(Ve_ligt-va_1)**2.0)/(np.nansum(np.abs(wl_t))-np.nansum(wl_t**2.0)/np.nansum(np.abs(wl_t))))
                spect_t[np.isnan(spect_t)]=0
                spect_i=inst_disp(wave,spect_t,sigma_inst)
                spect=interp1d(wave,spect_i,bounds_error=False,fill_value=0.)(wave_f)
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
                                spect_sfg=spect_sfg+ran.randn(nw_g)*np.median(spect_sfg)*0.01
                                spect_g=spect_sfg+spect_g    
                                va_1g.extend([v_rad_g[nt[k]]])
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
                    va_1g=np.array(va_1g)
                    wf_tg=np.array(wf_tg)
                    wl_tg=np.array(wl_tg)
                    Sig_flux_g=np.sqrt(np.nansum(np.abs(wf_tg)*(Veg_flux-va_1g)**2.0)/(np.nansum(np.abs(wf_tg))-np.nansum(wf_tg**2.0)/np.nansum(np.abs(wf_tg))))
                    Sig_ligt_g=np.sqrt(np.nansum(np.abs(wl_tg)*(Veg_ligt-va_1g)**2.0)/(np.nansum(np.abs(wl_tg))-np.nansum(wl_tg**2.0)/np.nansum(np.abs(wl_tg))))    
                spec_val[6,con]=np.sum(Av_g[nt_g])         
                spec_val[4,con]=Av_sg/len(nt_g)
                spect_g[np.isnan(spect_g)]=0
                spect_g_i=inst_disp(wave_g,spect_g,sigma_inst)
                spect_gf=interp1d(wave_g,spect_g_i,bounds_error=False,fill_value=0.)(wave_f)
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
            spec_val[5,con]=spec_val[0,con]+spec_val[4,con]
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
    ifu_v=np.zeros([30,nl,nl])
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
            spt_val=np.zeros(30)
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
    ifu_v[1,:,:]=np.log10(ifu_v[1,:,:]+1.0)
    ifu_v[10,:,:]=np.log10(ifu_v[10,:,:]+1.0)
    ifu_v[29,:,:]=np.log10(ifu_v[29,:,:]+1.0)
    ifu_v[15,:,:]=-2.5*np.log10(ifu_v[15,:,:]+0.0001)
    ifu_v[17,:,:]=-2.5*np.log10(ifu_v[17,:,:]+0.0001)
    h1t=pyf.PrimaryHDU(ifu_v)
    h=h1t.header
    h["NAXIS"]=3
    h["NAXIS3"]=30
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
    h['Type16']=('FLUX    ','1e-16 ergs/s/cm2')
    h['Type17']=('Av_fw   ','Mag')
    h['Type18']=('VEL_lw  ','km/s')
    h['Type19']=('VEL_fw  ','km/s')
    h['Type20']=('DIS_lw   ','km/s BETA')
    h['Type21']=('DIS_fw   ','km/s BETA')
    h['Type22']=('Av_lw_g  ','Mag')
    h['Type23']=('Av_fw_g  ','Mag')
    h['Type24']=('VEL_mw   ','km/s')
    h['Type25']=('VEL_lw   ','km/s')
    h['Type26']=('DIS_l_gas','km/s BETA')
    h['Type27']=('DIS_f_gas','km/s BETA')
    h['Type28']=('FLUX_gas','1e-16 ergs/s/cm2 Bolometric')
    h['Type29']=('LUM_gas ','log10(Lsun) Bolometric')
    h['RADECSYS']='ICRS    '
    h['SYSTEM']='FK5     '
    h['EQUINOX']=2000.00
    h['PSF']=sig
    h['FOV']=Rifu*2.0
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
    h['PSF']=sig
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
    the=the*180/np.pi*3600+ran.randn(len(rad))*1.43#/2.0
    phi=phi*180/np.pi*3600+ran.randn(len(rad))*1.43#/2.0
    phi_g=np.arcsin(x_g/radA_g)
    the_g=np.arcsin(y_g/(radA_g*np.cos(phi_g)))
    the_g=the_g*180/np.pi*3600+ran.randn(len(rad_g))*1.43#/2.0
    phi_g=phi_g*180/np.pi*3600+ran.randn(len(rad_g))*1.43#/2.0
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
                 
def mock_halo(idh,basename='artsp8-',template3="../home/sanchez/ppak/legacy/gsd61_156.fits",template5="../../Base_bc03/templete_bc03_5.fits",template2="templete_gas.fits",file=file,file2='',dir1='',fib_n=7,psf=0,ho=0.704,Lam=0.7274,SN=15.0,Fluxm=20.0,Om=0.2726,nl=110,fov=30.0,fov1=0.2,sig=2.5,thet=0.0,plots=1,rx=[0,0.5,1.0,2.0],observer=[0,0,0],ifutype="MaNGA"):
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
    mass_g=mass_g/ho
    mass=mass[nt]/ho
    mas0=mas0[nt]/ho
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
        cube_conv(cubef,x,y,z,vx,vy,vz,x_g,y_g,z_g,vx_g,vy_g,vz_g,age_s,meta,mass,meta_g,volm,dens,sfri,temp_g,Av_g,mass_g,template3=template3,template5=template5,template2=template2,psfi=psf,SNi=SN,Flux_m=Fluxm,dir_o=dir1,red_0=red_0,ho=ho,Lam=Lam,Om=Om,nl=fib_n,fov=fov,sig=sig,thet=thet,ifutype=ifutype)
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
    mass_g=mass_g/ho
    mass=mass[nt]/ho
    mas0=mas0[nt]/ho
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
        mass_g=mass_g/ho
        mass=mass[nt]/ho
        mas0=mas0[nt]/ho
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
        mass_g=mass_g/ho
        mass=mass[nt]/ho
        mas0=mas0[nt]/ho
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
        mass_g=mass_g/ho
        mass=mass[nt]/ho
        mas0=mas0[nt]/ho
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
        mass_g=mass_g/ho
        mass=mass[nt]/ho
        mas0=mas0[nt]/ho
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
    


def mock_sp(fib_n,ang,modt=0,template_1="libs/gsd61_156.fits",template_2="libs/templete_gas.fits",template_3="libs/templete_bc03_5.fits",template="libs/templete_bc03_2.fits",n_pix=440,fov_p=0,fov=0,rx=142135.5,Om=0.2726,Lam=0.7274,ho=0.704,cam=0,vx=[-0.7156755,-0.5130859,0.4738687],vy=[0.6984330,-0.5257526,0.4855672],vz=[0.0000000,0.6784741,0.7346244],base_name='artsp8-',typef1="MaNGA",id1='A2-0',psf=0,redo=0,SN=15.0,Fluxm=20.0,dir1='',file_red="sp8/star_out_0.dat",file_gas="sp8/Gas_out_0.dat",file_out='mock_mass_ill_0.out',file_out_f='mock_mass_ill.out'):
    sig=2.5
    thet=0.0
    plots=1
    nl=110
    if cam == 0:
        cam=(rx*2.0)*7.0/np.float_(fib_n)#/2.34 # 2.34 sp6#/3.0#5.5#/2.2#IF CALIFA 5.5
    cam=-cam
    if fov_p == 0:
        fov_p=np.round(rx/np.abs(cam)*91.0)
    if fov == 0:
        fov=np.round(rx/np.abs(cam)*62.0)#30
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
        mock_halo(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    if modt == 1:
        mock_photo(id1,basename=base_name,file=file_red,template2=template_2,template=template,file2=file_gas,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=n_pix,fov=fov_p,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1)
        mock_halo(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
        mock_halo_s(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype='SDSS')
    if modt == 2:
        mock_halo(id1,basename=base_name,file=file_red,template3=template_1,template5=template_3,template2=template_2,file2=file_gas,SN=SN,psf=psf,Fluxm=Fluxm,dir1=dir1,fib_n=fib_n,ho=ho,Lam=Lam,Om=Om,nl=nl,fov=fov,fov1=fov1,sig=sig,thet=thet,plots=plots,rx=rads,observer=obs_ang1,ifutype=typef1)
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
