#!/usr/bin/python
import sys
import numpy as np
import sin_ifu_mod as sifu
import warnings
warnings.filterwarnings("ignore")
sys.argv=filter(None,sys.argv)
if len(sys.argv) >= 2:
    var1=sys.argv[1].split(';')
    var1=filter(None,var1)
    if len(var1) == 1:
        if "Help" in var1[0] or "help" in var1[0]:
            sifu.sycall('echo USAGE:')
            sifu.sycall("ifu_mocks 'id_name;base_name;dir_outputs;mod;type_ifu' 'file_sim_stars;file_sim_gas' 'fib_n;angle;SN;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
            sifu.sycall("id_name = name of the mock")
            sifu.sycall("base_name = base name of the mock: base_name-id_name")
            sifu.sycall("dir_outputs = output directory")
            sifu.sycall("mod = Mode of the mocks: mod=0 implies full set of mocks, mod=1 implies photometry IFU and single fiber only, mod=2 implies only IFU mocks, mod=3 implies only SDSS single fiber mocks")
            sifu.sycall("type_ifu = Type of the IFU mock: 'MaNGA' for MaNGA type, 'CALIFA' for CALIFA type and 'MUSE' for muse type")
            sifu.sycall("file_sim_stars = File for the stellar particules of the desired hidrodynamical simulation")
            sifu.sycall("file_sim_gas = File for the gas cells of the desired hidrodynamical simulation")
            sifu.sycall("fib_n = number of fibers per octagon size")
            sifu.sycall("angle = angle between the reference vector e1 and the observer")
            sifu.sycall("SN = Signal to noise ratio")
            sys.exit()
            
            #;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
        else:
            id1=var1[0]
            basename='artsp8-'
            dir1='../'
            modt=1
            typef1="MaNGA"
    if len(var1) == 2:
        id1=var1[0]
        basename=var1[1]
        dir1='../'
        modt=1
        typef1="MaNGA"
    if len(var1) == 3:
        id1=var1[0]
        basename=var1[1]
        dir1=var1[2]
        modt=1
        typef1="MaNGA"
    if len(var1) == 4:
        id1=var1[0]
        basename=var1[1]
        dir1=var1[2]
        modt=np.int(var1[3])
        typef1="MaNGA"
    if len(var1) == 5:
        id1=var1[0]
        basename=var1[1]
        dir1=var1[2]
        modt=np.int(var1[3])
        typef1=var1[4]
else:
    basename='artsp8-'
    typef1="MaNGA"
    id1='A-0'
    dir1='../'
    modt=1
if len(sys.argv) >= 3:
    var2=sys.argv[2].split(';')
    var2=filter(None,var2)
    if len(var2) == 1:
        if 'none' in var2[0]:
            file_red="sp8/star_out_0.dat"
            file_gas="sp8/Gas_out_0.dat"
        else:
            file_red=var2[0]
            file_gas="sp8/Gas_out_0.dat"
    if len(var2) == 2:
        file_red=var2[0]
        file_gas=var2[1]
else:
    file_red="sp8/star_out_0.dat"
    file_gas="sp8/Gas_out_0.dat"
if len(sys.argv) >= 4:
    var3=sys.argv[3].split(';')
    var3=filter(None,var3)
    if len(var3) == 1:
        if 'none' in var3[0]:
            ang=0.0
            fib_n=7
            SN=15.0
            Fluxm=26.83
            psf=1.43
        else:
            fib_n=np.int(var3[0])
            ang=0.0
            SN=15.0
            Fluxm=26.83
            psf=1.43
    if len(var3) == 2:
        fib_n=np.int(var3[0])
        ang=np.float(var3[1])
        SN=15.0
        Fluxm=26.83
        psf=1.43
    if len(var3) == 3:
        fib_n=np.int(var3[0])
        ang=np.float(var3[1])
        SN=np.float(var3[2])
        Fluxm=26.83
        psf=1.43
    if len(var3) == 4:
        fib_n=np.int(var3[0])
        ang=np.float(var3[1])
        SN=np.float(var3[2])
        Fluxm=np.float(var3[3])
        psf=1.43
    if len(var3) == 5:   
        fib_n=np.int(var3[0])
        ang=np.float(var3[1])
        SN=np.float(var3[2])
        Fluxm=np.float(var3[3])
        psf=np.float(var3[4])
else:
    ang=0.0
    fib_n=7
    SN=15.0
    Fluxm=26.83
    psf=1.43
if len(sys.argv) >= 5:
    var4=sys.argv[4].split(';')
    var4=filter(None,var4)
    if len(var4) == 1:
        if "none" in var4[0]:
            fov_p=0
            fov=0
            cam=0
            n_pix=440
            rx=142135.5
        else:
            n_pix=np.int(var4[0])
            fov_p=0
            fov=0
            cam=0
            rx=142135.5
    if len(var4) == 2:
        n_pix=np.int(var4[0])
        rx=np.float(var4[1])
        fov_p=0
        fov=0
        cam=0
    if len(var4) == 3:
        n_pix=np.int(var4[0])
        rx=np.float(var4[1])
        fov_p=np.float(var4[2])
        fov=0
        cam=0
    if len(var4) == 4:
        n_pix=np.int(var4[0])
        rx=np.float(var4[1])
        fov_p=np.float(var4[2])
        fov=np.float(var4[3])
        cam=0
    if len(var4) == 5:
        n_pix=np.int(var4[0])
        rx=np.float(var4[1])
        fov_p=np.float(var4[2])
        fov=np.float(var4[3])
        cam=np.float(var4[4])
else:
    fov_p=0
    fov=0
    cam=0
    n_pix=440
    rx=142135.5
if len(sys.argv) >= 6:
    var5=sys.argv[5].split(';')
    var5=filter(None,var5)
    if len(var5) == 1:
        if "none" in var5[0]:
            vx=[-0.7156755,-0.5130859,0.4738687]
            vy=[0.6984330,-0.5257526,0.4855672]
            vz=[0.0000000,0.6784741,0.7346244]
        else:
            vx1=var5[0].split(',')
            vx=[]
            if len (vx1) == 3:
                for i in range(0,3):
                    vx.extend([np.float(vx1[i])])
            else:
                vx=[-0.7156755,-0.5130859,0.4738687]
            vy=[1,1,-(vx[0]+vx[1])/vx[2]]/np.sqrt(2.0+((vx[0]+vx[1])/vx[2])**2.0)
            vz=np.cross(vx,vy)
            vz=vz/np.sqrt(vz[0]**2.0+vz[1]**2.0+vz[2]**2.0)
    if len(var5) == 2:
        vx1=var5[0].split(',')
        vy1=var5[1].split(',')
        vx=[]
        vz=[]
        if len (vx1) == 3:
            for i in range(0,3):
                vx.extend([np.float(vx1[i])])
        else:
            vx=[-0.7156755,-0.5130859,0.4738687]
        if len (vy1) == 3:
            for i in range(0,3):
                vz.extend([np.float(vy1[i])])
        else:
            vz=[0.0000000,0.6784741,0.7346244]
        vy=np.cross(vx,vz)
        vy=vy/np.sqrt(vy[0]**2.0+vy[1]**2.0+vy[2]**2.0)
    if len(var5) == 3:
        vx1=var5[0].split(',')
        vy1=var5[1].split(',')
        vz1=var5[2].split(',')
        vx=[]
        vy=[]
        vz=[]
        if len (vx1) == 3:
            for i in range(0,3):
                vx.extend([np.float(vx1[i])])
        else:
            vx=[-0.7156755,-0.5130859,0.4738687]
        if len (vy1) == 3:
            for i in range(0,3):
                vy.extend([np.float(vy1[i])])
        else:
            vy=[0.6984330,-0.5257526,0.4855672]
        if len (vz1) == 3:
            for i in range(0,3):
                vz.extend([np.float(vz1[i])])
        else:
            vz=[0.0000000,0.6784741,0.7346244]
else:
    vx=[-0.7156755,-0.5130859,0.4738687]
    vy=[0.6984330,-0.5257526,0.4855672]
    vz=[0.0000000,0.6784741,0.7346244]
if len(sys.argv) >= 7:
    var6=sys.argv[6].split(';')
    var6=filter(None,var6)
    if len(var6) == 1:
        if "none" in var6[0]:
            template_0="libs/templete_bc03_2.fits"
            template_1="libs/gsd61_156.fits"
            template_2="libs/templete_gas.fits"
            template_3="libs/templete_bc03_5.fits"
        else:
            template_0=var6[0]        
            template_1="libs/gsd61_156.fits"
            template_2="libs/templete_gas.fits"
            template_3="libs/templete_bc03_5.fits"
    if len(var6) == 2:
        template_0=var6[0]
        template_1=var6[1]
        template_2="libs/templete_gas.fits"
        template_3="libs/templete_bc03_5.fits" 
    if len(var6) == 3:
        template_0=var6[0]
        template_1=var6[1]
        template_2=var6[2] 
        template_3="libs/templete_bc03_5.fits" 
    if len(var6) == 4:
        template_0=var6[0]
        template_1=var6[1]
        template_2=var6[2] 
        template_3=var6[3] 
else:
    template_0="libs/templete_bc03_2.fits"
    template_1="libs/gsd61_156.fits"
    template_2="libs/templete_gas.fits"
    template_3="libs/templete_bc03_5.fits"
if len(sys.argv) >= 8:
    var7=sys.argv[7].split(';')
    var7=filter(None,var7)
    if len(var7) == 1:
        if "none" in var7[0]:
            Om=0.2726
            Lam=0.7274
            ho=0.704
        else:
            Om=np.float(var7[0])
            Lam=1.0-Om
            ho=0.704
    if len(var7) == 2:
        Om=np.float(var7[0])
        Lam=np.float(var7[1])
        ho=0.704
    if len(var7) == 3:
        Om=np.float(var7[0])
        Lam=np.float(var7[1])
        ho=np.float(var7[2])
else:
    Om=0.2726
    Lam=0.7274
    ho=0.704

if len(sys.argv) == 9:
    redo=np.float(sys.argv[8])
else:
    redo=0
file_out='mock_mass_ill_0.out'
file_out_f='mock_mass_ill.out'
sifu.mock_sp(fib_n,ang,modt=modt,template_1=template_1,template_2=template_2,template_3=template_3,template=template_0,n_pix=n_pix,fov_p=fov_p,fov=fov,rx=rx,Om=Om,Lam=Lam,ho=ho,cam=cam,vx=vx,vy=vy,vz=vz,base_name=basename,typef1=typef1,id1=id1,psf=psf,redo=redo,SN=SN,Fluxm=Fluxm,dir1=dir1,file_red=file_red,file_gas=file_gas,file_out=file_out,file_out_f=file_out_f)