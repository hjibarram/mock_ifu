#!/usr/bin/python
import sys
import numpy as np
import sin_ifu_mod as sifu
import warnings
import os.path as pt
warnings.filterwarnings("ignore")
sys.argv=filter(None,sys.argv)
def_t=0
conf_t=0
if len(sys.argv) == 1:
    sifu.sycall('echo USAGE:')
    sifu.sycall("echo Running Using Docker:")
    sifu.sycall("echo docker pull hjibarram/ifu_mocks")
    sifu.sycall("echo docker run -p 8000:8000 -v '`pwd`'/dir_outputs:/dir_outputs  hjibarram/ifu_mocks 'id_name;base_name;dir_outputs;mod;type_ifu' 'file_sim_stars;file_sim_gas' 'fib_n;angle;SN;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
    sifu.sycall("echo Running Using Python:")
    sifu.sycall("echo python mocks.py 'id_name;base_name;dir_outputs;mod;type_ifu' 'file_sim_stars;file_sim_gas' 'fib_n;angle;SN;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
    sifu.sycall("echo  ")
    sifu.sycall("echo key words:")
    sifu.sycall("echo id_name = name of the mock")
    sifu.sycall("echo base_name = base name of the mock: base_name-id_name")
    sifu.sycall("echo dir_outputs = output directory")
    sifu.sycall("echo mod = Mode of the mocks: mod=0 implies full set of mocks, mod=1 implies photometry IFU and single fiber only, mod=2 implies only IFU mocks, mod=3 implies only SDSS single fiber mocks")
    sifu.sycall("echo type_ifu = Type of the IFU mock: 'MaNGA' for MaNGA type, 'CALIFA' for CALIFA type and 'MUSE' for muse type")
    sifu.sycall("echo file_sim_stars = File for the stellar particules of the desired hidrodynamical simulation")
    sifu.sycall("echo file_sim_gas = File for the gas cells of the desired hidrodynamical simulation")
    sifu.sycall("echo fib_n = number of fibers per octagon size")
    sifu.sycall("echo angle = angle between the reference vector e1 and the observer")
    sifu.sycall("echo SN = Signal to noise ratio")
    sifu.sycall("echo  ")
    sifu.sycall("echo Running with Configuration file:")
    sifu.sycall("echo docker run -p 8000:8000 -v '`pwd`'/dir_outputs:/dir_outputs hjibarram/ifu_mocks 'Config_name:dir_outputs/config_name'")
    sifu.sycall("echo python mocks.py 'Config_name:config_name'")
    sifu.sycall("echo You must put Config_name: before config_name")
    sifu.sycall("echo config_name = name of the configuration file")
    sifu.sycall("echo  ")
    sifu.sycall("echo Running with DEFAULT VALUES:")
    sifu.sycall("echo docker run -p 8000:8000 hjibarram/ifu_mocks 'Default'")
    sifu.sycall("echo python mocks.py 'Default'")
    sifu.sycall("echo  ")
    sifu.sycall("echo  ")
    sifu.sycall("echo  Default configuration file example:")
    sifu.sycall("echo  ")            
    sifu.sycall("echo basename    artsp8-")
    sifu.sycall("echo typef1      MaNGA")
    sifu.sycall("echo id1         A-0")
    sifu.sycall("echo dir1        ../")
    sifu.sycall("echo modt        1")
    sifu.sycall("echo file_red    sp8/star_out_0.dat")
    sifu.sycall("echo file_gas    sp8/Gas_out_0.dat")
    sifu.sycall("echo ang         0.0")
    sifu.sycall("echo fib_n       7")
    sifu.sycall("echo SN          15.0")
    sifu.sycall("echo Fluxm       26.83")
    sifu.sycall("echo psf         1.43")
    sifu.sycall("echo fov_p       0")
    sifu.sycall("echo fov         0")
    sifu.sycall("echo cam         0")
    sifu.sycall("echo n_pix       440")
    sifu.sycall("echo rx          142135.5")
    sifu.sycall("echo vx          [-0.7156755,-0.5130859,0.4738687]")
    sifu.sycall("echo vy          [0.6984330,-0.5257526,0.4855672]")
    sifu.sycall("echo vz          [0.0000000,0.6784741,0.7346244]")
    sifu.sycall("echo template_0  libs/templete_bc03_2.fits")
    sifu.sycall("echo template_1  libs/gsd61_156.fits")
    sifu.sycall("echo template_2  libs/templete_gas.fits")
    sifu.sycall("echo template_3  libs/templete_bc03_5.fits")
    sifu.sycall("echo Om          0.2726")
    sifu.sycall("echo Lam         0.7274")
    sifu.sycall("echo ho          0.704")
    sifu.sycall("echo redo        0")
    sys.exit()    
if len(sys.argv) >= 2:
    var1=sys.argv[1].split(';')
    var1=filter(None,var1)
    if len(var1) == 1:
        if "Default" in var1[0] or "default" in var1[0]:
            def_t=1
            sifu.sycall("echo running with the default values: The face-on Sp8D galaxy observed with the MaNGA like IFU of 127-fibers.")
        elif "Config" in var1[0] or "config" in var1[0]:
            conf_t=1
            cfile=var1[0].split(":")
            cfile=filter(None,cfile)
            c_file=cfile[1]
            sifu.sycall("echo running with the configuration file")
            sifu.sycall("echo The configuration file given is: "+c_file)
        elif "Help" in var1[0] or "help" in var1[0]:
            sifu.sycall('echo USAGE:')
            sifu.sycall("echo Running Using Docker:")
            sifu.sycall("echo docker pull hjibarram/ifu_mocks")
            sifu.sycall("echo docker run -p 8000:8000 -v '`pwd`'/dir_outputs:/dir_outputs  hjibarram/ifu_mocks 'id_name;base_name;dir_outputs;mod;type_ifu' 'file_sim_stars;file_sim_gas' 'fib_n;angle;SN;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
            sifu.sycall("echo Running Using Python:")
            sifu.sycall("echo python mocks.py 'id_name;base_name;dir_outputs;mod;type_ifu' 'file_sim_stars;file_sim_gas' 'fib_n;angle;SN;Flux_lim;psf' 'n_pix;r_0;fov_p;fov;cam' 'ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3' 'teplate1;template2;template3;template4' 'Om,Ol,ho' 'redo'")
            sifu.sycall("echo  ")
            sifu.sycall("echo key words:")
            sifu.sycall("echo id_name = name of the mock")
            sifu.sycall("echo base_name = base name of the mock: base_name-id_name")
            sifu.sycall("echo dir_outputs = output directory")
            sifu.sycall("echo mod = Mode of the mocks: mod=0 implies full set of mocks, mod=1 implies photometry IFU and single fiber only, mod=2 implies only IFU mocks, mod=3 implies only SDSS single fiber mocks")
            sifu.sycall("echo type_ifu = Type of the IFU mock: 'MaNGA' for MaNGA type, 'CALIFA' for CALIFA type and 'MUSE' for muse type")
            sifu.sycall("echo file_sim_stars = File for the stellar particules of the desired hidrodynamical simulation")
            sifu.sycall("echo file_sim_gas = File for the gas cells of the desired hidrodynamical simulation")
            sifu.sycall("echo fib_n = number of fibers per octagon size")
            sifu.sycall("echo angle = angle between the reference vector e1 and the observer")
            sifu.sycall("echo SN = Signal to noise ratio")
            sifu.sycall("echo  ")
            sifu.sycall("echo Running with Configuration file:")
            sifu.sycall("echo docker run -p 8000:8000 -v '`pwd`'/dir_outputs:/dir_outputs hjibarram/ifu_mocks 'Config_name:dir_outputs/config_name'")
            sifu.sycall("echo python mocks.py 'Config_name:config_name'")
            sifu.sycall("echo You must put Config_name: before config_name")
            sifu.sycall("echo config_name = name of the configuration file")
            sifu.sycall("echo  ")
            sifu.sycall("echo Running with DEFAULT VALUES:")
            sifu.sycall("echo docker run -p 8000:8000 hjibarram/ifu_mocks 'Default'")
            sifu.sycall("echo python mocks.py 'Default'")
            sifu.sycall("echo  ")
            sifu.sycall("echo  ")
            sifu.sycall("echo  Default configuration file example:")
            sifu.sycall("echo  ")            
            sifu.sycall("echo basename    artsp8-")
            sifu.sycall("echo typef1      MaNGA")
            sifu.sycall("echo id1         A-0")
            sifu.sycall("echo dir1        ../")
            sifu.sycall("echo modt        1")
            sifu.sycall("echo file_red    sp8/star_out_0.dat")
            sifu.sycall("echo file_gas    sp8/Gas_out_0.dat")
            sifu.sycall("echo ang         0.0")
            sifu.sycall("echo fib_n       7")
            sifu.sycall("echo SN          15.0")
            sifu.sycall("echo Fluxm       26.83")
            sifu.sycall("echo psf         1.43")
            sifu.sycall("echo fov_p       0")
            sifu.sycall("echo fov         0")
            sifu.sycall("echo cam         0")
            sifu.sycall("echo n_pix       440")
            sifu.sycall("echo rx          142135.5")
            sifu.sycall("echo vx          [-0.7156755,-0.5130859,0.4738687]")
            sifu.sycall("echo vy          [0.6984330,-0.5257526,0.4855672]")
            sifu.sycall("echo vz          [0.0000000,0.6784741,0.7346244]")
            sifu.sycall("echo template_0  libs/templete_bc03_2.fits")
            sifu.sycall("echo template_1  libs/gsd61_156.fits")
            sifu.sycall("echo template_2  libs/templete_gas.fits")
            sifu.sycall("echo template_3  libs/templete_bc03_5.fits")
            sifu.sycall("echo Om          0.2726")
            sifu.sycall("echo Lam         0.7274")
            sifu.sycall("echo ho          0.704")
            sifu.sycall("echo redo        0")
            sys.exit()
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
if def_t == 1:
    basename='artsp8-'
    typef1="MaNGA"
    id1='A-0'
    dir1='../'
    modt=1
    file_red="sp8/star_out_0.dat"
    file_gas="sp8/Gas_out_0.dat"
    ang=0.0
    fib_n=7
    SN=15.0
    Fluxm=26.83
    psf=1.43
    fov_p=0
    fov=0
    cam=0
    n_pix=440
    rx=142135.5
    vx=[-0.7156755,-0.5130859,0.4738687]
    vy=[0.6984330,-0.5257526,0.4855672]
    vz=[0.0000000,0.6784741,0.7346244]
    template_0="libs/templete_bc03_2.fits"
    template_1="libs/gsd61_156.fits"
    template_2="libs/templete_gas.fits"
    template_3="libs/templete_bc03_5.fits"
    Om=0.2726
    Lam=0.7274
    ho=0.704
    redo=0
if conf_t == 1:
    if pt.exists(c_file) == True:
        k1=0
        k2=0
        k3=0
        k4=0
        k5=0
        k6=0
        k7=0
        k8=0
        k9=0
        k10=0
        k11=0
        k12=0
        k13=0
        k14=0
        k15=0
        k16=0
        k17=0
        k18=0
        k19=0
        k20=0
        k21=0
        k22=0
        k23=0
        k24=0
        k25=0
        k26=0
        k27=0
        k28=0        
        f=open(c_file,"r")
        for line in f:
            if not "#" in line:
                data=line.replace("\n","").split(" ")
                data=filter(None,data)
                if "basename" in data[0]:
                    basename=data[1]
                    k1=1
                if "typef1" in data[0]:
                    typef1=data[1]
                    k2=1
                if "id1" in data[0]:
                    id1=data[1]
                    k3=1
                if "dir1" in data[0]:
                    dir1=data[1]
                    k4=1
                if "modt" in data[0]:
                    modt=np.int(data[1])
                    k5=1
                if "file_red" in data[0]:
                    file_red=data[1]
                    k6=1
                if "file_gas" in data[0]:
                    file_gas=data[1]
                    k7=1
                if "ang" in data[0]:
                    ang=np.float(data[1])
                    k8=1
                if "fib_n" in data[0]:
                    fib_n=np.int(data[1])
                    k9=1
                if "SN" in data[0]:
                    SN=np.float(data[1])
                    k10=1
                if "Fluxm" in data[0]:
                    Fluxm=np.float(data[1])
                    k11=1
                if "psf" in data[0]:
                    psf=np.float(data[1])
                    k12=1
                if "fov_p" in data[0]:
                    fov_p=np.float(data[1])
                    k13=1
                if "fov" in data[0] and not "fov_p" in data[0]:
                    fov=np.float(data[1])
                    k14=1
                if "cam" in data[0]:
                    cam=np.float(data[1])
                    k15=1
                if "n_pix" in data[0]:
                    n_pix=np.int(data[1])
                    k16=1
                if "rx" in data[0]:
                    rx=np.float(data[1])
                    k17=1
                if "vx" in data[0]:
                    st_t=data[1].replace("[","").replace("]","")
                    st_n=st_t.split(",")
                    st_n=filter(None,st_n)
                    vt1=np.float(st_n[0])
                    vt2=np.float(st_n[1])
                    vt3=np.float(st_n[2])
                    vx=[vt1,vt2,vt3]
                    k18=1
                if "vy" in data[0]:
                    st_t=data[1].replace("[","").replace("]","")
                    st_n=st_t.split(",")
                    st_n=filter(None,st_n)
                    vt1=np.float(st_n[0])
                    vt2=np.float(st_n[1])
                    vt3=np.float(st_n[2])
                    vy=[vt1,vt2,vt3]
                    k19=1
                if "vz" in data[0]:
                    st_t=data[1].replace("[","").replace("]","")
                    st_n=st_t.split(",")
                    st_n=filter(None,st_n)
                    vt1=np.float(st_n[0])
                    vt2=np.float(st_n[1])
                    vt3=np.float(st_n[2])
                    vz=[vt1,vt2,vt3]
                    k20=1
                if "template_0" in data[0]:
                    template_0=data[1]
                    k21=1
                if "template_1" in data[0]:
                    template_1=data[1]
                    k22=1
                if "template_2" in data[0]:
                    template_2=data[1]
                    k23=1
                if "template_3" in data[0]:
                    template_3=data[1]
                    k24=1
                if "Om" in data[0]:
                    Om=np.float(data[1])
                    k25=1
                if "Lam" in data[0]:
                    Lam=np.float(data[1])
                    k26=1
                if "ho" in data[0]:
                    ho=np.float(data[1])
                    k27=1
                if "redo" in data[0]:
                    redo=np.int(data[1])
                    k28=1
        f.close()
        if k1 == 0:
             basename='artsp8-'
        if k2 == 0:
            typef1="MaNGA"
        if k3 == 0:
            id1='A-0'
        if k4 == 0:
            dir1='../'
        if k5 == 0:
            modt=1
        if k6 == 0:
            file_red="sp8/star_out_0.dat"
        if k7 == 0:
            file_gas="sp8/Gas_out_0.dat"
        if k8 == 0:
            ang=0.0
        if k9 == 0:
            fib_n=7
        if k10 == 0:
            SN=15.0
        if k11 == 0:
            Fluxm=26.83
        if k12 == 0:
            psf=1.43
        if k13 == 0:
            fov_p=0
        if k14 == 0:
            fov=0
        if k15 == 0:
            cam=0
        if k16 == 0:
            n_pix=440
        if k17 == 0:
            rx=142135.5
        if k18 == 0:
            vx=[-0.7156755,-0.5130859,0.4738687]
        if k19 == 0:
            vy=[0.6984330,-0.5257526,0.4855672]
        if k20 == 0:
            vz=[0.0000000,0.6784741,0.7346244]
        if k21 == 0:
            template_0="libs/templete_bc03_2.fits"
        if k22 == 0:
            template_1="libs/gsd61_156.fits"
        if k23 == 0:
            template_2="libs/templete_gas.fits"
        if k24 == 0:
            template_3="libs/templete_bc03_5.fits"
        if k25 == 0:
            Om=0.2726
        if k26 == 0:
            Lam=0.7274
        if k27 == 0:
            ho=0.704
        if k28 == 0:
            redo=0   
    else:
        sifu.sycall("echo The configuration file does not exists")
        sys.exit()
file_out='mock_mass_ill_0.out'
file_out_f='mock_mass_ill.out'
sifu.mock_sp(fib_n,ang,modt=modt,template_1=template_1,template_2=template_2,template_3=template_3,template=template_0,n_pix=n_pix,fov_p=fov_p,fov_m=fov,rx=rx,Om=Om,Lam=Lam,ho=ho,cam=cam,vx=vx,vy=vy,vz=vz,base_name=basename,typef1=typef1,id1=id1,psf=psf,redo=redo,SN=SN,Fluxm=Fluxm,dir1=dir1,file_red=file_red,file_gas=file_gas,file_out=file_out,file_out_f=file_out_f)