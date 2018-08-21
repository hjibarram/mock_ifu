# mock_ifu
A tool for mocking IFU observations from hydrodynamical n-body simulations
# USAGE:
# Running Using Docker:
docker pull hjibarram/ifu_mocks
docker run -p 8000:8000 -v `pwd`/dir_outputs:/dir_outputs hjibarram/ifu_mocks id_name;base_name;dir_outputs;mod;type_ifu file_sim_stars;file_sim_gas fib_n;angle;SN;Flux_lim;psf n_pix;r_0;fov_p;fov;cam ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3 teplate1;template2;template3;template4 Om,Ol,ho redo
# Running Using Python:
python mocks.py id_name;base_name;dir_outputs;mod;type_ifu file_sim_stars;file_sim_gas fib_n;angle;SN;Flux_lim;psf n_pix;r_0;fov_p;fov;cam ex1,ey1,ez1;ex2,ey2,ez2;ex2,ey3,ez3 teplate1;template2;template3;template4 Om,Ol,ho redo

# key words:
id_name = name of the mock
base_name = base name of the mock: base_name-id_name
dir_outputs = output directory
mod = Mode of the mocks: mod=0 implies full set of mocks, mod=1 implies photometry IFU and single fiber only, mod=2 implies only IFU mocks, mod=3 implies only SDSS single fiber mocks
type_ifu = Type of the IFU mock: MaNGA for MaNGA type, CALIFA for CALIFA type and MUSE for muse type
file_sim_stars = File for the stellar particules of the desired hidrodynamical simulation
file_sim_gas = File for the gas cells of the desired hidrodynamical simulation
fib_n = number of fibers per octagon size
angle = angle between the reference vector e1 and the observer
SN = Signal to noise ratio

# Running with Configuration file:
docker run -p 8000:8000 -v '`pwd`'/dir_outputs:/dir_outputs hjibarram/ifu_mocks Config_name:dir_outputs/config_name
python mocks.py Config_name:config_name
You must put Config_name: before config_name
config_name = name of the configuration file

# Running with DEFAULT VALUES:
docker run -p 8000:8000 hjibarram/ifu_mocks Default
python mocks.py Default


# Default configuration file example:

basename     artsp8-
typef1       MaNGA
id1          A-0
dir1         ../
modt         1
file_red     sp8/star_out_0.dat
file_gas     sp8/Gas_out_0.dat
ang          0.0
fib_n        7
SN           15.0
Fluxm        26.83
psf          1.43
fov_p        0
fov          0
cam          0
n_pix        440
rx           142135.5
vx           [-0.7156755,-0.5130859,0.4738687]
vy           [0.6984330,-0.5257526,0.4855672]
vz           [0.0000000,0.6784741,0.7346244]
template_0   libs/templete_bc03_2.fits
template_1   libs/gsd61_156.fits
template_2   libs/templete_gas.fits
template_3   libs/templete_bc03_5.fits
Om           0.2726
Lam          0.7274
ho           0.704
redo         0