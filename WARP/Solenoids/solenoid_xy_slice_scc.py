from warp import *
import matplotlib.pyplot as plt
import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = """
Warp Script with a simple solenoid lattice to demonstrate space charge compensation.

USPAS Jan 2018
"""

def coulomblog(ne, q, vb):
    """
    returns the coulomb logarithm
    """
    b_90 = q * echarge**2 / 4 / np.pi / eps0 / emass / vb**2
    om_pe = np.sqrt(ne * echarge**2 / (emass * eps0))
    b_max = vb / om_pe
    return 4.0 * np.pi * np.log(b_max / b_90)

def set_compensated_beam():
    """
    Calculate the neutralization according to Gabovich's formula
    """
    global fe

    V = 30000  # beam energy in eV
    pressure = 1e-7  # Torr
    sigma_i = np.array([5e-20])  # ion production cross-section of 30 keV p on (m^2)
    sigma_e = np.array([1e-20])  # electron production cross-section of 30 keV p on H2 (m^2)
    n0 = 133.322 * pressure / boltzmann / 300  # residual gas density at room temp (300 K)
    Tr = 300 # Room temperature in K (~300 K)
    PHIi = 15.51
    Ti = 5  # secondary ion temperature in eV (guess)
    mi = 2  # ~Mass of H2 in amu  (assuming residual gas is H2)
    vi = np.sqrt(echarge * Ti * 2 / (mi * amu))
    I = np.array([i_beam / NParticles * protons.getn()])
    q = np.array([1])
    m = np.array([1.0072766])
    vb = np.array([top.vbeam])
    # print("vb = ", vb, "m/s")
    r0 = 2.0 * np.std(protons.getx())
    # print("r0 = ", 1000.0*r0, "mm")
    # print("I = ", 1000*I," mA")

    delta_phi_full = np.sum(I / (4 * np.pi * vb * eps0))
    #print("delta_phi_full = ", delta_phi_full)

    nb = I / (echarge * q * vb * np.pi * r0**2)
    ni = r0 * n0 * nb * vb * sigma_i / (2 * vi)
    ne = np.sum(q * nb + ni)
    #print("ne = ", ne)
    #print("nb = ", nb)

    L = coulomblog(ne, q, vb)
    #print("L = ", L)

    chi = np.sqrt(ne / n0 * np.sum(I * m * amu * L) / np.sum(I * sigma_e / q) * \
                  echarge**2.0 / (4.0 * np.pi * eps0)**2.0 * PHIi * 3.0 / emass / V) / delta_phi_full

    neut = (1.0 + 0.5 * chi**2.0 - 0.5 * chi * np.sqrt(4.0 + chi**2.0))
    
    if neut >= 0.999:
        neut = 0.999

    #print("fe = ", neut)

    fe.append(neut)

    # Reset weigths and beam current according to space charge compensation
    Weight = I * top.sp_fract / (top.vbeam_s * echarge * top.zion_s * protons.getn())
    top.ibeam_s = (1.0 - neut) * I
    top.pgroup.sw = (1.0 - neut) * Weight
    
def multiplot_init(beam_species):

    if not isinstance(beam_species, list):
        beam_species = [beam_species]

    ax1 = plt.subplot2grid((2, 4), (0, 0), colspan=4)
    ax2 = plt.subplot2grid((2, 4), (1, 0))
    ax3 = plt.subplot2grid((2, 4), (1, 1))
    ax4 = plt.subplot2grid((2, 4), (1, 2))
    ax5 = plt.subplot2grid((2, 4), (1, 3))

    for species in beam_species:

        ax2.scatter(1e3*species.getx(), 1e3*species.gety(),
                    s=0.5, label="{}".format(species.name))
        ax3.scatter(1e3*species.getx(), 1e3*species.getxp(),
                    s=0.5, label="{}".format(species.name))

    plt.sca(ax2)
    plt.legend
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.title("Init Cross-Section")
    ax2.set_aspect(1)

    plt.sca(ax3)
    plt.legend
    plt.xlabel("x (mm)")
    plt.ylabel("x' (mrad)")
    plt.title("Init XX' Phase Space")
    ax3.set_aspect(1)
    
    plt.sca(ax4)
    plt.legend
    plt.xlabel("x (mm)")
    plt.ylabel("y (mm)")
    plt.title("Final Cross-Section")
    ax4.set_aspect(1)
    
    plt.sca(ax5)
    plt.legend
    plt.xlabel("x (mm)")
    plt.ylabel("x' (mrad)")
    plt.title("Final XX' Phase Space")
    ax5.set_aspect(1)
    
    return plt.gcf(), [ax1, ax2, ax3, ax4, ax5]


def multiplot_fin(ax, beam_species):

    if not isinstance(beam_species, list):
        beam_species = [beam_species]

    for species in beam_species:

        ax[3].scatter(1e3*species.getx(), 1e3*species.gety(),
                    s=0.5, label="{}".format(species.name))
        ax[4].scatter(1e3*species.getx(), 1e3*species.getxp(),
                    s=0.5, label="{}".format(species.name))

    return 0

# Set informational labels included on all output cgm plots.
top.pline2   = "xy-slice simulation: solenoid channel"
top.pline1   = " "

setup()

top.runmaker = "USPAS"

# Define some user variables
e_kin = 30.0 * keV    # kinetic energy [eV]
emit = 10.0e-7        # rms edge emittance [m-rad]
i_beam = 10.0 * mA     # beam current per species [A]
r_x = 10.0 * mm       # beam edge in x [m]
r_y = 10.0 * mm       # beam edge in y [m]
r_xp = 10.0e-3        # initial beam divergence [rad]
r_yp = r_xp           # initial beam divergence [rad]
NParticles = 10000

# Define ion species
protons = Species(type=Proton, charge_state = +1, name = "Protons")
beam_species = [protons]  # Cave: Even though the space charge compensation can only handle one species,
                          # we have to make a list here so that the rest of the script works

for beam in beam_species:
    beam.ekin     = e_kin       # kinetic energy of beam particle [eV]
    beam.vbeam    = 0.0         # beam axial velocity [m/sec]
    beam.ibeam    = i_beam      # beam current [A]
                             
    beam.emitx    = emit        # beam x-emittance, rms edge [m-rad]
    beam.emity    = emit        # beam y-emittance, rms edge [m-rad]
    beam.vthz     = 0.0         # axial velocity spread [m/sec]

    # This routine will calculate vbeam and other quantities.
    derivqty()

    # Beam centroid and envelope initial conditions
    beam.x0  = 0.0   # initial x-centroid xc = <x> [m]
    beam.y0  = 0.0   # initial y-centroid yc = <y> [m]
    beam.xp0 = 0.0   # initial x-centroid angle xc' = <x'> = d<x>/ds [rad]
    beam.yp0 = 0.0   # initial y-centroid angle yc' = <y'> = d<y>/ds [rad]
    beam.a0  = r_x   # initial x-envelope edge a = 2*sqrt(<(x-xc)^2>) [m]
    beam.b0  = r_y   # initial y-envelope edge b = 2*sqrt(<(y-yc)^2>) [m]
    beam.ap0 = r_xp  # initial x-envelope angle ap = a' = d a/ds [rad]
    beam.bp0 = r_yp  # initial y-envelope angle bp = b' = d b/ds [rad]

# Set up solenoid lattice
run_length = 3.6
drift_length = 0.6
solenoid_length = 0.4
solenoid_radius = 7.5e-2

solenoid_zi = [drift_length + i * solenoid_length + i * drift_length for i in range(3)]
solenoid_ze = [drift_length + (i + 1) * solenoid_length + i * drift_length for i in range(3)]

# Solenoids are created here
addnewsolenoid(zi=solenoid_zi[0],
               zf=solenoid_ze[0],
               ri=solenoid_radius,
               maxbz=0.115)

addnewsolenoid(zi=solenoid_zi[1],
               zf=solenoid_ze[1],
               ri=solenoid_radius,
               maxbz=0.075)

addnewsolenoid(zi=solenoid_zi[2],
               zf=solenoid_ze[2],
               ri=solenoid_radius,
               maxbz=0.12)

# Pipe in the solenoid transport
pipe = ZCylinderOut(radius=solenoid_radius, zlower=0.0, zupper=run_length)
conductors=pipe

# set up particle termination at cylindrical wall
top.prwall = solenoid_radius

n_grid = 256  # number of grid cells in x and y direction
w3d.nx = n_grid
w3d.ny = n_grid

w3d.xmmax =  solenoid_radius  # x-grid max limit [m]
w3d.xmmin = -solenoid_radius  # x-grid min limit [m]
w3d.ymmax =  solenoid_radius  # y-grid max limit [m]
w3d.ymmin = -solenoid_radius  # y-grid min limit [m]

# Particle distribution options
w3d.distrbtn = "KV" # initial distribution "KV" for K-V, "WB" for Waterbag

top.npmax = NParticles

# Random number options to use in loading
w3d.xrandom  = "digitrev"    # load x,y,z  with digitreverse random numbers
w3d.vtrandom = "digitrev"    # load vx, vy with digitreverse random numbers
w3d.vzrandom = "digitrev"    # load vz     with digitreverse random numbers
w3d.cylinder = True          # load a cylinder

top.lrelativ   =  False    # turn off relativistic kinematics
top.relativity = 0         # turn off relativistic self-field correction
wxy.ds = 1.0e-3            # ds for part adv [m]
wxy.lvzchang = True        # Use iterative stepping, which is needed if the vz changes
top.ibpush   = 2           # magnetic field particle push: 0 - off, 1 - fast, 2 - accurate 

# Setup field solver using 2d multigrid field solver.
w3d.boundxy = 0              # Neuman boundary conditions on edge of grid.
w3d.solvergeom = w3d.XYgeom  # fieldsolve type to 2d multigrid 

# Generate the xy PIC code.  In the generate, particles are allocated and
# loaded consistent with initial conditions and load parameters
# set previously.  Particles are advanced with the step() command later
# after various diagnostics are setup.
package("wxy")
generate()

# Install conducting aperture on mesh
installconductors(conductors, dfill=largepos)

# Carry out explicit fieldsolve after generate to include conducing pipe
# with initial beam
fieldsolve()

# Some local runtime variables
nsteps = int(np.ceil(run_length/wxy.ds))

x_rms = []
for species in beam_species:
    x_rms.append(np.empty(nsteps+1))
    x_rms[-1][0] = 2.0 * np.std(species.getx())
    
zpos = 0.0

## ASSIGNMENT for WARP Task 2c: Uncomment to include space charge compensation
#fe = []
#installbeforestep(set_compensated_beam)

# Prepare plot axes and plot initial cross-section and phase space
fig, ax = multiplot_init(beam_species)

# Main simulation loop
for i in range(nsteps):

    step()  # do one step

    for j, species in enumerate(beam_species):
        x_rms[j][i+1] = 2.0 * np.std(species.getx())

    # zpos = np.round(top.zbeam, 8)

multiplot_fin(ax, beam_species)
plt.sca(ax[0])

for i, species in enumerate(beam_species):
    plt.plot(np.arange(nsteps+1) * wxy.ds, x_rms[i], label="{}, X".format(species.name))

for zi, ze in zip(solenoid_zi, solenoid_ze):
    plt.axvline(zi, ls="--", color='black')
    plt.axvline(ze, ls="--", color='black')

plt.xlabel("z (m)")
plt.ylabel("x (m)")
plt.title("2*RMS Envelopes")
plt.legend()
plt.xlim(0.0, run_length)
plt.ylim(0.0, None)

## ASSIGNMENT for WARP Task 2c: Uncomment to see fe plot
#ax[0].twinx()
#plt.plot(np.arange(nsteps) * wxy.ds, np.array(fe), color='green')

plt.tight_layout()
plt.show()
