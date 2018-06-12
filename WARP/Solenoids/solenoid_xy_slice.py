from warp import *
import matplotlib.pyplot as plt
import numpy as np

__author__ = "Daniel Winklehner"
__doc__ = """
Warp Script with a simple solenoid lattice to demonstrate space charge compensation.

USPAS Jan 2018
"""
    
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
emit = 30.0e-7        # rms edge emittance [m-rad]
i_beam = 10.0 * mA    # beam current per species [A]
r_x = 10.0 * mm       # beam edge in x [m]
r_y = 10.0 * mm       # beam edge in y [m]
r_xp = 10.0e-3        # initial beam divergence [rad]
r_yp = r_xp           # initial beam divergence [rad]
NParticles = 10000

# Define ion species
protons = Species(type=Proton, charge_state = +1, name = "Protons")
#hy2plus = Species(type=Dihydrogen, charge_state = +1, name = "H2+")

## ASSIGNMENT for LEBT Task 1: Uncomment H2+ above and add it to the list below
beam_species = [protons]

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

addnewsolenoid(zi=solenoid_zi[0],
               zf=solenoid_ze[0],
               ri=solenoid_radius,
               maxbz=0.115)  # 0.115 T for p+, 0.17 T for H2+

addnewsolenoid(zi=solenoid_zi[1],
               zf=solenoid_ze[1],
               ri=solenoid_radius,
               maxbz=0.075)  # 0.075 T for p+, 0.125 T for H2+

addnewsolenoid(zi=solenoid_zi[2],
               zf=solenoid_ze[2],
               ri=solenoid_radius,
               maxbz=0.12)  # 0.12 T for p+, 0.14 T for H2+

 
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
top.npmax = NParticles
# ASSIGNEMENT for LEBT Task 1 change type of distribution
w3d.distrbtn = "KV"          # initial distribution "KV" for K-V, "WB" for Waterbag

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

# Prepare plot axes and plot initial cross-section and phase space
fig, ax = multiplot_init(beam_species)

# Main simulation loop
for i in range(nsteps):

    step()  # do one step

    for j, species in enumerate(beam_species):
        x_rms[j][i+1] = 2.0 * np.std(species.getx())

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

plt.tight_layout()
plt.show()
