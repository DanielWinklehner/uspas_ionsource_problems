"""
Example Pierce diode calculation.
Hot plate source emitting singly ionized potassium
"""
from warp import *

# --- Set four-character run id, comment lines, user's name.
top.pline2   = "Pierce diode example"
top.pline1   = "Injected beam. Semi-Gaus."
top.runmaker = "J.-L. Vay"

# --- Invoke setup routine for the plotting
setup()

# --- Set the dimensionality
w3d.solvergeom = w3d.RZgeom

# ------------------------------------------------------------------------------
# setup conducting plates
# ------------------------------------------------------------------------------

# --- Basic parameters
channel_radius     = 15.*cm
extractor_voltage  = -93.*kV
diode_current      = 1.

# --- Setup source plate
source_radius = channel_radius

# --- Setup diode aperture plate
zplate = 8.*cm # --- plate location

# ------------------------------------------------------------------------------
# setup beam and injection
# ------------------------------------------------------------------------------

# --- Setup simulation species
beam = Species(type=Potassium, charge_state=+1, name='beam')

# --- Set basic beam parameters
beam.a0       = source_radius
beam.b0       = source_radius
beam.ibeam    = diode_current
derivqty()

# --- Specify injection of the particles
top.inject = 1                # 1 means constant; 2 means space-charge limited injection
top.npinject = 15              # Approximate number of particles injected each step
top.vinject = 0.

# --- If using the RZ geometry, set so injection uses the same geometry
w3d.l_inj_rz = (w3d.solvergeom == w3d.RZgeom)

# ------------------------------------------------------------------------------
# setup grid and boundary conditions
# ------------------------------------------------------------------------------

# --- Set boundary conditions
# ---   for field solve
w3d.bound0  = dirichlet
w3d.boundnz = neumann
w3d.boundxy = neumann
# ---   for particles
top.pbound0  = absorb
top.pboundnz = absorb
top.pboundxy = reflect

# --- Set field grid size
w3d.xmmin = w3d.ymmin = -channel_radius
w3d.xmmax = w3d.ymmax = +channel_radius
w3d.zmmin = 0.; w3d.zmmax = zplate*1.

# --- Field grid dimensions - note that nx and ny must be even.
w3d.nx = w3d.ny = 4
w3d.nz = 128

# --- Set the time step size. This needs to be small enough to satisfy the Courant limit.
dz = (w3d.zmmax - w3d.zmmin)/w3d.nz
vzfinal = sqrt(2.*abs(extractor_voltage)*jperev/beam.mass)
top.dt = 0.4*(dz/vzfinal)

# ------------------------------------------------------------------------------
# setup field solver, install conductors 
# ------------------------------------------------------------------------------

# --- Set up fieldsolver
f3d.mgtol = 1.e-1 # Multigrid solver convergence tolerance, in volts
solver = MultiGrid2D()
registersolver(solver)

# --- Create source conductors
source = ZPlane(zcent=w3d.zmmin,zsign=-1.,voltage=0.)
installconductor(source, dfill=largepos)

# --- Create ground plate
plate = ZPlane(voltage=extractor_voltage, zcent=zplate)
installconductor(plate,dfill=largepos)

# --- Setup the particle scraper
scraper = ParticleScraper([source, plate])

# ------------------------------------------------------------------------------
# load main package and initialize 
# ------------------------------------------------------------------------------

# --- Set pline1 to include appropriate parameters
top.pline1 = ("Injected beam. Parallel plates %dx%d. npinject=%d, dt=%d"%
             (w3d.nx, w3d.nz, top.npinject, top.dt))

# --- Generate the PIC code (allocate storage, load ptcls, t=0 plots, etc.)
package("w3d")
generate()

# ------------------------------------------------------------------------------
# setup plots
# ------------------------------------------------------------------------------

# --- Open up plotting windows
window(1, hcp='current.cgm', dump=1)
winon()

def plzprofiles(l_CL=True):

    # --- get phi averaged over x
    phi = ave(solver.phi)[0][1:-1]
    # --- get rho averaged over x
    rho = ave(solver.rho)[0]
    # --- get Ez from phi
    Ez   = -(phi[1:]-phi[:-1])/w3d.dz
    # --- compute location of Ez (between mesh nodes)
    zEz  = 0.5*w3d.dz+w3d.zmesh[:-1]

    if l_CL:
        # --- Child-Langmuir Law
        print("Write expression of Child-Langmuir law below") 
#        JCL =
        Program_breaks_here # Follow instructions on previous line and remove this line.
        # --- phi from Child-Langmuir as a function of z (w3d.zmesh)
        print("write expression of Child-Langmuir potential below")
#        phiCL = 
        Program_breaks_here # Follow instructions on previous line and remove this line.
        # --- velocity from Child-Langmuir as a function of phiCL
        vCL   = sqrt(-2.*beam.charge*phiCL/beam.mass)
        # --- computes charge density from Child-Langmuir from JCL and vCL
        rhoCL = zeros(w3d.nz+1)
        rhoCL[1:] = JCL/vCL[1:]
        # --- get Child-Langmuir Ez from phiCL
        EzCL = -(phiCL[1:]-phiCL[:-1])/w3d.dz
        # --- compute current from current density
        ICL = pi*top.a0**2*JCL
        print(('Maximum current from Child-Langmuir law = ',ICL))
    
    # --- plot currents
    plsys(3);
    I = diode_current
    plg([I,I], [w3d.zmmin, w3d.zmmax],width=4,color=blue,type='dot')
    if l_CL:plg([ICL,ICL], [w3d.zmmin, w3d.zmmax],width=4,color=red)
    pzcurr(width=3,linetype='dash',titles=0)
    if l_CL:
        limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., max(diode_current,ICL)*1.5)
    else:
        limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., diode_current*1.5)
    ptitles('Current [A]','Z [m]','')
    
    # --- plot charge density
    plsys(4)
    if l_CL:pla(rhoCL,w3d.zmesh,width=4,color=red)
    pla(rho,w3d.zmesh,width=3,type='dash')
    ptitles('Rho [e/m^3__]','Z [m]','')
    
    # --- plot electric field
    plsys(5)
    if l_CL:pla(EzCL,zEz,color=red,width=4)
    pla(Ez,zEz,type='dash',width=3)
    if l_CL:
        limits(0.,w3d.zmmax,0.,maxnd(EzCL))
    else:
        limits(0.,w3d.zmmax,0.,maxnd(Ez))
    ptitles('E_z^^ [V/m]','Z [m]','')
    
    # --- plot potential
    plsys(6)
    if l_CL:pla(phiCL,w3d.zmesh,width=4,color=red)
    pla(phi,w3d.zmesh,width=3,type='dash')
    ptitles('Phi [V]','Z [m]','')
        
def beamplots():
    # --- select window 0
    window(0)
    # --- clear frame
    fma()
    # --- plot potential
    pcphizr(titles=False,filled=1,contours=100)
    # --- plot particles
    ppzr(titles=False,msize=2)
    # --- sets plot limits
    limits(w3d.zmminglobal, w3d.zmmaxglobal, 0., channel_radius)
    # --- add title
    ptitles('Hot plate source', 'Z (m)', 'R (m)')
    # --- refresh screen for live plots
    refresh()

    # --- select window 0
    window(1)
    # --- clear frame
    fma()
    # --- plots profiles of current, Rho, Ez and Phi vs z
    plzprofiles(l_CL=False)
    # --- refresh screen for live plots
    refresh()

# --- Call beamplots after every 25 steps
@callfromafterstep
def makebeamplots():
    if top.it%25 == 0:
        beamplots()

# ------------------------------------------------------------------------------
# Main loop
# ------------------------------------------------------------------------------

step(2000)

# --- Make sure that last plot frames get sent to the cgm file
window(0)
hcp()
window(1)
hcp()
