#include <fstream>
#include <iomanip>
#include <limits>
#include "epot_bicgstabsolver.hpp"
#include "meshvectorfield.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "dxf_solid.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "random.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "transformation.hpp"
#include "convergence.hpp"
#include <stdlib.h>
#include <iostream>   // std::cout
#include <string>     // std::string, std::to_string
#include <math.h>
#include "stl_solid.hpp"
#include "stlfile.hpp"

using namespace std;

#define RZ_MODE

// *** Some constants *** //
const double echarge = 1.60217662e-19; // Elementary Charge (C)
const double emass_mev = 0.5109989461; // Electron Mass (MeV)
const double amu_kg = 1.660539040e-27; // Atomic Mass Unit (kg)
const double amu_mev = 931.4940954;    // Atomic Mass Unit (MeV)
const double m_p1p_amu = 1.0073;       // proton mass (amu)
const double m_h2p_amu = 2.0151;       // H2+ mass (amu)

// *** Initial values and setup of plasma parameters *** //
bool interactive = false;    // flag for interactive plotting mode, default is false
                             // (override by command line)
bool savepart = false;       // flag for saving final particle distribution
                             // (override by command line)                 
double fe = 0.0;             // scc factor (override by command line)
double scc_onset = 143.0e-3; // beginning of space charge compensation (m)

// ** Voltages (can be overridden by command line) ** //
double VSource = 8.0e3;     // V
double VPuller = -1.0e3;     // V
double VLens = 7.0e3;       // V
double VRepeller = 0.0;     // V

// ** Geometry ** //
//string fdxf_name = "mist1a_aperture6mm_distance8mm.dxf";  // DXF filename
//string fdxf_name = "mist1b_aperture6mm_distance8mm.dxf";  // DXF filename
string fdxf_name = "mist1k_aperture4mm_distance5mm.dxf";  // DXF filename

double rmax = 32.0e-3;   // m
double zmin = -2.0e-3;   // m
double zmax = 143.0e-3;  // m
double h = 0.2e-3;  // cell size (m)
double r_source = 2.0e-3;   // m
double r_plasma = 6.0e-3;  // m

double I_p1p = 10e-3;   // proton current (A)
double I_h2p = 20e-3;   // H2+ current (A)

double Te = 3.0;  // electron temperature (eV)
double Ti = 2.0;  // ion temperature (eV)

double J_p1p = I_p1p / (r_source * r_source * M_PI);  // A/m^2
double J_h2p = I_h2p / (r_source * r_source * M_PI);  // A/m^2

double v_p1p = sqrt(echarge * (Ti + Te) / (amu_kg * m_p1p_amu));
double v_h2p = sqrt(echarge * (Ti + Te) / (amu_kg * m_h2p_amu));

double cd_p1p = J_p1p / v_p1p;  // C/m^3
double cd_h2p = J_h2p / v_h2p;  // C/m^3

double n_p1p = cd_p1p / echarge;  // 1/m^3
double n_h2p = cd_h2p / echarge;  // 1/m^3

double NPart = 120000;  // number of macroparticles per species

//double E0 = 0.5 * (Te + Ti);
double E0 = Te + Ti;

//double Up = Te*(log(2*J)-log(J*pow(2*3.14159*(1+Tt/Te)*9.109e-31/(3.3474e-27) ,.5)+(J*pow(2*3.14159*(1+Tt/Te)*9.109e-31/(.5*3.3474e-27), .5) ))) ;
double Up = Te * (log(n_p1p + n_h2p) -  log(
					    (n_p1p * sqrt(2.0 * M_PI * emass_mev / (amu_mev * m_p1p_amu) * (1 + Ti / Te) )) +
					    (n_h2p * sqrt(2.0 * M_PI * emass_mev / (amu_mev * m_h2p_amu) * (1 + Ti / Te) ))
					    )
		  );

void simu(int *argc, char ***argv) {
  
  ibsimu.message(1) << "Calculated Plasma Parameters:\n";
  ibsimu.message(1) << "  E0 = " << E0 << " eV\n";
  ibsimu.message(1) << "  Up = " << Up << " V\n";
  
#ifdef RZ_MODE // RZ
  // *** Set Mesh Parameters, Simulation space size, etc. *** //  
  Vec3D origin(zmin, 0.0, 0.0);
  Vec3D sizereq(zmax -zmin + h, rmax, 0.0);

  Int3D size (floor(sizereq[0] / h) + 1,
              floor(sizereq[1] / h) + 1,
	      1);
  
  Geometry geom(MODE_CYL, size, origin, h);
#else  // 3D
  // *** Set Mesh Parameters, Simulation space size, etc. *** //  
  Vec3D origin(-rmax, -rmax, zmin);
  Vec3D sizereq(2 * rmax, 2 * rmax, zmax - zmin + h);

  Int3D size(floor(sizereq[0] / h) + 1,
             floor(sizereq[1] / h) + 1,
             floor(sizereq[2] / h) + 1);
  
  Geometry geom(MODE_3D, size, origin, h );
#endif

  // *** Transformations for the STL files*** //
  /*
  const double scaling = 1e-3;   

  Transformation Tr_Plasma;
  Tr_Plasma.scale( Vec3D( scaling, scaling, scaling ) );
  Tr_Plasma.translate( Vec3D( 0.0, 0.0, 0.0 ) );

  Transformation Tr_Puller;
  Tr_Puller.scale( Vec3D( scaling, scaling, scaling ) );
  Tr_Puller.translate( Vec3D( 0.0, 0.0, 0.0 ) );
  
  Transformation Tr_Lens;
  Tr_Lens.scale( Vec3D( scaling, scaling, scaling ) );
  Tr_Lens.translate( Vec3D( 0.0, 0.0, 0.0 ) );

  Transformation Tr_Repeller;
  Tr_Repeller.scale( Vec3D( scaling, scaling, scaling ) );
  Tr_Repeller.translate( Vec3D( 0.0, 0.0, 0.0 ) );

  Transformation Tr_Ground;
  Tr_Ground.scale( Vec3D( scaling, scaling, scaling ) );
  Tr_Ground.translate( Vec3D( 0.0, 0.0, 0.0 ) );

  // *** Open the STL files *** //
  STLFile *fplasma = new STLFile("PlasmaAperture_6mm.stl");
  STLFile *fpuller = new STLFile( "Puller_distance8mm.stl" );
  STLFile *flens = new STLFile( "Lens_distance8mm.stl" );
  STLFile *frepeller = new STLFile( "Repeller_distance8mm.stl" );
  STLFile *fground = new STLFile( "Ground_distance8mm.stl" );

  // *** Create the solid objects from the STL files *** //
  STLSolid *plasma = new STLSolid;
  plasma->set_transformation( Tr_Plasma );
  plasma->add_stl_file( fplasma );

  STLSolid *puller = new STLSolid;
  puller->set_transformation( Tr_Puller );
  puller->add_stl_file( fpuller );
  
  STLSolid *lens = new STLSolid;
  lens->set_transformation( Tr_Lens );
  lens->add_stl_file( flens );

  STLSolid *repeller = new STLSolid;
  repeller->set_transformation( Tr_Repeller );
  repeller->add_stl_file( frepeller );

  STLSolid *ground = new STLSolid;
  ground->set_transformation( Tr_Ground );
  ground->add_stl_file( fground );
  */

  // *** Open the DXF file with the electrode cross-sections *** //
  MyDXFFile *fdxf_mist1a = new MyDXFFile;
  fdxf_mist1a->set_warning_level(1);
  fdxf_mist1a->read(fdxf_name);

  // *** Create the solid objects from the DXF file *** //
  DXFSolid *plasma = new DXFSolid(fdxf_mist1a, "plasma");
  DXFSolid *puller = new DXFSolid(fdxf_mist1a, "puller");
  DXFSolid *lens = new DXFSolid(fdxf_mist1a, "lens");
  DXFSolid *repeller = new DXFSolid(fdxf_mist1a, "repeller");
  DXFSolid *ground = new DXFSolid(fdxf_mist1a, "ground");
  
  const double scaling = 1.0e-3; 
 
  plasma->scale(scaling);
  puller->scale(scaling);
  lens->scale(scaling);
  repeller->scale(scaling);
  ground->scale(scaling);

#ifndef RZ_MODE
  // *** mapping not necessary in RZ mode *** //
  plasma->define_2x3_mapping(DXFSolid::rotz);
  puller->define_2x3_mapping(DXFSolid::rotz);
  lens->define_2x3_mapping(DXFSolid::rotz);
  repeller->define_2x3_mapping(DXFSolid::rotz);
  ground->define_2x3_mapping(DXFSolid::rotz);
#endif

  // *** Set the new solids in the geometry (RZ and 3D identical) *** //
  geom.set_solid(7, plasma);
  geom.set_solid(8, puller);
  geom.set_solid(9, lens);
  geom.set_solid(10, repeller);
  geom.set_solid(11, ground);
  
  // *** Set boundary conditions (Dirichlet = fixed Voltage) *** //
  // *** All potentials are shifted such that the source is at 0 V *** //
  geom.set_boundary(3, Bound(BOUND_NEUMANN, 0.0));      // ymin/rmin
  geom.set_boundary(4, Bound(BOUND_NEUMANN, 0.0));      // ymax/rmax

#ifdef RZ_MODE
  geom.set_boundary(1, Bound(BOUND_DIRICHLET, Up));      // xmin
  geom.set_boundary(2, Bound(BOUND_DIRICHLET, 0.0 - VSource));      // xmax
#else
  geom.set_boundary(1, Bound(BOUND_NEUMANN, 0.0));      // xmin
  geom.set_boundary(2, Bound(BOUND_NEUMANN, 0.0));      // xmax
  geom.set_boundary(5, Bound(BOUND_DIRICHLET, Up));      // zmin
  geom.set_boundary(6, Bound(BOUND_DIRICHLET, 0.0 - VSource)); // zmax
#endif
    
  geom.set_boundary(7, Bound(BOUND_DIRICHLET,  VSource - VSource));
  geom.set_boundary(8, Bound(BOUND_DIRICHLET,  VPuller - VSource));
  geom.set_boundary(9, Bound(BOUND_DIRICHLET,  VLens - VSource));
  geom.set_boundary(10, Bound(BOUND_DIRICHLET, VRepeller - VSource));
  geom.set_boundary(11, Bound(BOUND_DIRICHLET, 0.0 - VSource));

  // *** Build the mesh and calculate the electric field *** //
  geom.build_mesh();

#ifndef RZ_MODE
  geom.build_surface();
#endif
  
  EpotField epot(geom);
  MeshScalarField scharge(geom);
  MeshVectorField bfield;
  EpotEfield efield(epot);

#ifdef RZ_MODE
  field_extrpl_e extrapl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                               FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
                               FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE};
  InitialPlasma init_plasma(AXIS_X, 0.0e-3);
  ParticleDataBaseCyl pdb(geom);

  bool pmirror[6] = {false, false, true, false, false, false};
#else
  field_extrpl_e extrapl[6] = {FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                               FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                               FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE};
  InitialPlasma init_plasma(AXIS_Z, 0.0e-3);
  ParticleDataBase3D pdb(geom);

  bool pmirror[6] = {false, false, false, false, false, false};
#endif
  
  efield.set_extrapolation(extrapl);
  pdb.set_mirror(pmirror);
  
  EpotBiCGSTABSolver solver(geom);

  solver.set_initial_plasma(Up, &init_plasma);

  //pdb.set_polyint( true );
  //solver.solve(epot, scharge);
  
  for( int a = 0; a < 10; a++ ) {

    ibsimu.message(1) << "Major cycle " << a << "\n";
    ibsimu.message(1) << "-----------------------\n";

    if( a == 1 ) {
	double rhoe = pdb.get_rhosum();
	solver.set_pexp_plasma(-rhoe, Te, Up);
    }

    solver.solve(epot, scharge);
    
    if( solver.get_iter() == 0) {
      ibsimu.message(1) << "No iterations, breaking major cycle\n";
      break;
    }

    efield.recalculate();

    pdb.clear();

#ifdef RZ_MODE
    pdb.add_2d_beam_with_energy(NPart, J_p1p, 1.0, m_p1p_amu,
			        E0, 0.0, Ti,
				geom.origo(0),
				0,
				geom.origo(0),
				r_plasma);

    pdb.add_2d_beam_with_energy(NPart, J_h2p, 1.0, m_h2p_amu,
				E0, 0.0, Ti,
				geom.origo(0),
				0,
				geom.origo(0),
				r_plasma);
#else
    // *** Protons: *** //
    pdb.add_cylindrical_beam_with_energy(NPart, J_p1p, 1.0, m_p1p_amu,
					 E0, 0.0, Ti, 
    				         Vec3D(0, 0, geom.origo(2)), // center
					 Vec3D(1, 0, 0), // dir1
					 Vec3D(0, 1, 0), // dir2,
				         r_plasma);
    // *** H2+ *** //
    pdb.add_cylindrical_beam_with_energy(NPart, J_h2p, 1.0, m_h2p_amu,
					 E0, 0.0, Ti, 
    				         Vec3D(0, 0, geom.origo(2)), // center
					 Vec3D(1, 0, 0), // dir1
					 Vec3D(0, 1, 0), // dir2,
				         r_plasma);
#endif
    
    pdb.iterate_trajectories(scharge, efield, bfield);

    if (fe > 0.0) {
      // space charge compensation at z > 125 mm
#ifdef RZ_MODE
      for(uint32_t i = 0; i < geom.size(0); i++ ) {
        double x = geom.origo(0) + geom.h() * i;
        if( x > scc_onset ) {
          for( uint32_t j = 0; j < geom.size(1); j++ ) {
            scharge(i,j) *= (1.0 - fe);
          }
        }
      }
#else
      for( uint32_t k = 0; k < geom.size(2); k++ ) {
        double z = geom.origo(2) + geom.h() * k;
        if( z > scc_onset ) {
    	  for( uint32_t i = 0; i < geom.size(0); i++ ) {
	    for( uint32_t j = 0; j < geom.size(1); j++ )
	      scharge(i,j,k) *= (1.0 - fe);
	  }
        }
      }    
#endif
    }

    /*
    TrajectoryDiagnosticData tdata;
    vector<trajectory_diagnostic_e> diag;
    diag.push_back( DIAG_X );
    diag.push_back( DIAG_XP );
    diag.push_back( DIAG_Y );
    diag.push_back( DIAG_YP );
    pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diag );
    Emittance emit_xxp( tdata(0).data(), tdata(1).data() );
    Emittance emit_yyp( tdata(2).data(), tdata(3).data() );
    */
  }

    
  if (interactive == true) {
    // *** Final PLotting in interactive mode *** //
    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_efield( &efield );
    plotter.set_bfield( &bfield );
    plotter.set_scharge( &scharge );
    plotter.new_geometry_plot_window();
    plotter.run();
  }
  // *** Plot a side view picture of the beam *** //
  GeomPlotter geomplotter(geom);
  geomplotter.set_size(4000,1000);
  geomplotter.set_font_size(20);
  geomplotter.set_epot(&epot);
  geomplotter.set_particle_database(&pdb);
#ifndef RZ_MODE
  string outfname = "MIST1_3D_VSource_" + to_string(VSource) + "V_VLens_" + to_string(VLens) + "V_fe_" + to_string(fe) + ".png";
  geomplotter.set_view(VIEW_ZY);
#else
  string outfname = "MIST1_RZ_VSource_" + to_string(VSource) + "V_VLens_" + to_string(VLens) + "V_fe_" + to_string(fe) + ".png";
#endif
  geomplotter.plot_png(outfname);

#ifdef RZ_MODE
  // *** Append VLens, Size and emittance to a text file *** //
  TrajectoryDiagnosticData tdata;
  vector<trajectory_diagnostic_e> diag;
  diag.push_back(DIAG_R);
  diag.push_back(DIAG_RP);
  diag.push_back(DIAG_AP);
  diag.push_back(DIAG_CURR);
  pdb.trajectories_at_plane(tdata, AXIS_X, geom.max(0), diag);
  EmittanceConv emit(100, 100,
		     tdata(0).data(), tdata(1).data(),
		     tdata(2).data(), tdata(3).data());

  ofstream dataout("variation_summary.csv", ios_base::app);
  dataout << VSource << ", "
          << VLens << ", "
	  << fe << ", "
          << sqrt(emit.beta() * emit.epsilon()) << ", "
          << emit.alpha() << ", "
	  << emit.beta() << ", "
	  << emit.epsilon() << "\n";
  dataout.close();
#endif
  
  // *** Save the particle distribution at the end of the simulation *** //
  if (savepart == true) {
    string poutfname = "MIST1_3D_VSource_" + to_string(VSource) + "V_VLens_" + to_string(VLens) + "V_fe_" + to_string(fe) + ".txt";
    ofstream fileOut(poutfname.c_str());
    for ( size_t k = 0; k < pdb.size(); k++) {
#ifdef RZ_MODE
      ParticleCyl &pp = pdb.particle(k);
      if (pp(1) < 0.99 * zmax)  // Omit particles that didn't make it (within a margin)
        continue;
      fileOut << setw(12) << pp.IQ() << " ";
      fileOut << setw(12) << pp.m() << " ";
      for (size_t j = 0; j < 6; j++)
#else
      Particle3D &pp = pdb.particle(k);
      if (pp(5) < 0.99 * zmax)  // Omit particles that didn't make it (within a margin)
        continue;
      fileOut << setw(12) << pp.IQ() << " ";
      fileOut << setw(12) << pp.m() << " ";
      for (size_t j = 0; j < 7; j++)
#endif
	fileOut << setw(12) << pp(j) << " ";
      fileOut << "\n";
    }
  }
}

int main(int argc, char **argv)
{ 
  // *** Parse command line for -h --help parameter *** //
  for (int i = 1; i < argc; i++){
    string arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) { // Display usage
      cerr << "Usage: " << argv[0] << " <option(s)>\n"
           << "Options:\n" 
           << "\t-h, --help\t\tShow this help message\n"
	   << "\t-i, --interactive\tToggle on interactive plotting mode\n"
	   << "\t-o, --output\t\tToggle saving of particle distribution at end\n"
	   << "\t-c, --compensation\tspace charge compensation factor fe [0.0; 1.0)\n"
	   << "\t-s, --vsource\t\tSource voltage (V)\n"
	   << "\t-p, --vpuller\t\tPuller voltage (V)\n"
	   << "\t-l, --vlens\t\tEinzel lens voltage (V)\n"
	   << "\t-r, --vrepeller\t\tRepeller voltage (V)\n"
	   << endl;
      return 0;
    }
  }

  // *** Parse other command line parameters ***//
  if (argc > 1) cout << "Command line parameters:" << endl;

  for (int i = 1; i < argc; i ++){
    string arg = argv[i];
    if ((arg == "-i") || (arg == "--interactive")) {
      interactive = true;
      cout << "\t--interactive" << endl;
    } else if ((arg == "-o") || (arg == "--output")) {
      savepart = true;
      cout << "\t--output" << endl;
    } else if ((arg == "-c") || (arg == "--compensation")) {
      if (i + 1 < argc) {  // make sure there is a value for the argument
	fe = atof(argv[++i]) / 100.0;
	if ((fe < 0) || (fe >= 100)) {
	  cerr << "--compensation must be in [0; 100)" << endl;
	  return 1;
	}
        cout << "\t--compensation = " << 100.0 * fe << endl;
      } else {
        cerr << "--compensation option requires one argument [0; 100)" << endl;
        return 1;
      }
    } else if ((arg == "-s") || (arg == "--vsource")){
      if (i + 1 < argc) {  // make sure there is a value for the argument
	VSource = atof(argv[++i]);
        cout << "\t--vsource = " << VSource << " V" << endl;
      } else {
        cerr << "--vsource option requires one argument" << endl;
        return 1;
      }
    } else if ((arg == "-p") || (arg == "--vpuller")){
      if (i + 1 < argc) {  // make sure there is a value for the argument
	VPuller = atof(argv[++i]);
        cout << "\t--vpuller = " << VPuller << " V" << endl;
      } else {
        cerr << "--vpuller option requires one argument" << endl;
        return 1;
      }
    } else if ((arg == "-l") || (arg == "--vlens")){
      if (i + 1 < argc) {  // make sure there is a value for the argument
	VLens = atof(argv[++i]);
        cout << "\t--vlens = " << VLens << " V" << endl;
      } else {
        cerr << "--vlens option requires one argument" << endl;
        return 1;
      }
    } else if ((arg == "-r") || (arg == "--vrepeller")){
      if (i + 1 < argc) {  // make sure there is a value for the argument
	VRepeller = atof(argv[++i]);
        cout << "\t--vrepeller = " << VRepeller << " V" << endl;
      } else {
        cerr << "--vrepeller option requires one argument" << endl;
        return 1;
      }
    }
  }
  
  try {
    ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
    ibsimu.set_thread_count( 4 );
    simu(&argc, &argv);
  } catch( Error e ) {
    e.print_error_message( ibsimu.message(0) );
  }

  return( 0 );
}
