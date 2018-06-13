#include <fstream>
#include <iostream>
#include <ibsimu.hpp>
#include <error.hpp>
#include <geometry.hpp>
#include <dxf_solid.hpp>
#include <mydxffile.hpp>
#include <epot_field.hpp>
#include <epot_efield.hpp>
#include <meshvectorfield.hpp>
#include <epot_bicgstabsolver.hpp>
#include <gtkplotter.hpp>
#include <trajectorydiagnostics.hpp>



using namespace std;


double Vplasma = 0;
double Vpuller = -15e3;
double Veinzel_1 = -24e3;
double Veinzel_2 = -24e3;
double Vgnd = -14e3;

int Nrounds = 50;
double h = 0.5e-3;
double nperh = 500;
double r0 = 10e-3;
double Npart = r0/h*nperh;
double Jtotal = 30.0;
double Te = 10.0;
double E0 = 10.0;
double Tt = 1.0;
double Up = 20.0;
double sc_alpha = 0.7;


void add_beam( Geometry &geom, ParticleDataBaseCyl &pdb, double q, double m, 
	       double Jtotal, double frac )
{
    pdb.add_2d_beam_with_energy( Npart*frac, Jtotal*frac, q, m,
				 E0, 0.0, Tt,
				 geom.origo(0), 0,
				 geom.origo(0), r0 );
}


void simu( int *argc, char ***argv )
{
    double q = 6;
    double m = 15;
    double B0 = 0.9;
    double r_aperture = 4.0e-3;
    double vz = sqrt(-2.0*q*CHARGE_E*Vgnd/(m*MASS_U));
    double Erms = q*CHARGE_E*B0*r_aperture*r_aperture/(8*m*MASS_U*vz);
    ibsimu.message(1) << "Erms = "<< Erms << " m rad\n";

    Vec3D origin( -5e-3, 0, 0 );
    Vec3D sizereq( 505e-3, 65e-3, 0 );
    Int3D size( floor(sizereq[0]/h)+1,
		floor(sizereq[1]/h)+1,
		1 );
    Geometry geom( MODE_CYL, size, origin, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 1 );
    dxffile->read( "geom.dxf" );

    DXFSolid *s1 = new DXFSolid( dxffile, "7" );
    s1->scale( 1e-3 );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "8" );
    s2->scale( 1e-3 );
    geom.set_solid( 8, s2 );
    DXFSolid *s3 = new DXFSolid( dxffile, "9" );
    s3->scale( 1e-3 );
    geom.set_solid( 9, s3 );
    DXFSolid *s4 = new DXFSolid( dxffile, "10" );
    s4->scale( 1e-3 );
    geom.set_solid( 10, s4 );
    DXFSolid *s5 = new DXFSolid( dxffile, "11" );
    s5->scale( 1e-3 );
    geom.set_solid( 11, s5 );
    DXFSolid *s6 = new DXFSolid( dxffile, "12" );
    s6->scale( 1e-3 );
    geom.set_solid( 12, s6 );
    DXFSolid *s7 = new DXFSolid( dxffile, "13" );
    s7->scale( 1e-3 );
    geom.set_solid( 13, s7 );
    DXFSolid *s8 = new DXFSolid( dxffile, "14" );
    s8->scale( 1e-3 );
    geom.set_solid( 14, s8 );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) ); // xmin
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) ); // xmax
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) ); // rmin
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) ); // rmax

    geom.set_boundary(  7, Bound(BOUND_DIRICHLET, Vplasma) );
    geom.set_boundary(  8, Bound(BOUND_DIRICHLET, Vpuller) );
    geom.set_boundary(  9, Bound(BOUND_DIRICHLET, Veinzel_1) );
    geom.set_boundary( 10, Bound(BOUND_DIRICHLET, Vgnd) );
    geom.set_boundary( 11, Bound(BOUND_DIRICHLET, Veinzel_2) );
    geom.set_boundary( 12, Bound(BOUND_DIRICHLET, Vgnd) );
    geom.set_boundary( 13, Bound(BOUND_DIRICHLET, Vgnd) );
    geom.set_boundary( 14, Bound(BOUND_DIRICHLET, Vgnd) );

    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshScalarField scharge_ave( geom );

    //ASSIGNMENT: IBSimu Task 1: Uncomment the next line and comment out the following
    //three lines to load an empty bfield
    //MeshVectorField bfield;

    bool fsel[3] = {true, true, false};
    MeshVectorField bfield( MODE_CYL, fsel, 1.0e-3, 1.0, "bfield.txt" );
    bfield.translate( Vec3D(-19.3e-3,0,0) );

    EpotEfield efield( epot );
    field_extrpl_e extrapl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
				  FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
				  FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( extrapl );
  
    EpotBiCGSTABSolver solver( geom );
    InitialPlasma init_plasma( AXIS_X, 0.5e-3 );
    solver.set_initial_plasma( Up, &init_plasma );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = {false, false,
		       true, false,
		       false, false};
    pdb.set_mirror( pmirror );

    for( int a = 0; a < Nrounds; a++ ) {

	ibsimu.message(1) << "Major cycle " << a << "\n";
	ibsimu.message(1) << "-----------------------\n";

	if( a == 1 ) {
	    double rhoe = pdb.get_rhosum();
	    solver.set_pexp_plasma( rhoe, Te, Up );
	}

	solver.solve( epot, scharge_ave );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
      
	efield.recalculate();

	pdb.clear();
	add_beam( geom, pdb, 1, 15, Jtotal, 0.050 );
	add_beam( geom, pdb, 2, 15, Jtotal, 0.100 );
	add_beam( geom, pdb, 3, 15, Jtotal, 0.200 );
	add_beam( geom, pdb, 4, 15, Jtotal, 0.310 );
	add_beam( geom, pdb, 5, 15, Jtotal, 0.250 );
	add_beam( geom, pdb, 6, 15, Jtotal, 0.085 );
	add_beam( geom, pdb, 7, 15, Jtotal, 0.005 );
	//add_beam( geom, pdb, 6, 15, Jtotal, 1.0 );
	pdb.iterate_trajectories( scharge, efield, bfield );

    
	TrajectoryDiagnosticData tdata;
	vector<trajectory_diagnostic_e> diag;
	diag.push_back( DIAG_R );
	diag.push_back( DIAG_RP );
	diag.push_back( DIAG_AP );
	diag.push_back( DIAG_CURR );
	pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0), diag );
	EmittanceConv emit( 100, 100, 
			    tdata(0).data(), tdata(1).data(), 
			    tdata(2).data(), tdata(3).data() );

	ofstream dataout( "emit.txt", ios_base::app );
	dataout << emit.alpha() << " "
		<< emit.beta() << " "
		<< emit.epsilon() << "\n";
	dataout.close();

	if( a == 0 ) {
	    scharge_ave = scharge;
	} else {
	    double sc_beta = 1.0-sc_alpha;
	    uint32_t nodecount = scharge.nodecount();
	    for( uint32_t b = 0; b < nodecount; b++ ) {
		scharge_ave(b) = sc_alpha*scharge(b) + sc_beta*scharge_ave(b);
	    }
	}
    }

    geom.save( "geom.dat" );
    epot.save( "epot.dat" );
    pdb.save( "pdb.dat" );

    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_particledatabase( &pdb );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_bfield( &bfield );
    plotter.set_scharge( &scharge );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    remove( "emit.txt" );

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 2 );
	simu( &argc, &argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message(0) );
    }

    return( 0 );
}
