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
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"


using namespace std;


double Veinzel = -15e3;
double E0 = 80e3;

double q = 8.0;
double m = 40.0;
double Tt = 0.01;
double r0 = 10e-3;
double I = 1e-3;
double J = I/(r0*r0*M_PI);


void simu( int *argc, char ***argv )
{
    double h = 1e-3;
    Vec3D origo( -35e-3, 
		 -35e-3, 
		 0 );
    double sizereq[3] = { 70e-3,
                          70e-3, 
                          170e-3 };
    Int3D meshsize( (int)floor(sizereq[0]/h)+1,
                    (int)floor(sizereq[1]/h)+1,
		    (int)floor(sizereq[2]/h)+1 );
    Geometry geom( MODE_3D, meshsize, origo, h );

    MyDXFFile *dxffile = new MyDXFFile;
    dxffile->set_warning_level( 2 );
    dxffile->read( "einzel3d.dxf" );
    DXFSolid *s1 = new DXFSolid( dxffile, "gnd" );
    s1->scale( 1e-3 );
    s1->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 7, s1 );
    DXFSolid *s2 = new DXFSolid( dxffile, "einzel" );
    s2->scale( 1e-3 );
    s2->define_2x3_mapping( DXFSolid::rotz );
    geom.set_solid( 8, s2 );

    geom.set_boundary( 1, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 2, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 5, Bound(BOUND_NEUMANN, 0.0) );
    geom.set_boundary( 6, Bound(BOUND_NEUMANN, 0.0) );

    geom.set_boundary( 7, Bound(BOUND_DIRICHLET, 0.0) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, Veinzel) );

    geom.build_mesh();
    geom.build_surface();

    EpotBiCGSTABSolver solver( geom );
    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE,
                                     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    ParticleDataBase3D pdb( geom );
    bool pmirror[6] = { false, false, false, false, false, false };
    pdb.set_mirror( pmirror );

    for( size_t i = 0; i < 10; i++ ) {

	ibsimu.message(1) << "Major cycle " << i << "\n";
	ibsimu.message(1) << "-----------------------\n";

	solver.solve( epot, scharge );
	if( solver.get_iter() == 0 ) {
	    ibsimu.message(1) << "No iterations, breaking major cycle\n";
	    break;
	}
	efield.recalculate();

        pdb.clear(); 
	pdb.add_cylindrical_beam_with_energy( 1000, J, q, m,
					      E0, 0, Tt,
					      Vec3D(0,0,0),
					      Vec3D(1,0,0),
					      Vec3D(0,1,0),
					      r0 );
	
        pdb.iterate_trajectories( scharge, efield, bfield );

	TrajectoryDiagnosticData tdata;
        std::vector<trajectory_diagnostic_e> diagnostics;
        diagnostics.push_back( DIAG_X );
        diagnostics.push_back( DIAG_XP );
        pdb.trajectories_at_plane( tdata, AXIS_Z, geom.max(2)-geom.h(), diagnostics );
        Emittance emit( tdata(0).data(), tdata(1).data() );

        // Output
        ofstream dout( "emit.txt", ios_base::app );
        dout << emit.alpha() << " "
             << emit.beta() << " "
             << emit.epsilon() << "\n";
        dout.close();

    }
    
    MeshScalarField tdens(geom);
    pdb.build_trajectory_density_field(tdens);

    GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_bfield( &bfield );
    plotter.set_trajdens( &tdens );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    remove( "emit.txt" );

    try {
	//ibsimu.set_message_output( "ibsimu.txt" );
        ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( &argc, &argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
