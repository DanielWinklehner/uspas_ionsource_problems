#include <fstream>
#include <iomanip>
#include <limits>
#include "meshvectorfield.hpp"
#include "dxf_solid.hpp"
#include "mydxffile.hpp"
#include "gtkplotter.hpp"
#include "geomplotter.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "meshcolormap.hpp"
#include "epot_efield.hpp"
#include "error.hpp"
#include "ibsimu.hpp"
#include "trajectorydiagnostics.hpp"
#include "particledatabase.hpp"
#include "particlediagplotter.hpp"



using namespace std;



void simu( int argc, char **argv )
{
    std::ifstream is_geom( argv[1] );
    if( !is_geom.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[1] + "\'" ) );
    Geometry geom( is_geom );
    is_geom.close();

    std::ifstream is_epot( argv[2] );
    if( !is_epot.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[2] + "\'" ) );
    EpotField epot( is_epot, geom );
    is_epot.close();

    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
				     FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    std::ifstream is_pdb( argv[3] );
    if( !is_pdb.good() )
	throw( Error( ERROR_LOCATION, (string)"couldn\'t open file \'" + argv[3] + "\'" ) );
    ParticleDataBaseCyl pdb( is_pdb, geom );
    is_pdb.close();

    // ------------------------------------------
    // Produce particle diagnostic data for emittance calculation and plotting for each charge
    TrajectoryDiagnosticData tdata;
    vector<trajectory_diagnostic_e> diag;
    diag.push_back( DIAG_R );
    diag.push_back( DIAG_RP );
    diag.push_back( DIAG_AP );
    diag.push_back( DIAG_CURR );
    diag.push_back( DIAG_CHARGE );
    pdb.trajectories_at_plane( tdata, AXIS_X, geom.max(0), diag );

    // Go through the charge states
    for( int a = 1; a <= 7; a++ ) {
	// Produce data vectors, which only contain data for charge state a.
	std::vector<double> r;
	std::vector<double> rp;
	std::vector<double> ap;
	std::vector<double> curr;
	for( uint32_t b = 0; b < tdata.traj_size(); b++ ) {
	    if( fabs(tdata(b,4)-a) < 0.1 ) {
		r.push_back( tdata(b,0) );
		rp.push_back( tdata(b,1) );
		ap.push_back( tdata(b,2) );
		curr.push_back( tdata(b,3) );
	    }
	}

	// Calculate and print emittance
	EmittanceConv emit( 100, 100, r, rp, ap, curr );
	ibsimu.message(1) << a << " "
			  << emit.epsilon() << "\n";

	// Plot emittance histogram using plotting tools inside IBSimu
	const Histogram2D &histo = emit.histogram();
	double range[4];
	histo.get_range( range );
	MeshColormap colormap( range, histo.n(), histo.m(), histo.get_data() );
	colormap.set_interpolation( INTERPOLATION_CLOSEST );
	std::vector<Palette::Entry> entries;
	entries.push_back( Palette::Entry(Vec3D(1,1,1),0) );
	entries.push_back( Palette::Entry(Vec3D(1,1,0),1) );
	entries.push_back( Palette::Entry(Vec3D(1,0,0),2) );
	entries.push_back( Palette::Entry(Vec3D(0,0,0),3) );
	Palette palette( entries );
	palette.normalize();
	colormap.set_palette( palette );

	Frame frame;
	frame.add_graph( PLOT_AXIS_X1, PLOT_AXIS_Y1, &colormap );
	cairo_surface_t *surface = cairo_image_surface_create( CAIRO_FORMAT_ARGB32,
							       1024, 768 );
	cairo_t *cairo = cairo_create( surface );
	frame.set_geometry( 1024, 768, 0, 0 );
	frame.draw( cairo );
	stringstream fn;
	fn << "emittance_" << a << ".png";
	cairo_surface_write_to_png( surface, fn.str().c_str() );
	cairo_destroy( cairo );
	cairo_surface_destroy( surface );
    }


    MeshScalarField tdens( geom );
    pdb.build_trajectory_density_field( tdens );

    GTKPlotter plotter( &argc, &argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_efield( &efield );
    plotter.set_trajdens( &tdens );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();
}


int main( int argc, char **argv )
{
    if( argc <= 3 ) {
	cerr << "Usage: analysis geom epot pdb\n";
	exit( 1 );
    }

    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
	simu( argc, argv );
    } catch( Error e ) {
	e.print_error_message( ibsimu.message( 0 ) );
        exit( 1 );
    }

    return( 0 );
}
