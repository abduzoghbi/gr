/*
 * gr.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#include "inc/gr.hpp"

int main( int argc , char* argv[] ) {

	if( argc != 2 ) {
		gr::input::print_usage();
		exit(1);
	}
	gr::input		in( argv[1] );

	if( in.mode == 1 ){ // flash //
		gr::flash		flash( in.src_pos , in.src_drdt , in.spin );
		flash.illum( in.flash_num , in.flash_file );
	}

	if( in.mode == 2 ){ // image //
		gr::image		image( in.spin , in.theta_o , in.imsize , in.npix );
		image.write_hdf5( in.image_file );
	}

	if( in.mode == 3 ){ // Emissivity //
		gr::disk 		disk( in.flash_file , in.nr , in.rL , in.nphi );
		disk.emissivity(true);
	}

	if( in.mode == 4 ){ // image_flux //
		gr::disk 		disk( in.flash_file , in.nr , in.rL , in.nphi );
		disk.image_flux( in.image_file );
	}

	if( in.mode == 5 ){ // image_flux_time //
		gr::disk 		disk( in.flash_file , in.nr , in.rL , in.nphi );
		disk.image_flux_time( in.image_file , in.ntime , in.tL );
	}

	if( in.mode == 6 ){ // Transfer Function //
		gr::disk 		disk( in.flash_file , in.nr , in.rL , in.nphi );
		disk.tf( in.image_file , in.ntime , in.tL , in.nenergy , in.enL );
	}

}

