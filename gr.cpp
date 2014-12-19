/*
 * gr.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#include "inc/gr.hpp"

int main() {
	/*
	gr::photon 		ph( 0.7 , 0 , 0 , .8 , -1 , -1 );
	gr::tetrad		tet( pos , drdt , 0.8 );
	ph.propagate( pos );
	*/
	/*
	double			pos[4] = { 0 , 10 , 1e-5 , 0 };
	double			drdt[3] = {0,.0,1e-4};
	gr::flash		fl( pos , drdt , 1 );
	fl.illum(1000);
	*/


	double		rlim[2] = {1.2,50};
	gr::disk	disk( "illum.h5" , 300  , rlim , 300  );
	//disk.emissivity();
	disk.tf();


	/*
	gr::image		im( 1 , 1.4 , 40. );
	im.project_image(200);
	*/
}
