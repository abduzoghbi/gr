/*
 * gr.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#include "inc/gr.hpp"

int main() {
	/*
	double			pos[4] = { 0 , 10 , 1e-5 , 0 };
	double			drdt[3] = {0,.0,1e-4};
	gr::photon 		ph( 0.7 , 0 , 0 , .8 , -1 , -1 );
	//gr::tetrad		tet( pos , drdt , 0.8 );
	ph.propagate( pos , drdt );
	*/

	/*
	double			pos[4] = { 0 , 20 , 1e-3 , 0 };
	double			drdt[3] = {0,.0,1e-4};
	gr::flash		fl( pos , drdt , 0.9 );
	fl.illum(4000);
	*/

	/*
	//gr::image		im( 0.9 , 1. , 40. , 100 );
	//im.write_hdf5("image.h5");
	gr::image		im( "image.h5" );
	*/



	double		rlim[2] = {1.2,100};
	gr::disk	disk( "illum.h5" , 200  , rlim , 20  );
	//disk.emissivity();

	double	tLim[2]		=	{ 0 , 100 };
	double	enLim[2]	=	{ 0.5 , 2 };
	disk.tf(400,tLim,400,enLim);


}
