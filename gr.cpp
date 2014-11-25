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

	double			pos[4] = { 0 , 3 , 0.5 , 0 };
	double			drdt[3] = {.3,.0,.05};
	gr::flash		fl( pos , drdt , 0.8 );
	fl.illum(4);

	/*
	gr::image		im( 0.8 , .8 , 20. );
	im.project_image(40);
	*/
}
