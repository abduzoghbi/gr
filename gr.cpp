/*
 * gr.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#include "inc/gr.hpp"

int main() {
	//gr::photon 		ph( 0.7 , 0 , 0 , .8 , -1 , -1 );
	double			pos[4] = { 0 , 10 , 1.5 , 0 };
	double			drdt[3] = {.3,.0,.05};
	gr::tetrad		tet( pos , drdt , 0.8 );
	//ph.propagate( pos );
}
