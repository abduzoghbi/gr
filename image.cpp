/*
 * image.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#include "inc/image.hpp"

namespace gr {

image::image( double spin , double th_o , double im_s ) {
	a		=	spin;
	imsize	=	im_s;

	rvec		=	new double[4];
	rvec[0]		=	0.0;
	rvec[1]		=	RINF;
	rvec[2]		=	th_o;
	rvec[3]		=	0.0;

	sin0	=	sin(th_o);
	cos02	=	cos(th_o)*cos(th_o);
	sin02	=	sin0*sin0;

}

image::~image() {
	delete[] rvec;
}


/******* Project image to disk using npix in each side *********/
void image::project_image( int npix ){
	int		i,j,nside;
	double	unit,x,y,E,L,Q,*vec;
	vec		=	new double[4];

	/* npix should be odd to have symmetry */
	if( npix%2 == 0) npix++;
	nside	=	(npix-1)/2;

	unit	=	imsize/npix;

	for( i=-nside ; i<=nside ; i++ ){
		for( j=-nside ; j<=nside ; j++ ){
			//i = -nside; j= -nside+1;
			x	=	unit * i;
			y	=	unit * j;

			if(fabs(x)<1e-6) { x = 1e-6;}
			if(fabs(y)<1e-6) { y = 1e-6;}

			E	=	1.0;
			L	=	x*sin0;
			Q	=	y*y + cos02*( -a*a*E*E + L*L/sin02);


			photon	ph( E, L , Q , a , -1 , (y<0)?1:-1 );
			vec	=	rvec;
			ph.propagate( vec );
			ph._print_xyz();
			//exit(0);

		}
	}
}
/*=============================================================*/


} /* namespace gr */
