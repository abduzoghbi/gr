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

	// ISCO and Horizon radii //
	double		z1,z2;
	z1		=	(1 + pow(1-a*a,1./3)*( pow(1+a,1/3) + pow(1-a,1/3) ));
	z2		=	sqrt( 3*a*a + z1*z1 );
	rms		=	3 + z2	- sqrt( (3-z1)*(3+z1+2*z2) );

}

image::~image() {
	delete[] rvec;
}


/******* Project image to disk using npix in each side *********/
void image::project_image( int npix ){
	int		i,j,nside;
	double	unit,x,y, out[4];

	/* npix should be odd to have symmetry */
	if( npix%2 == 0) npix++;
	nside	=	(npix-1)/2;

	unit	=	imsize/npix;

	for( i=-nside ; i<=nside ; i++ ){
		x	=	unit * i;
		for( j=-nside ; j<=nside ; j++ ){
			y	=	unit * j;
			proj_xy( x , y , out );
		}
	}
}
/*=============================================================*/


/******** Project a photon from x,y and get r,phi,g,t **********/
void image::proj_xy( double x , double y , double* out ){
	double			E,L,Q,vel[4],g;
	int				thsign;

	if(fabs(x)<1e-6) { x = (x<0)?-1e-6:1e-6;}
	if(fabs(y)<1e-6) { y = (y<0)?-1e-6:1e-6;}


	E		=	1.0;
	L		=	x*sin0;
	Q		=	y*y + cos02*( -a*a*E*E + L*L/sin02);
	thsign	=	(y<0)?1:-1;
	photon		ph( E , L , Q , a , -1 , thsign );
	ph.propagate( rvec );
	disk::disk_velocity( ph.rvec[1] , vel , a , rms );
	g		=	-disk::p_dot_v( ph.rvec[1] , ph.rvec[2] , ph.rdot , vel , a );
	out[0]	=	ph.rvec[0];		// t
	out[1]	=	ph.rvec[1];		// r
	out[2]	=	ph.rvec[3];		// phi
	out[3]	=	g;				// g
}
/*=============================================================*/

} /* namespace gr */
