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
	int		i,j,ii,jj,nside;
	double	unit,x,y,dx,dy,dr,dphi,dA, out[4],outp[4],outn[4],*flux;


	/* npix should be odd to have symmetry */
	if( npix%2 == 0) npix++;
	nside	=	(npix-1)/2;
	flux	=	new double[npix*npix];

	unit	=	imsize/npix;

	printf("xcent ");for( i=-nside;i<=nside;i++) {printf("%g ",unit*i);}printf("\n");
	printf("ycent ");for( i=-nside;i<=nside;i++) {printf("%g ",unit*i);}printf("\n");

	for( i=-nside,ii=0 ; i<=nside ; i++ ){
		x	=	unit * i;
		dx	=	fabs(unit*(i+1) - x);
		for( j=-nside,jj=0 ; j<=nside ; j++ ){
			y	=	unit * j;
			dy	=	fabs(unit*(j+1) - y);
			proj_xy( x , y , out );
			/*
			proj_xy( x-dx/2 , y-dy/2 , outn );
			proj_xy( x+dx/2 , y+dy/2 , outp );
			dr		=	fabs(outp[1] - outn[1]);
			dphi	=	fabs(outp[2] - outn[2]);
			*/
			dA		=	disk::proper_disk_area( out[1] , a , rms ) ;//* dphi * dr;
			flux[ii*npix + jj]	=	pow(out[3],3)/dA;

			printf("%g ", flux[ ii*npix + jj]);
			jj++;
		}
		printf("\n");
		ii++;
	}
	delete[] flux;
}
/*=============================================================*/


/******** Project a photon from x,y and get r,phi,g,t **********/
void image::proj_xy( double x , double y , double* out ){
	double			E,L,Q,vel[4],g,twopi=2*M_PI;;
	int				thsign;

	if(fabs(x)<1e-6) { x = (x<0)?-1e-6:1e-6;}
	if(fabs(y)<1e-6) { y = (y<0)?-1e-6:1e-6;}


	E		=	1.0;
	L		=	x*sin0;
	Q		=	y*y + cos02*( -a*a*E*E + L*L/sin02);
	thsign	=	(y<0)?1:-1;
	photon		ph( E , L , Q , a , -1 , thsign );
	ph.propagate( rvec );
	if( ph.rh_stop ){
		g	=	0;
	}else{
		disk::disk_velocity( ph.rvec[1] , vel , a , rms );
		g		=	-disk::p_dot_v( ph.rvec[1] , ph.rvec[2] , ph.rdot , vel , a );
	}
	out[0]	=	ph.rvec[0];		// t
	out[1]	=	ph.rvec[1];		// r
	out[2]	=	ph.rvec[3];		// phi
	out[3]	=	g;				// g
}
/*=============================================================*/



} /* namespace gr */
