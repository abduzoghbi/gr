/*
 * flash.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#include "inc/flash.hpp"

namespace gr {

flash::flash( double s[] , double drdt[] , double spin ) {
	a		=	spin;
	src		=	new double[4];
	for( int i=0 ; i<4 ; i++ ){ src[i] = s[i]; }

	tet		=	new tetrad( s , drdt , spin );
	consts	=	new double[4];
	sign	=	new int[2];
	rvec	=	new double[4];
	rdot	=	new double[4];
}

flash::~flash() {
	delete[] src;
	delete[] consts;
	delete[] sign;
	delete[] rvec;
	delete[] rdot;
	delete tet;
}

/******** Const of motion for a given alpha,beta in source frame **********/
void flash::const_of_motion ( double alpha, double beta ){
	double		pt[4],rd[4];
	int			i,j;
	double		dum;

	/* Momentum in emitting source frame */
	pt[0]		=	1.0;
	pt[1]		=	cos(alpha);
	pt[2]		=	sin(alpha)*cos(beta);
	pt[3]		=	sin(alpha)*sin(beta);

	/* Momentum in coordinate frame */
	for( i=0 ; i<4 ; i++ ){
		dum		=	0;
		for( j=0 ; j<4 ;j++ ){ dum +=	tet->operator()(j,i) * pt[j];}
		rd[i]	=	dum;
	}

	/* Now get constants of motion given rdot */
	double		r,S,sin2,cos2,E,L,Q;
	r		=	src[1];
	sin2	=	sin(src[2]);sin2 = sin2*sin2;
	cos2	=	1. - sin2;
	S		=	r*r + a*a*cos2;

	E		=	(1-2*r/S) * rd[0] + (2*a*r*sin2/S) * rd[3];
	L		=	(-2*a*r*sin2/S) * rd[0] + (r*r+a*a+2*a*a*r*sin2/S)*sin2 * rd[3];
	Q		=	S*S*rd[2]*rd[2] + cos2*(a*a*(-E*E) + L*L/sin2);

	consts[0]	=	E;
	consts[1]	=	L;
	consts[2]	=	Q;

    sign[0]		=	(rd[1]<0)?-1:1;
    sign[1]		=	(rd[2]<0)?-1:1;

}
/* ====================================================================== */


/******** illuminate the source with a photon at given alpha, beta **********/
void flash::illum ( double alpha , double beta ){
	const_of_motion( alpha , beta );
	photon	ph( consts[0] , consts[1] , consts[2] , a , sign[0] , sign[1] );
	ph.propagate( src );
	tau		=	ph.Tau;
	rvec	=	ph.rvec;
	rdot	=	ph.rdot;
}
/* ======================================================================== */


/******** flash the source with isotropic num of photons **********/
void flash::illum( int num ){
	int		i,j;
	double	alpha,beta;

	for( i=0 ; i<num ; i++ ){
		alpha		=	acos( i*2.0/(num-1) - 1);
		for( j=0 ; j<num ; j++ ){
			beta	=	2*M_PI*j/(num-1);
			illum( alpha , beta );
		}
	}

}
/* ============================================================== */
} /* namespace gr */
