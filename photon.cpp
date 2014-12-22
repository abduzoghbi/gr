/*
 * photon.cpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#include "inc/photon.hpp"

namespace gr {

photon::photon(double e,double l,double q,double spin,int rs, int ths) {
	E		=	e;
	L		=	l;
	Q		=	q;
	a		=	spin;
	rsign	=	rs;
	thsign	=	ths;
	rh_stop		=	false;
	disk_stop	=	false;
	inf_stop	=	false;

}

photon::~photon() {}



/******* Calculate rdot vector *********/
void photon::calc_rdot( const double* rvec , double* rdot ){
	double		rdot2,thdot2;

	// ---- Setup vars ---- //
	double		r,r2,sin2,cos2,D,S,DS,A,a2=a*a;
	r		=	rvec[1]; r2 = r*r;
	sin2	=	sin( rvec[2] ); sin2 = sin2*sin2;
	cos2 	= 	1-sin2;
	D		=	r2 - 2*r + a2;
	S		=	r2 + a2*cos2;
	DS		=	D*S;
	A		=	(r2+a2)*(r2+a2) - a2*D*sin2;
	// ------------------- //

	rdot[0]		=	(A*E - 2*a*r*L) / DS;
	rdot[3]		=	( 2*a*r*E + L*(D/sin2 - a2) ) / DS;

	thdot2		=	( Q - cos2*(a2*(-E*E) + L*L/sin2) ) / (S*S);
	rdot2		=	(D/S) * ( E*rdot[0] - L*rdot[3] - S*thdot2 );


	if( fabs(thdot2)<(1e-5*RDOT_ZERO) ) thdot2 = ((thdot2<0)?-1:1)*(1e-5*RDOT_ZERO);

	if( rdot2 < 0 ){
		if(fabs(rdot2)<RDOT_ZERO and ir_change>NCHANGE){ rsign *= -1; ir_change=0; }
		rdot[0]	=	-1; // will be caught outside
		return;
	}
	if( thdot2 < 0 ){
		if(fabs(thdot2)<RDOT_ZERO and ith_change>NCHANGE){thsign *= -1; ith_change=0;}
		rdot[0]	=	-1; // will be caught outside
		return;
	}

	rdot[1]		=	rsign  * sqrt( rdot2  );
	rdot[2]		=	thsign * sqrt( thdot2 );


	//printf("--- %3.3e %3.3e %3.3e %3.3e %d %d\n",rdot[0],rdot[1],rdot[2],rdot[3],rsign,thsign);
}
/* =================================== */


/******* Propagate a photon from a position pos *********/
void photon::propagate( double* rvec, double* rdot ) {

	int			iter,i,status;
	double		diff,err,tau,tau_c,dtau = 0.001,halfpi=M_PI/2.,twopi=M_PI*2,rh;
	double		*rvec_c,*rvec_tmp,*rvec_err;
	bool		stopping;

	gsl_odeiv2_system				sys = { int_func , NULL , 4, (void*)(this) };
	const gsl_odeiv2_step_type 		*step = gsl_odeiv2_step_rkf45;
	gsl_odeiv2_step					*stepper = gsl_odeiv2_step_alloc ( step , 4 );


	/* INITIALIZE SOME VARIALES */
	rh			=	1 + sqrt(1-a*a);
	rvec_c		=	new double[4];
	rvec_tmp	=	new double[4];
	rvec_err	=	new double[4];
	for( i=0 ; i<4 ; i++ ){ rvec_c[i] = rvec[i]; rvec_tmp[i] = rvec[i];}
	tau			=	0;
	tau_c		=	0;
	ir_change	=	0;
	ith_change	=	0;
	stopping	=	false;

	/* ------------------------ */



	/* START THE INTEGRATION LOOP */
	for ( iter=1 ; iter<= NLOOP ; iter++ ){

		status	=	gsl_odeiv2_step_apply ( stepper , tau , dtau, rvec_tmp, rvec_err, NULL, NULL , &sys);
		if( status ){ // something is wrong
			tau		=	tau_c;
			for( i=0 ; i<4 ; i++){ rvec_tmp[i] = rvec_c[i]; }
			dtau	/=	DTAU_FAC;
			continue;
		}
		tau		+=	dtau;


		// Get maximum integration error and change dtau accordingly //
		err		=	-1e6;
		for( i=0 ; i<4 ; i++ ){ if(err<fabs(rvec_err[i])) err = fabs(rvec_err[i]);}
		if( err>INTG_ERROR ){
			dtau	/=	DTAU_FAC;
		}else {
			dtau	*=	DTAU_FAC;
		}


		/* If at RINF, stop */
		if ( rvec_tmp[1] > RINF ){
			diff		=	rvec_tmp[1] - RINF;
			tau			-=	(diff * (tau-tau_c)/(rvec_tmp[1] - rvec_c[1]));
			rvec_tmp[0]	-=	diff;
			rvec_tmp[1]	=	RINF;
			stopping	=	true;
			inf_stop	=	true;
		}

		/* If at the Horizon, stop */
		if ( (rvec_tmp[1]-rh) < STP_ZERO ) {stopping	=	true; rh_stop	=	true;}


		/* Are we at the disk at halfpi ? */
		if( fabs( rvec_tmp[2] - halfpi ) < STP_ZERO ) { stopping = true; disk_stop = true;}

		// Or are we close to HALF_PI? //
		if( ( rvec_tmp[2] - halfpi ) >= STP_ZERO and not stopping ){
			tau		=	tau_c;
			dtau	/=	(DTAU_FAC*DTAU_FAC);
			for( i=0 ; i<4 ; i++){ rvec_tmp[i] = rvec_c[i]; }
			continue;
		}

		// is theta negative? //
		if ( rvec_tmp[2] < 0 ) {
			rvec_tmp[2]		*=	-1;
			rvec_tmp[3]		+=	M_PI;
			thsign *= -1; ith_change=0;
		}


		/* If we get here, things are ok. Update saved rvec_c and tau_c */
		for( i=0 ;i<4 ; i++ ) { rvec_c[i] = rvec_tmp[i]; rvec[i] = rvec_tmp[i];}
		tau_c		=	tau;
		ir_change++; ith_change++;


		//_print_xyz(rvec);
		//calc_rdot( rvec , rdot );_const_of_motion(rvec,rdot);
		//for( i=0 ;i<4 ; i++ ) { printf("%3.3e ",rvec[i]);} printf("%3.3e %3.3e %3d\n",tau,dtau,iter);

		if ( stopping ) {
			calc_rdot( rvec , rdot );
			if( rdot[0] != -1 ){
				while( rvec[3] > twopi ) rvec[3] -= twopi;
				while( rvec[3] < 0 ) rvec[3] += twopi;
			}else{
				diff = sqrt(-1);
				rdot[0] = diff;rdot[1] = diff;rdot[2] = diff; rdot[3] = diff;
			}
			break;
		}
	}
	/* -------------------------- */


	/* CLEAR ALLOCATED MEMORY */
	delete[] rvec_c;
	delete[] rvec_tmp;
	delete[] rvec_err;
	gsl_odeiv2_step_free( stepper );
	/* ----------------------- */


}
/* ==================================================== */


/******** check constants of motion for debugging **********/
void photon::_const_of_motion(const double* rvec , double* rdot ){
	double	e,l,q,int1;
	calc_rdot( rvec , rdot );
	double		r,r2,sin2,cos2,D,S,a2=a*a;
	r		=	rvec[1]; r2 = r*r;
	sin2	=	sin( rvec[2] ); sin2 = sin2*sin2;
	cos2 	= 	1-sin2;
	D		=	r2 - 2*r + a2;
	S		=	r2 + a2*cos2;

	e		=	(1-2*r/S)*rdot[0] + (2*a*r*sin2/S)*rdot[3];
	l		=	(-2*a*r*sin2/S)*rdot[0] + (r2+a2+2*a2*r*sin2/S)*sin2*rdot[3];
	q		=	(S*S*rdot[2]*rdot[2]) + cos2*( a2*(-e*e) + l*l/sin2);
	int1	=	-(1-2*r/S)*rdot[0]*rdot[0] -(4*a*r*sin2/S)*rdot[0]*rdot[3] + \
			(S/D)*rdot[1]*rdot[1] + S*rdot[2]*rdot[2] + (r2+a2+2*a2*r*sin2/S)*sin2*rdot[3]*rdot[3];


    printf("%3.3e %3.3e %3.3e %3.3e\n",e,l,q,int1);
}
// ========================================================= /


/******* INTEGRATION FUNCTION ********/
int photon::int_func (double t, const double* r, double* drdt, void *params) {
	photon	*ph = &(*(photon *)params);
	ph->calc_rdot( r , drdt );
	if( drdt[0] == -1 ) return -1; // something is wrong
	return GSL_SUCCESS;
}
/* ================================= */


/******** Print coordinates in x,t,z for debugging **********/
void photon::_print_xyz( const double* rvec ){
	double	x,y,z,r,th,phi;
	r = rvec[1]; th=rvec[2]; phi=rvec[3];
    x		=	r*cos(phi)*sin(th);
    y		=	r*sin(phi)*sin(th);
    z		=	r*cos(th);
    printf("%3.3e %3.3e %3.3e\n",x,y,z);
}
// ========================================================= /


} /* namespace gr */
