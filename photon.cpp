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
	a2		=	a*a;
	rh		=	1 + sqrt(1-a2);

	rvec	=	new double[4];
	rdot	=	new double[4];
}

photon::~photon() {}


/******* Useful Variables *********/
void photon::_setup_vars(){
	r		=	rvec[1]; r2 = r*r;
	sin1	=	sin( rvec[2] ); sin2 = sin1*sin1;
	cos1	=	cos( rvec[2] ); cos2 = cos1*cos1;
	D		=	r2 - 2*r + a2;
	S		=	r2 + a2*cos2;
	DS		=	D*S;
	A		=	(r2+a2)*(r2+a2) - a2*D*sin2;
}
/* =============================== */


/******* Calculate rdot vector *********/
void photon::calc_rdot(){
	double		rdot2,thdot2;
	_setup_vars();

	rdot[0]		=	(A*E - 2*a*r*L) / DS;
	rdot[3]		=	( 2*a*r*E + L*(D/sin2 - a2) ) / DS;

	thdot2		=	( Q - cos2*(a2*(-E*E) + L*L/sin2) ) / (S*S);
	rdot2		=	(D/S) * ( E*rdot[0] - L*rdot[3] - S*thdot2 );

	if( rdot2 < 0 ){
		if(fabs(rdot2)<RDOT_ZERO and ir_change>NCHANGE){ rsign *= -1; ir_change=0; }
		throw(1);
	}
	if( thdot2 < 0 ){
		if(fabs(thdot2)<RDOT_ZERO and ith_change>NCHANGE){thsign *= -1; ith_change=0;}
		throw(2);
	}

	rdot[1]		=	rsign  * sqrt( rdot2  );
	rdot[2]		=	thsign * sqrt( thdot2 );
}
/* =================================== */


/******* Propagate a photon from a position pos *********/
void photon::propagate( double src[] ) {

	int			iter,i;
	double		diff,err,tau,tau_c,dtau = 0.01,halfpi=M_PI/2.;
	double		*rvec_c,*rvec_tmp,*rvec_err;
	bool		stopping;

	gsl_odeiv2_system				sys = { func , NULL , 4, (void*)(this) };
	const gsl_odeiv2_step_type 		*step = gsl_odeiv2_step_rk8pd;
	gsl_odeiv2_driver 				*driv =
		gsl_odeiv2_driver_alloc_standard_new (&sys, step, dtau, INTG_ERROR, INTG_ERROR, 1.0, 1.0);
	gsl_odeiv2_step					*stepper = gsl_odeiv2_step_alloc ( step , 4 );


	/* INITIALIZE SOME VARIALES */

	rvec_c		=	new double[4];
	rvec_tmp	=	new double[4];
	rvec_err	=	new double[4];
	for( i=0 ; i<4 ; i++ ){ rvec_c[i] = src[i]; rvec[i] = src[i]; rvec_tmp[i] = src[i];}
	tau			=	0;
	tau_c		=	0;
	ir_change	=	0;
	ith_change	=	0;
	stopping	=	false;

	/* ------------------------ */



	/* START THE INTEGRATION LOOP */
	for ( iter=1 ; iter<= NLOOP ; iter++ ){

		try{
			gsl_odeiv2_step_apply ( stepper , tau , dtau, rvec_tmp, rvec_err, NULL, NULL , &sys);
		}catch( int ie ){
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
		}

		/* If at the Horizon, stop */
		if ( (rvec_tmp[1]-rh) < STP_ZERO ) {stopping	=	true;}


		/* Are we at the disk at halfpi ? */
		if( fabs( rvec_tmp[2] - halfpi ) < STP_ZERO ) { stopping = true;}

		// Or are we close to HALF_PI? //
		if( ( rvec_tmp[2] - halfpi ) >= STP_ZERO and not stopping ){
			tau		=	tau_c;
			dtau	/=	(DTAU_FAC*DTAU_FAC);
			for( i=0 ; i<4 ; i++){ rvec_tmp[i] = rvec_c[i]; }
			continue;
		}


		/* If we get here, things are ok. Update saved rvec_c and tau_c */
		for( i=0 ;i<4 ; i++ ) { rvec_c[i] = rvec_tmp[i]; rvec[i] = rvec_tmp[i];}
		tau_c		=	tau;
		ir_change++; ith_change++;


		_print_xyz();
		//_const_of_motion();
		//for( i=0 ;i<4 ; i++ ) { printf("%3.3e ",rvec[i]);} printf("%3.3e %3.3e %3d\n",tau,dtau,iter);

		if ( stopping ) {
			try{
				calc_rdot();
				Tau		=	tau_c;
			}catch( int ie ){
				diff = sqrt(-1);
				rdot[0] = diff;rdot[1] = diff;rdot[2] = diff; rdot[3] = diff;Tau = 0;
			}
			break;
		}
	}
	/* -------------------------- */


	/* CLEAR ALLOCATED MEMORY */
	gsl_odeiv2_driver_free( driv );
	delete[] rvec_c;
	delete[] rvec_tmp;
	delete[] rvec_err;
	/* ----------------------- */


}
/* ==================================================== */


/******** check constants of motion for debugging **********/
void photon::_const_of_motion(){
	double	e,l,q,int1;
	calc_rdot();
	e		=	(1-2*r/S)*rdot[0] + (2*a*r*sin2/S)*rdot[3];
	l		=	(-2*a*r*sin2/S)*rdot[0] + (r2+a2+2*a2*r*sin2/S)*sin2*rdot[3];
	q		=	(S*S*rdot[2]*rdot[2]) + cos2*( a2*(-e*e) + l*l/sin2);
	int1	=	-(1-2*r/S)*rdot[0]*rdot[0] -(4*a*r*sin2/S)*rdot[0]*rdot[3] + \
			(S/D)*rdot[1]*rdot[1] + S*rdot[2]*rdot[2] + (r2+a2+2*a2*r*sin2/S)*sin2*rdot[3]*rdot[3];


    printf("%3.3e %3.3e %3.3e %3.3e\n",e,l,q,int1);
}
// ========================================================= /


/******* INTEGRATION FUNCTION ********/
int func (double t, const double r[], double drdt[], void *params) {
	int		i;
	photon	*ph = &(*(photon *)params);
	for( i=0 ; i<4 ; i++ ){ ph->rvec[i] = r[i];}
	ph->calc_rdot();
	for( i=0 ; i<4 ; i++ ){ drdt[i] = ph->rdot[i];}

	return GSL_SUCCESS;
}
/* ================================= */


/******** Print coordinates in x,t,z for debugging **********/
void photon::_print_xyz(){
	double	x,y,z,r,th,phi;
	r = rvec[1]; th=rvec[2]; phi=rvec[3];
    x		=	r*cos(phi)*sin(th);
    y		=	r*sin(phi)*sin(th);
    z		=	r*cos(th);
    printf("%3.3e %3.3e %3.3e\n",x,y,z);
}
// ========================================================= /


} /* namespace gr */
