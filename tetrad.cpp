/*
 * tetrad.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#include "inc/tetrad.hpp"

namespace gr {

tetrad::tetrad( double rvec_[] , double drdt[] , double spin) {
	a		=	spin; a2 = a*a;
	rvec	=	new double[4];
	tetrad_vec	=	new double*[4];
	for( int i=0 ; i<4 ; i++ ){ rvec[i] = rvec_[i]; tetrad_vec[i] = new double[4];}
	mu		=	1.0;
	_metric_elements();
	_drdt_limits( drdt );
	td		=	sqrt(-1/(gtt + 2*gtp*drdt[2] + gpp*drdt[2]*drdt[2] +
				grr*drdt[0]*drdt[0] + gthth*drdt[1]*drdt[1]));
	rd		=	td * drdt[0];
	thd		=	td * drdt[1];
	pd		=	td * drdt[2];


	_calc_consts_of_motion();
	calc_tetrad();
}

tetrad::~tetrad() {
	delete[] rvec;
	for(int i=0;i<4;i++ ){delete[] tetrad_vec[i];} delete[] tetrad_vec;
}

/******* Useful Variables *********/
void tetrad::_setup_vars(){
	r		=	rvec[1]; r2 = r*r;
	sin2	=	sin( rvec[2] ); sin2 = sin2*sin2;
	cos2	=	1. - sin2;
	D		=	r2 - 2*r + a2;
	S		=	r2 + a2*cos2;
	DS		=	D*S;
	A		=	(r2+a2)*(r2+a2) - a2*D*sin2;
}
// =============================== //


/******* Metric Elements *********/
void tetrad::_metric_elements(){
	_setup_vars();
	gtt			=	-(1-2*r/S);
	gtp			=	-2*a*r*sin2/S;
	gpp			=	(r2+a2+2*a2*r*sin2/S)*sin2;
	grr			=	S/D;
	gthth		=	S;
}
// =============================== /


/******* Limits on drdt given rvec *********/
void tetrad::_drdt_limits( double drdt[] ){
	/* Assumes metric elements have been calculated */
	/* aa = gtt + 2*gtp*drdt[2] + gpp*drdt[2]^2 + grr*drdt[0]^2 + gthth*drdt[1]^2 has to be negative */

	// requirement for dpdt: gpp*dpdt^2 + 2gtp*dpdt +gtt>0 //
	double	del,cond1;
	double	Phi1,Phi2;
	del		=	gtp*gtp - gpp*gtt;
	Phi1	=	(-gtp - sqrt(del))/gpp;
	Phi2	=	(-gtp + sqrt(del))/gpp;
	cond1	=	gpp*drdt[2]*drdt[2] + 2*gtp*drdt[2] + gtt;
	if( drdt[2]<Phi1 or drdt[2]>Phi2){
		printf("dp_dt value not allowed, it should be between [%g,%g]. Value given: %g\n",
				Phi1,Phi2,drdt[2]);
		exit(1);
	}

	/* Now given dpdt, what are the requirements on drdt, dthdt
	 * gpp*dpdt^2 + 2gtp*dpdt +gtt +grr*drdt^2 + gthth*dthdt^2>0 */
	double	del2,dum1,dum2,cond2,cond3;
	dum1		=	gpp*grr*drdt[0]*drdt[0];
	dum2		=	gpp*gthth*drdt[1]*drdt[1];
	del2		=	del - dum1 - dum2;
	if( del2<0 ){ // we don't have solutions to the equation, so value always positive //
		cond2	=	(del - dum2)/(gpp*grr);
		if(cond2>0){
			cond2	=	sqrt(cond2);
			printf("*** dr_dt not allowed, for the given dth_dt, it should be between [%g,%g]. Value given: %g ***\n",
						-cond2,cond2,drdt[0]);
		}else{
			printf("*** No dr_dt is available for the given dth_dt, please change it. ***\n");
		}

		cond3	=	(del - dum1)/(gpp*gthth);
		if(cond3>0){
			cond3 = sqrt(cond3);
			printf("*** dth_dt not allowed, for the given dr_dt, it should be between [%g,%g]. Value given: %g ***\n",
						-cond3,cond3,drdt[1]);
		}else{
			printf("*** No dth_dt is available for the given dr_dt, please change it.***\n");
		}
		exit(1);
	}
	/* There is one more requirement for dpdt that makes it fall outside the range of two solutions
	 * so that gpp*dpdt^2 + 2gtp*dpdt +gtt +grr*drdt^2 + gthth*dthdt^2 > 0
	 * and this is satisfied most of the time
	 */

}
// ========================================= /

/******** Calculate const_of_motion given rdot **********/
void tetrad::_calc_consts_of_motion(){
	// _metric_elements should have been called in the initializer //
	E		=	-gtt * td - gtp*pd;
	L		=	gtp * td + gpp * pd;
	Q		=	S*S*thd*thd + cos2*( a2*(mu*mu - E*E) + L*L/sin2 );
}
// ===================================================== /

/******** Calculate tetrad elements **********/
void tetrad::calc_tetrad(){

	/* e_t */
	tetrad_vec[0][0]	=	td;
	tetrad_vec[0][1]	=	rd;
	tetrad_vec[0][2]	=	thd;
	tetrad_vec[0][3]	=	pd;

	/*
	double dum;
	dum = gtt*td*td + 2*gtp*td*pd + gpp*pd*pd + grr*rd*rd + gthth*thd*thd;
	printf("%g\n",dum);exit(0);
	*/


	/* e_phi */
	double		x,y,xy = sqrt( gtt*L*L + 2*gtp*L*E + gpp*E*E );
	x	=	L/xy;
	y	=	E/xy;
	tetrad_vec[3][0]	=	x;
	tetrad_vec[3][1]	=	0;
	tetrad_vec[3][2]	=	0;
	tetrad_vec[3][3]	=	y;

	/*
	double dum1,dum2;
	dum1 = gtt*td*x + gtp*td*y + gtp*pd*x + gpp*y*pd ;
	dum2 = gtt*x*x + 2*gtp*x*y + gpp*y*y ;
	printf("%g %g\n",dum1,dum2);exit(0);
	*/



	/* e_r */
	double	A,B,C,denom;
	denom	=	sqrt(
				E*E * gtp*gtp * x*x + 2*E*gtp*gtt*L*x*x + gtt*gtt*L*L*x*x -
				grr*gtp*gtp*gtt*rd*rd*x*x + gpp*grr*gtt*gtt*rd*rd*x*x +
				2*E*E*gpp*gtp*x*y + 2*E*gtp*gtp*L*x*y + 2*E*gpp*gtt*L*x*y +
				2*gtp*gtt*L*L*x*y - 2*grr*gtp*gtp*gtp*rd*rd*x*y +
				2*gpp*grr*gtp*gtt*rd*rd*x*y + E*E*gpp*gpp*y*y + 2*E*gpp*gtp*L*y*y +
				gtp*gtp*L*L*y*y - gpp*grr*gtp*gtp*rd*rd*y*y + gpp*gpp*grr*gtt*rd*rd*y*y
				);
	A	=	(1/(gtt*x + gtp*y)) * ((	(sqrt(grr) * gtp * gtt * rd * x * x)  +
										(sqrt(grr) * gtp * gtp * rd * x * y)  +
										(sqrt(grr) * gtt * gpp * rd * x * y)  +
										(sqrt(grr) * gtp * gpp * rd * y * y)
									)/denom);
	B	=	(1/(grr*gtt*x + grr*gtp*y)) * ((
				(sqrt(grr) * E * gtp * gtt * x*x)   + (sqrt(grr) * L * gtt * gtt * x*x) +
				(sqrt(grr) * E * gtp * gtp * x*y)   + (sqrt(grr) * E * gtt * gpp * x*y) +
				(sqrt(grr) * L * gtp * gtt * x*y*2) + (sqrt(grr) * E * gtp * gpp * y*y) +
				(sqrt(grr) * L * gtp * gtp * y*y)
									)/denom);
	C	=	-(sqrt(grr) * rd * (gtt*x + gtp*y)) / denom;


	tetrad_vec[1][0]	=	A;
	tetrad_vec[1][1]	=	B;
	tetrad_vec[1][2]	=	0;
	tetrad_vec[1][3]	=	C;

	/*
	double dum1,dum2,dum3;
	dum1 = gtt*A*A + 2*gtp*A*C + gpp*C*C + grr*B*B ;
	dum2 = gtt*A*td + gtp*A*pd + gtp*C*td + gpp*C*pd + grr*B*rd;
	dum3 = gtt*A*x + gtp*A*y + gtp*C*x + gpp*C*y;
	printf("%g %g %g\n",dum1,dum2,dum3);exit(0);
	*/


	/* e_th */
	double		alpha,beta,gamma,delta;
	denom		=	sqrt(B*B*E*E*grr*gtp*gtp*x*x + 2*B*B*E*grr*gtp*gtt*L*x*x +
		    		B*B*grr*gtt*gtt*L*L*x*x + 2*B*C*E*grr*gtp*gtp*gtp*rd*x*x -
					2*B*C*E*gpp*grr*gtp*gtt*rd*x*x + 2*B*C*grr*gtp*gtp*gtt*L*rd*x*x -
					2*B*C*gpp*grr*gtt*gtt*L*rd*x*x + C*C*grr*gtp*gtp*gtp*gtp*rd*rd*x*x -
					2*C*C*gpp*grr*gtp*gtp*gtt*rd*rd*x*x + C*C*gpp*gpp*grr*gtt*gtt*rd*rd*x*x +
					C*C*gthth*gtp*gtp*gtp*gtp*thd*thd*x*x - 2*C*C*gpp*gthth*gtp*gtp*gtt*thd*thd*x*x -
					B*B*grr*gthth*gtp*gtp*gtt*thd*thd*x*x +
					C*C*gpp*gpp*gthth*gtt*gtt*thd*thd*x*x +
					B*B*gpp*grr*gthth*gtt*gtt*thd*thd*x*x + 2*B*B*E*E*gpp*grr*gtp*x*y +
					2*B*B*E*grr*gtp*gtp*L*x*y + 2*B*B*E*gpp*grr*gtt*L*x*y +
					2*B*B*grr*gtp*gtt*L*L*x*y + 2*B*C*E*gpp*grr*gtp*gtp*rd*x*y -
					2*A*B*E*grr*gtp*gtp*gtp*rd*x*y - 2*B*C*E*gpp*gpp*grr*gtt*rd*x*y +
					2*A*B*E*gpp*grr*gtp*gtt*rd*x*y + 2*B*C*grr*gtp*gtp*gtp*L*rd*x*y -
					2*B*C*gpp*grr*gtp*gtt*L*rd*x*y - 2*A*B*grr*gtp*gtp*gtt*L*rd*x*y +
					2*A*B*gpp*grr*gtt*gtt*L*rd*x*y - 2*A*C*grr*gtp*gtp*gtp*gtp*rd*rd*x*y +
					4*A*C*gpp*grr*gtp*gtp*gtt*rd*rd*x*y -
					2*A*C*gpp*gpp*grr*gtt*gtt*rd*rd*x*y -
					2*B*B*grr*gthth*gtp*gtp*gtp*thd*thd*x*y - 2*A*C*gthth*gtp*gtp*gtp*gtp*thd*thd*x*y +
					2*B*B*gpp*grr*gthth*gtp*gtt*thd*thd*x*y +
					4*A*C*gpp*gthth*gtp*gtp*gtt*thd*thd*x*y -
					2*A*C*gpp*gpp*gthth*gtt*gtt*thd*thd*x*y + B*B*E*E*gpp*gpp*grr*y*y +
					2*B*B*E*gpp*grr*gtp*L*y*y + B*B*grr*gtp*gtp*L*L*y*y -
					2*A*B*E*gpp*grr*gtp*gtp*rd*y*y + 2*A*B*E*gpp*gpp*grr*gtt*rd*y*y -
					2*A*B*grr*gtp*gtp*gtp*L*rd*y*y + 2*A*B*gpp*grr*gtp*gtt*L*rd*y*y +
					A*A*grr*gtp*gtp*gtp*gtp*rd*rd*y*y - 2*A*A*gpp*grr*gtp*gtp*gtt*rd*rd*y*y +
					A*A*gpp*gpp*grr*gtt*gtt*rd*rd*y*y - B*B*gpp*grr*gthth*gtp*gtp*thd*thd*y*y +
					A*A*gthth*gtp*gtp*gtp*gtp*thd*thd*y*y + B*B*gpp*gpp*grr*gthth*gtt*thd*thd*y*y -
					2*A*A*gpp*gthth*gtp*gtp*gtt*thd*thd*y*y + A*A*gpp*gpp*gthth*gtt*gtt*thd*thd*y*y);

	alpha	=	(1/(gtt*x + gtp*y)) * ((	(B*sqrt(grr)*sqrt(gthth)*gtp*gtt*thd*x*x) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gtp*thd*x*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtt*gpp*thd*x*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gpp*thd*y*y)

			) / denom);
	beta	=	(1/(grr*gtt*x + grr*gtp*y)) * ((
										   -(C*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtt*thd*x*x) +
											(C*sqrt(grr)*sqrt(gthth)*gpp*gtt*gtt*thd*x*x) -
											(C*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtp*thd*x*y) +
											(C*sqrt(grr)*sqrt(gthth)*gtp*gtt*gpp*thd*x*y) +
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtt*thd*x*y) -
											(A*sqrt(grr)*sqrt(gthth)*gtt*gtt*gpp*thd*x*y) +
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtp*thd*y*y) -
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtt*gpp*thd*y*y)
			) / denom );
	gamma	=	(1/(gthth*gtt*x + gthth*gtp*y)) * (
											(B*sqrt(grr)*sqrt(gthth)*gtp*gtt*E*x*x) +
											(B*sqrt(grr)*sqrt(gthth)*gtt*gtt*L*x*x) +
											(C*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtt*rd*x*x) -
											(C*sqrt(grr)*sqrt(gthth)*gtt*gtt*gpp*rd*x*x) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gtp*E*x*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtt*gpp*E*x*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gtt*L*x*y*2) +
											(C*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtp*rd*x*y) -
											(C*sqrt(grr)*sqrt(gthth)*gtp*gtt*gpp*rd*x*y) -
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtt*rd*x*y) +
											(A*sqrt(grr)*sqrt(gthth)*gtt*gtt*gpp*rd*x*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gpp*E*y*y) +
											(B*sqrt(grr)*sqrt(gthth)*gtp*gtp*L*y*y) -
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtp*gtp*rd*y*y) +
											(A*sqrt(grr)*sqrt(gthth)*gtp*gtt*gpp*rd*y*y)

			) / denom;
	delta	=	-(B*sqrt(grr)*sqrt(gthth)*thd*(gtt*x + gtp*y)) / denom;

	tetrad_vec[2][0]	=	alpha;
	tetrad_vec[2][1]	=	beta;
	tetrad_vec[2][2]	=	gamma;
	tetrad_vec[2][3]	=	delta;

	/*
	double dum1,dum2,dum3,dum4;
	dum1 = gtt*alpha*alpha + 2*gtp*alpha*delta + gpp*delta*delta + grr*beta*beta + gthth*gamma*gamma ;
	dum2 = gtt*alpha*td + gtp*alpha*pd + gtp*delta*td + gpp*delta*pd + grr*beta*rd + gthth*gamma*thd;
	dum3 = gtt*alpha*x + gtp*alpha*y + gtp*delta*x + gpp*delta*y;
	dum4 = gtt*alpha*A + gtp*alpha*C + gtp*delta*A + gpp*delta*C + gtt*beta*B;
	printf("%g %g %g %g\n",dum1,dum2,dum3,dum4);exit(0);
	*/

}
// ========================================== /


} /* namespace gr */
