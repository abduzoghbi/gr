/*
 * photon.hpp
 *
 *  Created on: Nov 18, 2014
 *      Author: abzoghbi
 */

#ifndef INC_PHOTON_HPP_
#define INC_PHOTON_HPP_

#include <cmath>
#include <iostream>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>


using namespace std;

namespace gr {

const static double INTG_ERROR 	=	1E-8;
const static double	DTAU_FAC	=	2;
const static double	RDOT_ZERO	=	1e-4;
const static double RINF		=	20;
const static double STP_ZERO	=	1e-4;
const static int	NLOOP		=	400;
const static int	NCHANGE		=	5;

class photon {

	/**** Private Variables ****/
	double		a2,rh;
	double		r,r2,sin1,sin2,cos1,cos2,D,S,DS,A;
	int			ir_change,ith_change;


	/**** Private Functions ****/
	void		_setup_vars();
	void		_const_of_motion();
	void		_print_xyz();

public:
	/**** Public Variables ****/
	double		E,L,Q,a;
	int			rsign,thsign;
	double		*rvec,*rdot;


	/**** Public Functions ****/
	photon(double,double,double,double,int,int);
	virtual ~photon();
	void 		calc_rdot();
	void		propagate( double* src );
};


// INTEGRATION FUNCTION //
int func (double t, const double r[], double rdot[], void *params);

} /* namespace gr */

#endif /* INC_PHOTON_HPP_ */
