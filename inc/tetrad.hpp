/*
 * tetrad.hpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#ifndef INC_TETRAD_HPP_
#define INC_TETRAD_HPP_

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <gsl/gsl_linalg.h>

namespace gr {

class tetrad {
	// -------- Private Variables -------- //
	double		a2,r,r2,D,S,DS,sin2,cos2,A;
	double		gtt,gtp,gpp,grr,gthth;


	// -------- Private Functions -------- //
	void _setup_vars();
	void _drdt_limits( double[] );
	void _metric_elements();
	void _calc_consts_of_motion();

public:
	// -------- Public Variables -------- //
	double		E,L,Q,a,mu;
	double		*rvec;
	double		**tetrad_vec,**itetrad_vec;
	double		td,rd,thd,pd; /* rdot elements */


	// -------- Public Functions -------- //
	/* Initializer takes, position vector, drdt 3-vector, spin */
	tetrad( double[] , double[] , double );
	tetrad(tetrad&);
	void calc_tetrad();
	void inv_tetrad();
	virtual ~tetrad();
	double operator()(int i, int j) {return tetrad_vec[i][j];}
	double operator()(int i, int j , int dum ) {return itetrad_vec[i][j];}
};

} /* namespace gr */

#endif /* INC_TETRAD_HPP_ */
