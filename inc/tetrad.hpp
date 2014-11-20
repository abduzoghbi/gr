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

namespace gr {

class tetrad {
	// -------- Private Variables -------- //
	double		a2,r,r2,D,S,DS,sin2,cos2,A;
	double		gtt,gtp,gpp,grr,gthth;
	double		td,rd,thd,pd; /* rdot elements */


	// -------- Private Functions -------- //
	void _setup_vars();
	void _drdt_limits( double[] );
	void _metric_elements();
	void _calc_consts_of_motion();

public:
	// -------- Public Variables -------- //
	double		E,L,Q,a,mu;
	double		*rvec;
	double		**tetrad_vec;


	// -------- Public Functions -------- //
	/* Initializer takes, position vector, drdt 3-vector, spin */
	tetrad( double[] , double[] , double );
	void calc_tetrad();
	virtual ~tetrad();
	double operator()(int i, int j) {return tetrad_vec[i][j];}
};

} /* namespace gr */

#endif /* INC_TETRAD_HPP_ */
