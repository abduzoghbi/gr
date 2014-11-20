/*
 * flash.hpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#ifndef INC_FLASH_HPP_
#define INC_FLASH_HPP_

#include "tetrad.hpp"
#include "photon.hpp"


namespace gr {

class flash {

	// -------- Private Variables -------- //
	double		a;

	// -------- Private Functions -------- //
	void 		const_of_motion( double , double );

public:

	// -------- Public Variables  -------- //
	tetrad*		tet;
	double		*src,*consts,*rvec,*rdot,tau;
	int			*sign;


	// -------- Public Functions -------- //
	flash( double[] , double[] , double );
	virtual ~flash();
	void		illum(double ,double );
	void		illum( int );
};

} /* namespace gr */

#endif /* INC_FLASH_HPP_ */
