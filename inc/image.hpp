/*
 * image.hpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#ifndef INC_IMAGE_HPP_
#define INC_IMAGE_HPP_

#include "flash.hpp"

namespace gr {

class image {
	// -------- Private Variables -------- //
	double		sin0,cos02,sin02;;


	// -------- Private Function -------- //

public:
	// -------- Public Variables -------- //
	double		a,theta_o,imsize,rms;
	double		*rvec;


	// -------- Public Functions -------- //
	image( double , double , double );
	virtual ~image();
	void project_image( int );
	void proj_xy( double x, double y, double* out );
};

} /* namespace gr */

#endif /* INC_IMAGE_HPP_ */
