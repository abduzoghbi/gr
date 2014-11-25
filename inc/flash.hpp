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
#include <omp.h>
#include <H5Cpp.h>


namespace gr {

class flash {

	// -------- Private Variables -------- //
	double		a;

	// -------- Private Functions -------- //
	void 		const_of_motion( double , double , double* , int* );
	void		write_hdf5( double* ,int , int , const string );
	void		read_hdf5(const string , double* );

public:

	// -------- Public Variables  -------- //
	tetrad*		tet;
	double		*src;


	// -------- Public Functions -------- //
	flash( double[] , double[] , double );
	virtual ~flash();
	void		illum(double ,double , double* );
	void		illum( int , const string fname="illum.h5" );
};

} /* namespace gr */

#endif /* INC_FLASH_HPP_ */
