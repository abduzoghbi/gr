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
#include <gsl/gsl_histogram2d.h>




namespace gr {

const static H5std_string	FLASH_DATASET( "DATA" ),FLASH_ATTR_LABEL("src4_drdt4_a");
const static int FLASH_ATTR_SIZE = 9;

class flash {

	// -------- Private Variables -------- //
	double		a;

	// -------- Private Functions -------- //
	void 		const_of_motion( double , double , double* , int* , double *src , tetrad& );
	void		write_hdf5( double* ,int , int , const string );

public:

	// -------- Public Variables  -------- //
	tetrad*		tet;
	double		*src;


	// -------- Public Functions -------- //
	flash( double[] , double[] , double );
	virtual ~flash();
	void		illum(double ,double , double* , double *src, tetrad&);
	void		illum( int , const string fname="illum.h5" );
	static void	read_hdf5( const string fname , double *&data , double* &attr , int *dim );
};

} /* namespace gr */

#endif /* INC_FLASH_HPP_ */
