/*
 * image.hpp
 *
 *  Created on: Dec 22, 2014
 *      Author: abzoghbi
 */

#ifndef INC_IMAGE_HPP_
#define INC_IMAGE_HPP_

#include "disk.hpp"

namespace gr {

const static H5std_string	IMAGE_DATASET( "DATA" ),IMAGE_ATTR_LABEL("imsize_npix_a_th0");
const static int IMAGE_ATTR_SIZE = 4;

class image {

	/******* Private Variables *********/
	double		*rvec,*data;

	/******* Private Functions *********/
	void image2disk();
	void proj_xy( double x, double y , double*out , double* rvec, double a , double rms );
public:
	/******* Public Variables *********/
	double		imsize,a,rms;
	int			npix;


	/******* Public Functions *********/
	image( double spin, double th0 , double imsize , int npix );
	image( const string fname );
	virtual ~image();
	double operator()(int i, int j, int k) {return data[i*(npix*4)+j*4+k];}
	void write_hdf5( const string fname );
	void read_hdf5( const string fname , double *&data , double* &attr , int* dim );
};

} /* namespace psdlag */

#endif /* INC_IMAGE_HPP_ */
