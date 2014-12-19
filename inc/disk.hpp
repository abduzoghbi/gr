/*
 * disk.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: abzoghbi
 */

#ifndef SRC_INC_DISK_HPP_
#define SRC_INC_DISK_HPP_

#include "flash.hpp"
#include <gsl/gsl_histogram2d.h>

namespace gr {

class disk {

	// -------- Private Variables  -------- //
	double		*data,*attr,a,rms,rh,*rL,*phiL;
	int			nph,ncol,nr,nphi;
	bool		rlog;

	gsl_histogram2d		*ph_count,*redshift;
	double				**area;


	// -------- Private Functions  -------- //

public:

	// -------- Public Variables  -------- //


	// -------- Public Functions  -------- //
	disk( const string fname , int nr , double* rlim , int np=1, bool rlog=true );
	virtual ~disk();
	static double proper_disk_area( double r , double a , double rms );
	static void	disk_velocity( double r , double* vel , double a, double rms );
	static double p_dot_v( double r, double th , double* p , double* v , double a );
	static void image2disk( int npix , double imsize,double *rvec , double a , double rms , double **ret );
	static void proj_xy( double x , double y , double *rvec , double a , double rms,double*out );
	void emissivity();
	void tf();
};

} /* namespace gr */

#endif /* SRC_INC_DISK_HPP_ */
