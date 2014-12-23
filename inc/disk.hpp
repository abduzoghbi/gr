/*
 * disk.hpp
 *
 *  Created on: Dec 19, 2014
 *      Author: abzoghbi
 */

#ifndef SRC_INC_DISK_HPP_
#define SRC_INC_DISK_HPP_

#include "flash.hpp"
#include "image.hpp"
#include <gsl/gsl_histogram2d.h>

namespace gr {

class disk {

	// -------- Private Variables  -------- //
	double		*data,*attr,a,rms,rh,*rL,*phiL;
	int			nph,ncol,nr,nphi;
	bool		rlog;

	gsl_histogram2d		*r_phi_count,*r_phi_redshift,*r_phi_time;
	double				**area;


	// -------- Private Functions  -------- //
	void flash_to_r_phi( double& tsource );

public:

	// -------- Public Variables  -------- //


	// -------- Public Functions  -------- //
	disk( const string flashfile , int nr , double* rlim , int np=1, bool rlog=true );
	virtual ~disk();
	static double proper_disk_area( double r , double a , double rms );
	static void	disk_velocity( double r , double* vel , double a, double rms );
	static double p_dot_v( double r, double th , double* p , double* v , double a );
	void emissivity();
	void tf( int ntime , double* tLim, int nenergy , double* enLim );
};

} /* namespace gr */

#endif /* SRC_INC_DISK_HPP_ */
