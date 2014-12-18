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

const static H5std_string	DATASET_NAME( "DATA" ),ATTR("src4_drdt4_a");
const static int ATTR_SIZE = 9;

class flash {

	// -------- Private Variables -------- //
	double		a;

	// -------- Private Functions -------- //
	void 		const_of_motion( double , double , double* , int* );
	void		write_hdf5( double* ,int , int , const string );

public:

	// -------- Public Variables  -------- //
	tetrad*		tet;
	double		*src,*drdt;


	// -------- Public Functions -------- //
	flash( double[] , double[] , double );
	virtual ~flash();
	void		illum(double ,double , double* );
	void		illum( int , const string fname="illum.h5" );
	static void	read_hdf5( const string fname , double *&data , double* &attr , int *dim );
};


class disk {


	// -------- Private Variables  -------- //
	double		*data,*attr,rLim[2],a,rms,rh,*rbinL,*rbinc,*area,*redshift,*emiss;
	int			nph,ncol,nr,*rcount;
	bool		islog;

	// -------- Private Functions  -------- //
	void		radial_bins();
	void 		mtm_at_position(double E,double L, double Q, double* pos , double* rdot );
	static double disk_area_int_func(double r, void * params){
		double		res,a2 = *(double*) params,r2=r*r;
		// srqt( grr * gpp )
		res		=	sqrt( (r2/(r2-2*r+a2)) * (r2+a2+2*a2/r) );
		return res;
	}

public:
	// -------- Public Variables  -------- //


	// -------- Public Functions -------- //
	disk( const string fname , int nr , double* rlim , bool islog=true );
	virtual ~disk();
	void		emissivity();
	static void	disk_velocity( double r , double* vel , double a, double rms );
	static double p_dot_v( double r, double th , double* p , double* v , double a );
};

} /* namespace gr */

#endif /* INC_FLASH_HPP_ */
