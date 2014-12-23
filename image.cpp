/*
 * image.cpp
 *
 *  Created on: Dec 22, 2014
 *      Author: abzoghbi
 */

#include "inc/image.hpp"

namespace gr {

image::image( double spin, double th0 , double imsize_ , int npix_  ) {
	a		=	spin;
	imsize	=	imsize_;
	npix	=	npix_;
	if( npix%2 == 0) npix++; // odd to include zero

	rvec	=	new double[4];
	rvec[0]	=	0; rvec[1] = 1000; rvec[2] = th0; rvec[3] = 0;

	// ISCO //
	double		z1,z2;
	z1		=	(1 + pow(1-a*a,1./3)*( pow(1+a,1/3) + pow(1-a,1/3) ));
	z2		=	sqrt( 3*a*a + z1*z1 );
	rms		=	3 + z2	- sqrt( (3-z1)*(3+z1+2*z2) );

	// data container: npix * npix * 4, flattened //
	data	=	new double[ npix*npix*4 ];

	image2disk();
}

image::image( const string fname ){
	double		*attr;
	int			dim[1];


	read_hdf5( fname , data , attr , dim );
	imsize		=	attr[0];
	npix		=	attr[1];
	a			=	attr[2];
	rvec	=	new double[4];
	rvec[0]	=	0; rvec[1] = 1000; rvec[2] = attr[3]; rvec[3] = 0;


	// ISCO //
	double		z1,z2;
	z1		=	(1 + pow(1-a*a,1./3)*( pow(1+a,1/3) + pow(1-a,1/3) ));
	z2		=	sqrt( 3*a*a + z1*z1 );
	rms		=	3 + z2	- sqrt( (3-z1)*(3+z1+2*z2) );
	delete[] attr;
}

image::~image() {
	delete[] rvec; delete[] data;
}

/******** Project image frame to disk  **********/
void image::image2disk(){

#pragma omp parallel
{
	int			nside;
	double		unit,x,y,rvec_tmp[4],A,RMS;
	A		=	a;
	RMS		=	rms;

	for( int i=0 ;i<4;i++) rvec_tmp[i] = rvec[i];

	nside	=	(npix-1)/2;
	unit	=	imsize/npix;
	#pragma omp for schedule(guided)
	for( int i=-nside ; i<=nside ; i++ ){
		x	=	unit * i;
		for( int j=-nside; j<=nside ; j++ ){
			y	=	unit * j;
			proj_xy( x , y , &(data[(i+nside)*(npix*4)+(j+nside)*4]) , rvec_tmp , A, RMS );
		}
	}
}// end pragma parallel
}
/*=================================================*/

/******** Project a photon from x,y on image r,phi,g,t on disk **********/
void image::proj_xy( double x , double y , double* out , double *rvec , double a, double rms ){
	double			E,L,Q,vel[4],g,sino,sino2,coso2,rdot[4],rvec_tmp[4];
	int				thsign;

	if(fabs(x)<1e-6) { x = (x<0)?-1e-6:1e-6;}
	if(fabs(y)<1e-6) { y = (y<0)?-1e-6:1e-6;}

	for( int i=0 ; i<4 ; i++ ) rvec_tmp[i] = rvec[i];

	sino	=	sin( rvec[2] ); sino2 = sino*sino; coso2 = 1-sino2;


	E		=	1.0;
	L		=	x*sino;
	Q		=	y*y + coso2*( -a*a*E*E + L*L/sino2);
	thsign	=	(y<0)?1:-1;
	photon		ph( E , L , Q , a , -1 , thsign );
	ph.propagate( rvec_tmp , rdot );
	if( ph.rh_stop ){
		g	=	0;
	}else{
		disk::disk_velocity( rvec_tmp[1] , vel , a , rms );
		g		=	-disk::p_dot_v( rvec_tmp[1] , rvec_tmp[2] , rdot , vel , a );
	}
	out[0]	=	rvec_tmp[0];		// t
	out[1]	=	rvec_tmp[1];		// r
	out[2]	=	rvec_tmp[3];		// phi
	out[3]	=	g;					// g
}
/*=========================================================================*/


/******** Write the output of illum to an hdf5 file **********/
void image::write_hdf5( const string fname ){
	H5std_string  	FILE_NAME( fname );
	hsize_t     	dim[1],dims[]={IMAGE_ATTR_SIZE};;
	dim[0]		=	npix*npix*4;

	try{

		H5::Exception::dontPrint();

		H5::H5File 		file( FILE_NAME, H5F_ACC_TRUNC );
		H5::DataSpace 	dataspace( 1, dim );


		H5::DataSet 	dataset = file.createDataSet( IMAGE_DATASET, H5::PredType::NATIVE_DOUBLE, dataspace );
		dataset.write( data, H5::PredType::NATIVE_DOUBLE );

		double attr[IMAGE_ATTR_SIZE] = {imsize,1.*npix,a,rvec[2]};
		dataset.createAttribute( IMAGE_ATTR_LABEL , H5::PredType::NATIVE_DOUBLE ,
				H5::DataSpace ( 1, dims )).write(H5::PredType::NATIVE_DOUBLE, attr);


	} catch ( H5::Exception &err ) {
		err.printError();
	}
}
/* ========================================================= */

/******** Read a saved hdf5 file **********/
void image::read_hdf5( const string fname , double *&data , double* &attr , int* dim ){

	H5std_string  	FILE_NAME( fname );
	hsize_t			dims[1];


	H5::H5File		file( FILE_NAME, H5F_ACC_RDONLY );
	H5::DataSet		dataset 	= 	file.openDataSet( IMAGE_DATASET );
	H5::DataSpace 	filespace 	= 	dataset.getSpace();
	H5::Attribute	attribute	=	dataset.openAttribute(IMAGE_ATTR_LABEL);
	attr	=	new double[IMAGE_ATTR_SIZE];
	attribute.read( H5::PredType::NATIVE_DOUBLE , attr );

	filespace.getSimpleExtentDims( dims );
	data	=	new double[dims[0]];
	dim[0]	=	dims[0];

	H5::DataSpace	dataspace( 1 , dims );
	dataset.read( data , H5::PredType::NATIVE_DOUBLE , dataspace , filespace );
}

/* ====================================== */

} /* namespace psdlag */
