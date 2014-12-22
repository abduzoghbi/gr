/*
 * flash.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#include "inc/flash.hpp"

namespace gr {

flash::flash( double s[] , double drdt[] , double spin ) {
	a		=	spin;
	src		=	new double[4];
	for( int i=0 ; i<4 ; i++ ){ src[i]  = s[i]; }

	tet		=	new tetrad( src , drdt , spin );
}

flash::~flash() {
	delete[] src;
	delete tet;
}

/******** Const of motion for a given alpha,beta in source frame **********/
void flash::const_of_motion ( double alpha, double beta , double* consts , int* sign , double* src , tetrad& tet){
	double		pt[4],rd[4];
	int			i,j;

	/* Momentum in emitting source frame */
	pt[0]		=	1.0;
	pt[1]		=	sin(alpha)*sin(beta);
	pt[2]		=	sin(alpha)*cos(beta);
	pt[3]		=	cos(alpha);

	/* Momentum in coordinate frame */
	for( i=0 ; i<4 ; i++ ){
		rd[i]		=	0;
		for( j=0 ; j<4 ;j++ ){ rd[i] +=	tet(j,i) * pt[j];}
	}

	/* Now get constants of motion given rdot */
	double		r,S,sin2,cos2,E,L,Q;
	r		=	src[1];
	sin2	=	sin(src[2]);sin2 = sin2*sin2;
	cos2	=	1. - sin2;
	S		=	r*r + a*a*cos2;

	E		=	(1-2*r/S) * rd[0] + (2*a*r*sin2/S) * rd[3];
	L		=	(-2*a*r*sin2/S) * rd[0] + (r*r+a*a+2*a*a*r*sin2/S)*sin2 * rd[3];
	Q		=	S*S*rd[2]*rd[2] + cos2*(a*a*(-E*E) + L*L/sin2);

	consts[0]	=	E;
	consts[1]	=	L;
	consts[2]	=	Q;

    sign[0]		=	(rd[1]<0)?-1:1;
    sign[1]		=	(rd[2]<0)?-1:1;

}
/* ====================================================================== */


/******** illuminate the source with a photon at given alpha, beta **********/
void flash::illum ( double alpha , double beta , double* out , double* src, tetrad& tet){
	double		consts[3];
	int			sign[2];
	for(int i=0; i<4 ; i++ ) out[i] = src[i];
	const_of_motion( alpha , beta , consts , sign , src , tet );
	photon	ph( consts[0] , consts[1] , consts[2] , a , sign[0] , sign[1] );
	ph.propagate( &(out[0]) , &(out[4]) );
}
/* ======================================================================== */


/******** flash the source with isotropic num of photons **********/
void flash::illum( int num , const string fname ){
	int		i,j,k,nph = num*num, ncol=8;
	double	*data,*Alpha,*Beta;

	Alpha		=	new double[nph];
	Beta		=	new double[nph];
	data		=	new double[nph*ncol];

	k = 0;
	for( i=0 ; i<num ; i++ ){
		for( j=0 ; j<num ; j++ ){
			Alpha[k]	=	acos( 1- i*2.0/(num-1) );
			Beta[k]		=	2.*M_PI*j/(num);
			k++;
		}
	}
	/*
	double		SRC[4];
	for( int i=0 ; i<4 ; i++ ) SRC[i] = src[i];
	tetrad		TET(*tet);
	k = 1;
	illum( Alpha[k] , Beta[k] , &data[k*ncol] , SRC, TET );
	exit(0);
	*/
#pragma omp parallel
{
	// Some variables to avoid calling globals in parallel //
	double		SRC[4];
	for( int i=0 ; i<4 ; i++ ) SRC[i] = src[i];
	tetrad		TET(*tet);


	#pragma omp for schedule(dynamic,num)
	for( int k=0 ; k<nph ; k++ ){
		illum( Alpha[k] , Beta[k] , &data[k*ncol] , SRC, TET );
	}
} // end pragma omp parallel


	write_hdf5( data , nph , ncol , fname );

	delete[] Alpha;
	delete[] Beta;
	delete[] data;
}
/* ============================================================== */



/******** Write the output of illum to an hdf5 file **********/
void flash::write_hdf5( double *data , int nx, int ny , const string fname ){
	H5std_string  	FILE_NAME( fname );
	hsize_t     	dim[2],dims[]={ATTR_SIZE};
	dim[0] = nx; dim[1] = ny;
	try{

		H5::Exception::dontPrint();

		H5::H5File 		file( FILE_NAME, H5F_ACC_TRUNC );
		H5::DataSpace 	dataspace( 2, dim );


		H5::DataSet 	dataset = file.createDataSet( DATASET_NAME, H5::PredType::NATIVE_DOUBLE, dataspace );
		dataset.write( data, H5::PredType::NATIVE_DOUBLE );

		double attr[ATTR_SIZE] = {src[0],src[1],src[2],src[3],tet->td,tet->rd,tet->thd,tet->pd,a};
		dataset.createAttribute( ATTR , H5::PredType::NATIVE_DOUBLE ,
				H5::DataSpace ( 1, dims )).write(H5::PredType::NATIVE_DOUBLE, attr);


	} catch ( H5::Exception &err ) {
		err.printError();
	}
}
/* ========================================================= */


/******** Read a saved hdf5 file **********/
void flash::read_hdf5( const string fname , double *&data , double* &attr , int dim[] ){

	H5std_string  	FILE_NAME( fname );
	hsize_t			dims[2];


	H5::H5File		file( FILE_NAME, H5F_ACC_RDONLY );
	H5::DataSet		dataset 	= 	file.openDataSet( DATASET_NAME );
	H5::DataSpace 	filespace 	= 	dataset.getSpace();
	H5::Attribute	attribute	=	dataset.openAttribute(ATTR);
	attr	=	new double[ATTR_SIZE];
	attribute.read( H5::PredType::NATIVE_DOUBLE , attr );

	filespace.getSimpleExtentDims( dims );
	data	=	new double[dims[0]*dims[1]];
	dim[0]	=	dims[0]; dim[1] = dims[1];

	H5::DataSpace	dataspace( 2 , dims );
	dataset.read( data , H5::PredType::NATIVE_DOUBLE , dataspace , filespace );
}

/* ====================================== */

} /* namespace gr */
