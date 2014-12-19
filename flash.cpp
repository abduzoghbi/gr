/*
 * flash.cpp
 *
 *  Created on: Nov 19, 2014
 *      Author: abzoghbi
 */

#include "inc/flash.hpp"

namespace gr {

flash::flash( double s[] , double drt[] , double spin ) {
	a		=	spin;
	src		=	new double[4];
	drdt	=	new double[3];
	for( int i=0 ; i<4 ; i++ ){ src[i]  = s[i]; }
	for( int i=0 ; i<3 ; i++ ){ drdt[i] = drt[i]; }

	tet		=	new tetrad( src , drdt , spin );
}

flash::~flash() {
	delete[] src;
	delete[] drdt;
	delete tet;
}

/******** Const of motion for a given alpha,beta in source frame **********/
void flash::const_of_motion ( double alpha, double beta , double* consts , int* sign ){
	double		pt[4],rd[4];
	int			i,j;
	double		dum;

	/* Momentum in emitting source frame */
	pt[0]		=	1.0;
	pt[1]		=	cos(alpha);
	pt[2]		=	sin(alpha)*cos(beta);
	pt[3]		=	sin(alpha)*sin(beta);

	/* Momentum in coordinate frame */
	for( i=0 ; i<4 ; i++ ){
		dum		=	0;
		for( j=0 ; j<4 ;j++ ){ dum +=	tet->tetrad_vec[j][i] * pt[j];}
		rd[i]	=	dum;
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
void flash::illum ( double alpha , double beta , double* out ){
	double		consts[3];
	int			sign[2];
	const_of_motion( alpha , beta , consts , sign );
	photon	ph( consts[0] , consts[1] , consts[2] , a , sign[0] , sign[1] );
	ph.propagate( src );
	out[0]		=	ph.Tau;
	for( int i=0 ; i<4 ; i++ ){ out[i+1] = ph.rvec[i]; out[i+5] = ph.rdot[i]; }
	//out[9] = consts[0]; out[10] = consts[1]; out[11] = consts[2];
}
/* ======================================================================== */


/******** flash the source with isotropic num of photons **********/
void flash::illum( int num , const string fname ){
	int		i,j,k,nph = num*num, ncol=9;
	double	*data,*Alpha,*Beta;

	Alpha		=	new double[nph];
	Beta		=	new double[nph];
	data		=	new double[nph*ncol];

	k = 0;
	for( i=0 ; i<num ; i++ ){
		for( j=0 ; j<num ; j++ ){
			Alpha[k]	=	acos( i*2.0/(num-1) - 1);
			Beta[k]		=	2*M_PI*j/(num-1);
			k++;
		}
	}

#pragma omp parallel
{
	#pragma omp for schedule(dynamic,num) private(k)
	for( k=0 ; k<nph ; k++ ){
		illum( Alpha[k] , Beta[k] , &data[k*ncol] );
	}
} // end pragma omp parallel


	write_hdf5( data , nph , ncol , fname );
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



/************************************************/
/**************** Disk CLASS ********************/
/************************************************/

disk::disk( const string fname , int nr_ , double* rlim , int np ,bool islog_ ) {

	// Read the flash file //
	int		dim[2];
	flash::read_hdf5( fname , data , attr , dim );
	nph		=	dim[0];
	ncol	=	dim[1];
	a		=	attr[8];


	// ISCO and Horizon radii //
	double		z1,z2;
	z1		=	(1 + pow(1-a*a,1./3)*( pow(1+a,1/3) + pow(1-a,1/3) ));
	z2		=	sqrt( 3*a*a + z1*z1 );
	rms		=	3 + z2	- sqrt( (3-z1)*(3+z1+2*z2) );
	rh		=	1 + sqrt(1-a*a);



	// Work out radial bins and their areas //
	nr		=	nr_;
	rLim[0]	=	rlim[0]; rLim[1] = rlim[1];
	if(rLim[0]<rms) rLim[0] = rms;
	islog	=	islog_;


	// Phi bins //
	if( np != 0 ) nphi = np;

	radial_bins();

}

disk::~disk(){
	delete[] data; delete[] attr; delete[] rbinL; delete[] rbinc;delete[] rcount;delete[] area;
	delete[] redshift; delete[] emiss;delete[] phiL;
}


/******* Radial bins and areas ********/
void disk::radial_bins(){

	rbinL		=	new double[nr+1];
	rbinc		=	new double[nr];
	rcount		=	new int[nr];
	area		=	new double[nr];
	redshift	=	new double[nr];
	emiss		=	new double[nr];
	phiL		=	new double[nphi+1];

	// Work out the radial bins //

	for( int ir=0 ; ir<nr ; ir++ ){
		if( islog ){
			rbinL[ir]	=	exp(log(rLim[0]) + ir*(log(rLim[1] / rLim[0])	)/nr);
		}else{
			rbinL[ir]	=	rLim[0] + ir*(rLim[1] - rLim[0])/nr;
		}
		rcount[ir]		=	0;
	}
	rbinL[nr]		=	rLim[1];


	// The coordinate area of each radial bin //
	// integral of sqrt(grr*gpp)dr dp //
	for( int ir=0 ; ir<nr ; ir++ ){
		if( islog ){
			rbinc[ir]	=	exp( (log(rbinL[ir+1]) + log(rbinL[ir]))/2.0 );
		}else{
			rbinc[ir]	=	(rbinL[ir+1] + rbinL[ir])/2.0;
		}
		// -------- 2*pi*sqrt[grr*gpp]dr --------- //
		area[ir]		=	2*M_PI*proper_disk_area( rbinc[ir] , a , rms )*( rbinL[ir+1]-rbinL[ir] );
	}

	for( int ip=0 ; ip<=nphi ; ip++ ){
		phiL[ip]	=	(M_PI*2)*ip/nphi;
	}
}
/*====================================*/


/******** Proper area of disk element **********/
/* actual area is proper_disk_area(r)*dr*dphi */
double disk::proper_disk_area( double r , double a , double rms ){
	double			area,a2=a*a,v_disk[4],omega,src[4],drdt[3],v_lnrf[4],gamma;
	area		=	disk_area_int_func( r , &a2 );


	// ---- Getting the area boost factor as seen by lnrf observer ---- //
	disk_velocity( r , v_disk , a , rms );

	omega			=	2*a*r/( (r*r+a2)*(r*r+a2) - a2*(r*r-2*r+a2) );
	src[0] 			= 	0; src[1] = r; src[2] = M_PI/2.; src[3] = 0;
	drdt[0] 		= 	0; drdt[1] = 0; drdt[2] = omega;
	tetrad			lnrf( src , drdt , a );

	// V in lnrf frame //
	for( int i=0 ; i<4 ; i++ ){
		v_lnrf[i]	=	0;
		for( int j=0 ; j<4 ;j++ ){ v_lnrf[i] +=	lnrf(j,i,0) * v_disk[j];}
	}
	gamma		=	1/sqrt(1-((v_lnrf[1]*v_lnrf[1])/(v_lnrf[0]*v_lnrf[0])+
							  (v_lnrf[3]*v_lnrf[3])/(v_lnrf[0]*v_lnrf[0]) ));
	area		=	area*gamma;
	return area;
}
/*==============================================*/


/********* Disk velocity inside and outside the isco *********/
void disk::disk_velocity( double in_r , double* vel , double a , double rms ){

	double		B,r,r2,a2=a*a,E,L,D,S,A;
	r		=	( in_r<rms )? rms:in_r;
	B		=	sqrt( (3/r + 2*a*pow(r,-1.5) - 1) / ( 4*a2 - 9*r + 6*r*r - r*r*r ) );

	vel[2]	=	0;
	if( in_r<rms ){
		E		=	(a - 2*pow(r,0.5) + pow(r,1.5)) * B;
		L		=	(a2 - 2*a*pow(r,0.5) + r*r ) * B;


		r		=	in_r;
		r2		=	r*r;
		D		=	r2 - 2*r + a2;
		S		=	r2;
		A		=	(r2+a2)*(r2+a2) - a2*D;
		vel[0]	=	(A*E - 2*a*r*L)/(S*D);
		vel[3]	=	(2*a*r*E + L*(D-a2) )/(S*D);

		vel[1]	=	-sqrt((D/S) * fabs( -1 + E*vel[0] - L*vel[3]));

	}else{
		vel[0]	=	(a+pow(r,1.5)) * B;
		vel[1]	=	0;
		vel[3]	=	B;
	}
}
/*===========================================================*/


/************* Calculate disk emissivity **********/
void disk::emissivity(){

	double		r,th;
	double		p_at_disk[4],v_disk[4];
	int			ir;

	// Count the number of photons in each radial bin //
	for ( int n=0 ; n<nph ; n++ ){
		r	=	data[n*ncol + 2 ];
		th	=	data[n*ncol + 3 ];

		if( fabs( th-M_PI/2. )<1e-2 ){
			if( islog ){
				ir	=	floor(log(r/rLim[0])*nr/log(rLim[1] / rLim[0]));
			}else{
				ir	=	floor((r-rLim[0])*nr/(rLim[1] - rLim[0]));
			}
			if( ir >= 0 and ir<nr){
				rcount[ir]++;
			}else{
				continue;
			}

			// photon momentum at source (=1) and at disk to get redshift factor //

			// photon Momentum at disk //
			for( int i=0 ; i<4 ; i++ ) p_at_disk[i] = data[n*ncol+5+i];
			disk_velocity( r , v_disk , a , rms );

			// get redshift //
			double		dum;
			// p_dot_v at source is -1
			dum		=	-p_dot_v( r , th , p_at_disk , v_disk , a );
			redshift[ir]	+=	dum;
		}
	}
	for( ir=0;ir<nr;ir++){
		redshift[ir]	/=	rcount[ir];
		emiss[ir]		=	rcount[ir] * redshift[ir] * redshift[ir] / area[ir];
		printf("%g %g\n",rbinc[ir],emiss[ir]);
	}
}

/***** Calculate photon momentum at a given position given its consts of motion *****/
void disk::mtm_at_position(double E,double L, double Q, double* pos , double* rdot ){
	double		r,r2,a2,D,A,S,sin2,cos2,thdot2,rdot2;
	r		=	pos[1]; r2 = r*r; a2 = a*a;
	sin2	=	sin(pos[2])*sin(pos[2]);
	cos2	=	1. - sin2;
	S		=	r2 + a2*cos2;
	D		=	r2 - 2*r + a2;
	A		=	(r2+a2)*(r2+a2) - a2*D*sin2;

	rdot[0]		=	(A*E - 2*a*r*L )/(S*D);
	rdot[3]		=	( 2*a*r*E + L*(D/sin2 - a2) ) / (S*D);
	thdot2		=	( Q - cos2*(a2*(-E*E) + L*L/sin2) ) / (S*S);
	rdot2		=	(D/S) * ( E*rdot[0] - L*rdot[3] - S*thdot2 );
	rdot[1]		=	sqrt(rdot2);
	rdot[2]		=	sqrt(thdot2);
}
/*==================================================================================*/


/***** Porject mtm onto a frame moving at v *****/
double disk::p_dot_v( double r, double th , double* p , double* v , double a ){
	double		S,D,r2,a2,sin2,gtt,gtp,grr,gthth,gpp,dum1;
	r2 = r*r; a2 = a*a; sin2 = sin(th)*sin(th);
	S			=	r2+a2*(1-sin2); D = r2-2*r+a2;
	gtt			=	-(1-2*r/S);
	gtp			=	-(2*a*r*sin2/S);
	grr			=	S/D;
	gthth		=	S;
	gpp			=	(r2+a2+2*a2*r*sin2/S)*sin2;

	dum1		=	gtt * p[0] *v[0]  + gtp*p[0]*v[3] + gtp*p[3]*v[0] +
					gpp * p[3]*v[3] + grr * p[1]*v[1] + gthth * p[2]*v[2];
	return dum1;
}


/***** TF from source to disk and to observer *****/
void disk::tf(){
	int		npix 		= 	300;
	int		nhist		=	50;
	double	lim[4]		=	{ 0.5 , 2 , 0 , 200 };
	double	imsize		=	200.;
	double	theta_o		=	1.4;


	// read flash data //
	double		r,th,phi,redshift,tsource = 0;
	double		p_at_disk[4],v_disk[4],*Area,*Count,*Redshift,*Time;
	int			it=0;
	gsl_histogram2d		*r_phi_count 		= gsl_histogram2d_alloc( nr , nphi );
	gsl_histogram2d		*r_phi_redshift 	= gsl_histogram2d_alloc( nr , nphi );
	gsl_histogram2d		*r_phi_time 		= gsl_histogram2d_alloc( nr , nphi );

	Area		=	new double[nr*nphi];
	Count		=	new double[nr*nphi];
	Redshift	=	new double[nr*nphi];
	Time		=	new double[nr*nphi];

	gsl_histogram2d_set_ranges ( r_phi_count, rbinL, nr+1, phiL, nphi+1);
	gsl_histogram2d_set_ranges ( r_phi_redshift, rbinL, nr+1, phiL, nphi+1);
	gsl_histogram2d_set_ranges ( r_phi_time, rbinL, nr+1, phiL, nphi+1);

	// Count the number of photons in each radial bin //
	for ( int n=0 ; n<nph ; n++ ){
		r	=	data[n*ncol + 2 ];
		th	=	data[n*ncol + 3 ];
		phi	=	data[n*ncol + 4 ];
		if( fabs( th-M_PI/2. )<1e-2 and r>=rbinL[0] and r<=rbinL[nr] ){
			gsl_histogram2d_increment ( r_phi_count , r , phi );

			// photon momentum at source (=1) and at disk to get redshift factor //
			// photon Momentum at disk //
			for( int i=0 ; i<4 ; i++ ) p_at_disk[i] = data[n*ncol+5+i];
			disk_velocity( r , v_disk , a , rms );

			// get redshift //
			// p_dot_v at source is -1
			redshift		=	-p_dot_v( r , th , p_at_disk , v_disk , a );
			gsl_histogram2d_accumulate (r_phi_redshift , r , phi , redshift );
			gsl_histogram2d_accumulate (r_phi_time , r , phi , data[n*ncol + 1] );
		}
		if( r == 1000 ){
			tsource += data[ n*ncol + 1 ];
			it++;
		}
	}
	tsource		/=	it;


	double		lo,hi,dr,dphi,count,dA;
	for( int i=0 ; i<nr ; i++ ){
		gsl_histogram2d_get_xrange( r_phi_count , i , &lo , &hi );
		r		=	(hi+lo)/2;
		dr		=	hi-lo;
		for( int j=0 ; j<nphi ; j++ ){
			gsl_histogram2d_get_yrange( r_phi_count , j , &lo , &hi );
			dphi	=	hi-lo;

			count		=	gsl_histogram2d_get( r_phi_count , i , j );
			redshift	=	gsl_histogram2d_get( r_phi_redshift , i , j ) / count;
			dA			=	proper_disk_area( r , a , rms ) * dr * dphi;
			Area[i*nr + j  ]	=	dA;
			Count[i*nr + j ]	=	count;
			Redshift[i*nr + j ]	=	redshift;
			Time[i*nr + j]		=	gsl_histogram2d_get( r_phi_time , i , j ) / count;
		}
	}
	// END FLASH SECTION //



	// IMAGE SECTION //
	int		i,j,nside,idx;
	long unsigned int	iir,iip;
	double	unit,x,y, out[4],flux,illum_flux,g,time;
	gsl_histogram2d		*g_t	=	gsl_histogram2d_alloc( nhist,nhist );
	gsl_histogram2d_set_ranges_uniform ( g_t, lim[0], lim[1], lim[2], lim[3]);
	image		im( a , theta_o , imsize );


	/* npix should be odd to have symmetry */
	if( npix%2 == 0) npix++;
	nside	=	(npix-1)/2;

	unit	=	imsize/npix;
	for( i=-nside ; i<=nside ; i++ ){
		x	=	unit * i;
		for( j=-nside ; j<=nside ; j++ ){
			y	=	unit * j;
			im.proj_xy( x , y , out );
			r		=	out[1]; phi = out[2]; g = out[3]; time = out[0];

			if(r<rbinL[0] or r>rbinL[nr]) continue;
			gsl_histogram2d_find ( r_phi_count , r, phi, &iir, &iip );
			idx			=	iir*nr + iip;
			dA			=	Area[ idx ];
			redshift	=	Redshift[ idx ];
			count		=	Count[ idx ];
			illum_flux	=	count * redshift * redshift / dA ;

			flux		=	illum_flux * pow(g,3)/dA;
			gsl_histogram2d_accumulate ( g_t, g, time+Time[ idx ] - tsource, flux);
		}
	}

	printf("ycent ");for( i=0 ; i<nhist ; i++ ){gsl_histogram2d_get_xrange ( g_t, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	printf("xcent ");for( i=0 ; i<nhist ; i++ ){gsl_histogram2d_get_yrange ( g_t, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	for( i=0 ; i<nhist ; i++ ){
		for( j=0 ; j<nhist ; j++ ){
			printf("%g ",gsl_histogram2d_get ( g_t, i, j ));
		}
		printf("\n");
	}
	gsl_histogram2d_free(g_t);
	gsl_histogram2d_free(r_phi_count);
	gsl_histogram2d_free(r_phi_redshift);
	gsl_histogram2d_free(r_phi_time);
	delete[] Area;delete[] Count;delete[] Redshift;delete[] Time;
}
/*================================================*/

} /* namespace gr */
