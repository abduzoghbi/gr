/*
 * disk.cpp
 *
 *  Created on: Dec 19, 2014
 *      Author: abzoghbi
 */

#include "inc/disk.hpp"

namespace gr {

disk::disk( const string flashfile , int nr_ , double* rlim , int np ,bool rlog_ ) {
	// Read the flash file //
	int		dim[2];
	flash::read_hdf5( flashfile , data , attr , dim );
	nph		=	dim[0];
	ncol	=	dim[1];
	a		=	attr[8];


	// ISCO and Horizon radii //
	double		z1,z2;
	z1		=	(1 + pow(1-a*a,1./3)*( pow(1+a,1/3) + pow(1-a,1/3) ));
	z2		=	sqrt( 3*a*a + z1*z1 );
	rms		=	3 + z2	- sqrt( (3-z1)*(3+z1+2*z2) );
	rh		=	1 + sqrt(1-a*a);


	// ----- r bins ------- //
	nr		=	nr_;
	rL		=	new double[nr+1];
	if(rlim[0]<rms) rlim[0] = rms;
	rlog	=	rlog_;

	for( int ir=0 ; ir<=nr ; ir++ ){
		if( rlog ){
			rL[ir]	=	exp(log(rlim[0]) + ir*(log(rlim[1] / rlim[0])	)/nr);
		}else{
			rL[ir]	=	rlim[0] + ir*(rlim[1] - rlim[0])/nr;
		}
	}

	// ------ phi bins -------- //
	nphi = np;
	phiL	=	new double[nphi+1];
	for( int ip=0 ; ip<=nphi ; ip++ ){
		phiL[ip]	=	(2*M_PI)*ip/nphi;
	}


	// ------- containers -------- //
	r_phi_count	=	gsl_histogram2d_alloc( nr , nphi );
	r_phi_redshift	=	gsl_histogram2d_alloc( nr , nphi );
	r_phi_time		=	gsl_histogram2d_alloc( nr , nphi );
	gsl_histogram2d_set_ranges ( r_phi_count, rL, nr+1, phiL, nphi+1);
	gsl_histogram2d_set_ranges ( r_phi_redshift, rL, nr+1, phiL, nphi+1);
	gsl_histogram2d_set_ranges ( r_phi_time, rL, nr+1, phiL, nphi+1);

	illum_flux		=	new double*[ nr ];
	src_time		=	new double*[ nr ];
	for( int ir = 0 ; ir<nr ; ir++ ){
		illum_flux[ir]	=	new double[nphi];
		src_time[ir]	=	new double[nphi];
	}


	// ------ Proper area of disk elements ------ //
	area		=	new double*[nr];
	double		dum,r,dr,dphi;
	for( int ir = 0 ; ir<nr ; ir++ ){
		r		=	(rL[ir+1] + rL[ir])/2.0;
		dr		=	(rL[ir+1] - rL[ir]);
		dum		=	proper_disk_area( r , a , rms ) * dr ;
		area[ir] = new double[nphi];
		for( int ip = 0 ; ip<nphi ; ip++ ){
			dphi			=	(phiL[ip+1] - phiL[ip]);
			area[ir][ip]	= 	dum * dphi;
		}
	}
}

disk::~disk() {
	delete[] rL; delete[] phiL;
	gsl_histogram2d_free( r_phi_count );
	gsl_histogram2d_free( r_phi_redshift );
	gsl_histogram2d_free( r_phi_time );
	for( int ir=0;ir<nr;ir++ ) {delete[] area[ir];delete[] illum_flux[ir]; delete[] src_time[ir];}
	delete[] area;delete[] illum_flux;delete[] src_time;
	delete[] data; delete[] attr;
}



/******** Proper area of disk element, sqrt(grr*gpp) **********/
/******** NOTE: dr and dphi are to multiplied outside *********/
double disk::proper_disk_area( double r , double a , double rms ){
	double			area,a2=a*a,r2=r*r,v_disk[4],omega,src[4],drdt[3],v_lnrf[4],gamma;
	area		=	sqrt( (r2/(r2-2*r+a2)) * (r2+a2+2*a2/r) );;


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
/*=================================================================*/






/********* Disk velocity inside and outside the isco *********/
void disk::disk_velocity( double in_r , double* vel , double a , double rms ){
	// Keplarian outside rms and ballisitc fall inside taking constant of
	// motion at the rms
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


/********* bins flash data into r_phi bins **********/
void disk::flash_to_r_phi( double& tsource ){
	double		r,th,phi,dum,count;
	double		p_at_disk[4],v_disk[4];
	int			it=0;
	tsource 	= 0;


	// Count the number of photons in each radial bin //
	for ( int n=0 ; n<nph ; n++ ){
		r	=	data[n*ncol + 1 ];
		th	=	data[n*ncol + 2 ];
		phi	=	data[n*ncol + 3 ];
		if( cos( th )<1e-2 ){
			gsl_histogram2d_increment ( r_phi_count , r , phi );

			// photon mtm at source (=1) and at disk to get redshift factor //
			// photon mtm at disk //
			for( int i=0 ; i<4 ; i++ ) p_at_disk[i] = data[n*ncol+4+i];
			disk_velocity( r , v_disk , a , rms );
			// p_dot_v at source is -1
			dum		=	-p_dot_v( r , th , p_at_disk , v_disk , a );
			gsl_histogram2d_accumulate ( r_phi_redshift , r , phi , dum );
			gsl_histogram2d_accumulate ( r_phi_time , r , phi , data[n*ncol ] );
		}

		if( r == 1000 ){
			tsource += data[ n*ncol ];
			it++;
		}
	}
	tsource		/=	it;

	for( int ir=0;ir<nr;ir++){
		r				=	(rL[ir+1] + rL[ir+1])/2;
		for( int ip=0 ; ip<nphi ; ip++ ){
			count		=	gsl_histogram2d_get( r_phi_count , ir , ip );
			dum			=	gsl_histogram2d_get( r_phi_redshift , ir , ip ) / count;
			illum_flux[ir][ip]		=	count * dum * dum / area[ir][ip];
			src_time[ir][ip]		=	gsl_histogram2d_get( r_phi_time , ir , ip ) / count;
		}
	}
}
/*==================================================*/


/************* Calculate disk emissivity **********/
void disk::emissivity( bool avg_phi ){

	double		dum;
	flash_to_r_phi( dum );

	double		r;
	if( avg_phi ){
		for( int ir=0;ir<nr;ir++){
			r		=	(rL[ir+1] + rL[ir+1])/2;
			dum		=	0;
			for( int ip=0 ; ip<nphi ; ip++ ){ dum += illum_flux[ir][ip];}
			printf("%g %g\n",r,dum);
		}
	}else{
		for( int ir=0;ir<nr;ir++){for( int ip=0 ; ip<nphi ; ip++ ){
			r		=	(rL[ir+1] + rL[ir+1])/2;
			printf("%g %g\n",r,illum_flux[ir][ip]);
		}}
	}

}
/*===================================================*/

/***** Project mtm onto a frame moving at v *****/
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
/*===============================================*/


/***** TF from source to disk and to observer *****/
void disk::tf( const string image_file, int ntime , double* tLim, int nenergy , double* enLim ){

	// SETUP FLASH: Read emissivity -> illum_flux and src-disk time -> src_time //
	double		tsource,r,dum,phi;
	flash_to_r_phi( tsource );


	// -------------------------------------------------- //
	// -------------------------------------------------- //
	// -------------------------------------------------- //

	// IMAGE SECTION //

	gsl_histogram2d		*tf	=	gsl_histogram2d_alloc( ntime , nenergy );
	gsl_histogram2d_set_ranges_uniform ( tf, tLim[0], tLim[1], enLim[0], enLim[1]);

	// Image to disk //
	image		im2disk( image_file );
	int	npix	=	im2disk.npix;


	// Image dimensions //
	unsigned long 	iir,iip;
	double			t,g;
	for( int ii=0 ; ii<npix ; ii++ ){for(int jj=0 ; jj<npix ; jj++ ){
		t		=	im2disk(ii,jj,0);
		r		=	im2disk(ii,jj,1);
		phi		=	im2disk(ii,jj,2);
		g		=	im2disk(ii,jj,3);

		if(r<rL[0] or r>rL[nr]) continue;
		gsl_histogram2d_find ( r_phi_count , r, phi, &iir, &iip );
		dum		=	illum_flux[iir][iip] * pow(g,3)/area[iir][iip];
		gsl_histogram2d_accumulate ( tf, t + src_time[iir][iip] - tsource, g , dum );
	}}


	// ------- Print tf ------ //
	double	lo,hi;
	printf("# time: ");for( int i=0 ; i<ntime   ; i++ ){gsl_histogram2d_get_xrange ( tf, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	printf("# g: ");for( int i=0 ; i<nenergy ; i++ ){gsl_histogram2d_get_yrange ( tf, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	for( int ie=0 ; ie<nenergy ; ie++ ){
	for( int it=0 ; it<ntime ; it++ ){printf("%g ",gsl_histogram2d_get ( tf, it, ie ));}
		printf("\n");
	}
	gsl_histogram2d_free( tf );
}
/*================================================*/


/***** Flux from illuminated disk at each image pixel *****/
void disk::image_flux( const string image_file ){

	// SETUP FLASH: Read emissivity -> illum_flux and src-disk time -> src_time //
	double		tsource,r,dum,phi;
	flash_to_r_phi( tsource );


	// -------------------------------------------------- //

	// Image to disk //
	image		im2disk( image_file );
	int	npix	=	im2disk.npix;
	double		*image_flux;
	image_flux	=	new double[npix*npix];


	// Image dimensions //
	unsigned long 	iir,iip;
	double			g;
	for( int ii=0 ; ii<npix ; ii++ ){for(int jj=0 ; jj<npix ; jj++ ){
		r		=	im2disk(ii,jj,1);
		phi		=	im2disk(ii,jj,2);
		g		=	im2disk(ii,jj,3);

		if(r<rL[0] or r>rL[nr]) continue;
		gsl_histogram2d_find ( r_phi_count , r, phi, &iir, &iip );
		dum		=	illum_flux[iir][iip] * pow(g,3)/area[iir][iip];
		image_flux[ii*npix + jj]		=	dum;
	}}


	// ------- Print image flux ------ //
	int		nside = (npix-1)/2;
	double	unit = im2disk.imsize/npix;
	printf("# xcent: ");for( int i=-nside ; i<=nside ; i++ ){printf("%g ",unit*i);}printf("\n");
	for( int j=0 ; j<npix ; j++ ){
		for( int i=0 ; i<npix ; i++ ){printf("%g ",image_flux[i*npix+j]);}
		printf("\n");
	}

	delete[] image_flux;
}
/*================================================*/


/***** flux at each pixel as a function of time *****/
void disk::image_flux_time( const string image_file, int ntime , double* tLim ){

	// SETUP FLASH: Read emissivity -> illum_flux and src-disk time -> src_time //
	double		tsource,r,dum,phi;
	flash_to_r_phi( tsource );

	// -------------------------------------------------- //

	// IMAGE SECTION //

	gsl_histogram		*tf	=	gsl_histogram_alloc( ntime  );
	gsl_histogram_set_ranges_uniform ( tf, tLim[0], tLim[1] );

	// Image to disk //
	image		im2disk( image_file );
	int	npix	=	im2disk.npix;
	double		**image_flux;
	image_flux	=	new double*[ntime];
	for( int it=0;it<ntime;it++) {
		image_flux[it]	=	new double[npix*npix];
		for(int i=0;i<npix*npix;i++) image_flux[it][i] = 0;
	}


	// Image dimensions //
	unsigned long 	iir,iip,iit;
	double			t,g;
	for( int ii=0 ; ii<npix ; ii++ ){for(int jj=0 ; jj<npix ; jj++ ){
		t		=	im2disk(ii,jj,0);
		r		=	im2disk(ii,jj,1);
		phi		=	im2disk(ii,jj,2);
		g		=	im2disk(ii,jj,3);

		if(r<rL[0] or r>rL[nr]) continue;
		gsl_histogram2d_find ( r_phi_count , r, phi, &iir, &iip );
		t		=	t + src_time[iir][iip] - tsource;
		if( t<tLim[0] or t>tLim[1]) continue;
		gsl_histogram_find   ( tf , t, &iit );
		dum		=	illum_flux[iir][iip] * pow(g,3)/area[iir][iip];
		image_flux[iit][ii*npix + jj]	=	dum;
		gsl_histogram_accumulate ( tf, t , dum );
	}}

	// ------- Print image flux ------ //
	int		nside = (npix-1)/2;
	double	unit = im2disk.imsize/npix, tunit = (tLim[1]-tLim[0])/ntime;

	// use split -p '^xcen' out.dat out_ to split into veus
	printf("# xcent: ");for( int i=-nside ; i<=nside ; i++ ){printf("%g ",unit*i);}printf("\n");
	printf("# time: ");for( int i=0 ; i<ntime ; i++ ){printf("%g ", tLim[0]+tunit*(i+0.5));}printf("\n");
	for( int it=0 ; it<ntime ; it ++ ){
		for( int j=0 ; j<npix ; j++ ){
			for( int i=0 ; i<npix ; i++ ){printf("%g ",image_flux[it][i*npix+j]);}
			printf("\n");
		}
	}
	gsl_histogram_free( tf );
	for( int it=0;it<ntime;it++) delete[] image_flux[it];
	delete[] image_flux;
}
/*================================================*/

} /* namespace gr */
