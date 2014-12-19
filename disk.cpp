/*
 * disk.cpp
 *
 *  Created on: Dec 19, 2014
 *      Author: abzoghbi
 */

#include "inc/disk.hpp"

namespace gr {

disk::disk( const string fname , int nr_ , double* rlim , int np ,bool rlog_ ) {
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
	ph_count	=	gsl_histogram2d_alloc( nr , nphi );
	redshift	=	gsl_histogram2d_alloc( nr , nphi );
	gsl_histogram2d_set_ranges ( ph_count, rL, nr+1, phiL, nphi+1);
	gsl_histogram2d_set_ranges ( redshift, rL, nr+1, phiL, nphi+1);


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
			area[ir][ip]	= 	dum * dphi;;
		}
	}
}

disk::~disk() {
	delete[] rL; delete[] phiL;
	gsl_histogram2d_free( ph_count );
	gsl_histogram2d_free( redshift );
	for( int ir=0;ir<nr;ir++ ) {delete[] area[ir];} delete[] area;
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


/************* Calculate disk emissivity **********/
void disk::emissivity(){

	double		r,th,phi,dum;
	double		p_at_disk[4],v_disk[4];
	//int			ir;

	// Count the number of photons in each radial bin //
	for ( int n=0 ; n<nph ; n++ ){
		r	=	data[n*ncol + 2 ];
		th	=	data[n*ncol + 3 ];
		phi	=	data[n*ncol + 4 ];

		if( cos( th )<1e-2 ){
			gsl_histogram2d_increment ( ph_count , r , phi );


			// photon mtm at source (=1) and at disk to get redshift factor //
			// photon mtm at disk //
			for( int i=0 ; i<4 ; i++ ) p_at_disk[i] = data[n*ncol+5+i];
			disk_velocity( r , v_disk , a , rms );
			// p_dot_v at source is -1
			dum		=	-p_dot_v( r , th , p_at_disk , v_disk , a );
			gsl_histogram2d_accumulate ( redshift , r , phi , dum );
		}
	}

	double		count,emiss;
	for( int ir=0;ir<nr;ir++){for( int ip=0 ; ip<nphi ; ip++ ){
		r			=	(rL[ir+1] + rL[ir+1])/2;
		count		=	gsl_histogram2d_get( ph_count , ir , ip );
		dum			=	gsl_histogram2d_get( redshift , ir , ip ) / count;
		emiss		=	count * dum * dum / area[ir][ip];
		printf("%g %g\n",r,emiss);
	}}

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
void disk::tf(){

	// read flash data //
	double		r,th,phi,dum,tsource = 0;
	double		p_at_disk[4],v_disk[4];
	int			it=0;
	gsl_histogram2d		*r_phi_time = gsl_histogram2d_alloc( nr , nphi );
	gsl_histogram2d_set_ranges ( r_phi_time, rL, nr+1, phiL, nphi+1);


	// Count the number of photons in each radial bin //
	for ( int n=0 ; n<nph ; n++ ){
		r	=	data[n*ncol + 2 ];
		th	=	data[n*ncol + 3 ];
		phi	=	data[n*ncol + 4 ];
		if( cos( th )<1e-2 ){
			gsl_histogram2d_increment ( ph_count , r , phi );

			// photon mtm at source (=1) and at disk to get redshift factor //
			// photon mtm at disk //
			for( int i=0 ; i<4 ; i++ ) p_at_disk[i] = data[n*ncol+5+i];
			disk_velocity( r , v_disk , a , rms );
			// p_dot_v at source is -1
			dum		=	-p_dot_v( r , th , p_at_disk , v_disk , a );
			gsl_histogram2d_accumulate ( redshift , r , phi , dum );
			gsl_histogram2d_accumulate ( r_phi_time , r , phi , data[n*ncol + 1] );
		}

		if( r == 1000 ){
			tsource += data[ n*ncol + 1 ];
			it++;
		}
	}
	tsource		/=	it;


	double		count,**illum_flux,**src_time;
	illum_flux		=	new double*[nr];
	src_time		=	new double*[nr];
	for( int ir=0;ir<nr;ir++){
		illum_flux[ir]	=	new double[nphi];
		src_time[ir]	=	new double[nphi];
		for( int ip=0 ; ip<nphi ; ip++ ){
			r			=	(rL[ir+1] + rL[ir+1])/2;
			count		=	gsl_histogram2d_get( ph_count , ir , ip );
			dum			=	gsl_histogram2d_get( redshift , ir , ip ) / count;
			illum_flux[ir][ip]		=	count * dum * dum / area[ir][ip];
			src_time[ir][ip]		=	gsl_histogram2d_get( r_phi_time , ir , ip ) / count;
		}
	}
	// END FLASH SECTION //



	// -------------------------------------------------- //
	// -------------------------------------------------- //
	// -------------------------------------------------- //

	// IMAGE SECTION //

	int		npix 		= 	300;
	int		ntime		=	100;
	int		nenergy		=	100;
	double	tLim[2]		=	{ 0 , 200 };
	double	enLim[2]	=	{ 0.5 , 2 };
	double	imsize		=	200.;
	double	theta_o		=	1.4;

	gsl_histogram2d		*tf	=	gsl_histogram2d_alloc( ntime , nenergy );
	gsl_histogram2d_set_ranges_uniform ( tf, tLim[0], tLim[1], enLim[0], enLim[1]);

	// Location of the observer and relevant variables //
	double 		rvec[4];
	rvec[0]		=	0; rvec[1] = 1000; rvec[2] = theta_o ; rvec[3] = 0;

	// Image dimensions //
	if( npix%2 == 0) npix++;
	double		**im2disk = new double*[npix*npix];
	for( int i=0;i<(npix*npix);i++) im2disk[i] = new double[4];
	// convert from pixels at x,y to t,r,phi,g //
	image2disk( npix , imsize , rvec , a , rms ,im2disk );

	unsigned long 	iir,iip;
	double			t,g;
	for( int ix =0 ; ix<(npix*npix) ; ix++ ){

		t		=	im2disk[ ix ][0];
		r		=	im2disk[ ix ][1];
		phi		=	im2disk[ ix ][2];
		g		=	im2disk[ ix ][3];


		if(r<rL[0] or r>rL[nr]) continue;
		gsl_histogram2d_find ( ph_count , r, phi, &iir, &iip );
		dum		=	illum_flux[iir][iip] * pow(g,3)/area[iir][iip];
		gsl_histogram2d_accumulate ( tf, t + src_time[iir][iip] - tsource, g , dum );
	}


	// PRINT //
	double	lo,hi;
	printf("xcent ");for( int i=0 ; i<ntime   ; i++ ){gsl_histogram2d_get_xrange ( tf, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	printf("ycent ");for( int i=0 ; i<nenergy ; i++ ){gsl_histogram2d_get_yrange ( tf, i, &lo, &hi);printf("%g ",(lo+hi)/2);}printf("\n");
	for( int ie=0 ; ie<nenergy ; ie++ ){
		for( int it=0 ; it<ntime ; it++ ){
			printf("%g ",gsl_histogram2d_get ( tf, it, ie ));
		}
		printf("\n");
	}

	gsl_histogram2d_free(r_phi_time);
	for( int ir=0;ir<nr;ir++){delete[] illum_flux[ir]; delete[] src_time[ir];}
	delete[] illum_flux;delete[] src_time;
	for( int i=0;i<(npix*npix);i++) delete[] im2disk[i];
	delete[] im2disk;
	gsl_histogram2d_free( tf );
}
/*================================================*/


/******** Project image frame to disk  **********/
void disk::image2disk( int npix , double imsize ,double *rvec , double a , double rms, double**ret ){

	int			nside;
	double		unit,x,y,out[4];

	nside	=	(npix-1)/2;
	unit	=	imsize/npix;

	for( int i=-nside, ii=0 ; i<=nside ; i++ ){
		x	=	unit * i;
		for( int j=-nside, jj=0 ; j<=nside ; j++ ){
			y	=	unit * j;
			proj_xy( x , y , rvec , a , rms , out );
			for( int n=0;n<4;n++) ret[ii*npix+jj][n] = out[n];
			jj++;
		}
		ii++;
	}
}
/*=================================================*/


/******** Project a photon from x,y on image r,phi,g,t on disk **********/
void disk::proj_xy( double x , double y , double *rvec , double a, double rms, double* out ){
	double			E,L,Q,vel[4],g,sino,sino2,coso2;
	int				thsign;

	if(fabs(x)<1e-6) { x = (x<0)?-1e-6:1e-6;}
	if(fabs(y)<1e-6) { y = (y<0)?-1e-6:1e-6;}

	sino	=	sin( rvec[2] ); sino2 = sino*sino; coso2 = 1-sino2;


	E		=	1.0;
	L		=	x*sino;
	Q		=	y*y + coso2*( -a*a*E*E + L*L/sino2);
	thsign	=	(y<0)?1:-1;
	photon		ph( E , L , Q , a , -1 , thsign );
	ph.propagate( rvec );
	if( ph.rh_stop ){
		g	=	0;
	}else{
		disk_velocity( ph.rvec[1] , vel , a , rms );
		g		=	-disk::p_dot_v( ph.rvec[1] , ph.rvec[2] , ph.rdot , vel , a );
	}
	out[0]	=	ph.rvec[0];		// t
	out[1]	=	ph.rvec[1];		// r
	out[2]	=	ph.rvec[3];		// phi
	out[3]	=	g;				// g
}
/*=========================================================================*/


} /* namespace gr */
