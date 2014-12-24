/*
 * input.cpp
 *
 *  Created on: Dec 23, 2014
 *      Author: abzoghbi
 */

#include "inc/input.hpp"

namespace gr {

input::input( char* fname ) {
	// Some definitions //
	ifstream		fp( fname );
	string			line;
	stringstream	ss,ls;

	if( not fp ){
		raise_error( "I cannot read the input file" );
	}
	// Defaults //
	src_init	=	false;
	img_init	=	false;
	disk_init	=	false;
	time_init	=	false;
	en_init		=	false;
	a_init		=	false;
	incl_init	=	false;

	// Loop through the input file lines //
	while ( getline( fp , line ) ){

		// Skip empty and comment lines //
		if ( line.empty() or line[0] == '#' or line[0] == '*' ) continue;


		// Remove white space //
		line.erase( remove( line.begin(),line.end(),' ' ) , line.end() );

		if( line[0] == 'f' or line[0] == 'F' ){
			// Flash: either a file or a set of parameters
			// {src_r,src_th,src_p,drdt_r,drdt_th,drdt_p,num_ab,flash_wfile}
			// or {flash_rfile} to read //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			string		tmps;
			double		dum;

			ss.clear(); ss.str( line );
			getline( ss , tmps , ',' );
			dum		=	atof( tmps.c_str() );
			if( dum ){
				// src_position //
				src_pos[0]	=	0;
				src_pos[1]	=	dum;
				getline( ss , tmps , ',' ); src_pos[2] = atof( tmps.c_str() );
				getline( ss , tmps , ',' ); src_pos[3] = atof( tmps.c_str() );

				// src_drdt //
				for( int i=0 ; i< 3; i++ ) {getline( ss , tmps , ',' ); src_drdt[i] = atof( tmps.c_str() );}

				// flash_num: num of alpha,beta values. nph = flash_num*flash_num //
				getline( ss , tmps , ',' ); flash_num = atoi( tmps.c_str() );

				// flash_file to write to //
				getline( ss , tmps , ',' ); flash_file = tmps;
				src_init	=	true;
			}else{
				// flash_file to read //
				flash_file	=	tmps;
			}
		}
		else if( line[0] == 's' or line[0] == 'S' ){
			// spin //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			spin	=	atof( line.c_str() );
			a_init	=	true;
		}
		else if( (line[0] == 't' or line[0] == 'T') and (line[1] == 'h' or line[1] == 'H') ){
			// theta //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			theta_o		=	atof( line.c_str() );
			incl_init	=	true;
		}
		else if( line[0] == 'I' or line[0] == 'i' ){
			// Image parameters {imsize,npix,image_file} or {image.h5} file
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			string		tmps;
			double		dum;

			ss.clear(); ss.str( line );
			getline( ss , tmps , ',' );
			dum		=	atof( tmps.c_str() );
			if( dum ){
				// read params {imsize,npix}
				imsize		=	dum;
				getline( ss , tmps , ',' ); npix   = atoi( tmps.c_str() );
				getline( ss , tmps , ',' ); image_file   = tmps;
				img_init	=	true;
			}else{
				// there is image file //
				image_file	=	tmps;
			}
		}

		else if( line[0] == 'D' or line[0] == 'd' ){
			// Disk bins parameters {rl1,rl2,nr,nphi} //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			string		tmps;

			ss.clear(); ss.str( line );
			getline( ss , tmps , ',' );rL[0] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );rL[1] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );nr = atoi( tmps.c_str() );
			getline( ss , tmps , ',' );nphi = atoi( tmps.c_str() );
			disk_init	=	true;
		}

		else if( line[0] == 't' or line[0] == 'T' ){
			// time parameters {tLim1,tLim2,ntime,enL1,enL2,nenergy} //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			string		tmps;

			ss.clear(); ss.str( line );
			getline( ss , tmps , ',' );tL[0] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );tL[1] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );ntime = atoi( tmps.c_str() );
			time_init	=	true;
		}
		else if( line[0] == 'e' or line[0] == 'E' ){
			// energy parameters {enL1,enL2,nenergy} //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			string		tmps;

			ss.clear(); ss.str( line );
			getline( ss , tmps , ',' );enL[0] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );enL[1] = atof( tmps.c_str() );
			getline( ss , tmps , ',' );nenergy = atoi( tmps.c_str() );
			en_init		=	true;
		}
		else if( line[0] == 'm' or line[0] == 'M' ){
			// spin //
			line.erase( 0 , line.find("{")+1 );
			line.erase( line.find("}") , line.size()-1 );
			mode	=	atoi( line.c_str() );
		}
	}
	check_input();

}

input::~input() {}


/******** Print error message *******/
void input::raise_error( string err ){
	static char		chars[] = {'*','#',' ','\0'};
	cerr << endl;
	cerr << setfill(chars[0]) << setw(50) << chars[0] << setfill(' ') << endl;
	print_msg( cerr , "ERROR" , &(chars[0]) );
	print_msg( cerr , err , &(chars[0]) );
	cerr << setfill(chars[0]) << setw(50) << chars[0] << setfill(' ') << endl;
	cerr << endl;
	exit(1);
}
/* ================================= */

void input::print_msg( ostream& out , const string msg , char* sep){
	string	s(3,*sep);
	out << left	<< s << " " << setw(42) << msg << right << " " << s << "\n";
}


/******* Print Usage *********/
void input::print_usage(){
	static char		chars[] = {'*','#',' ','\0'};
	cerr << setfill(chars[0]) << setw(50) << chars[0] << setfill(' ') << endl;
	print_msg( cerr , "ERROR: Wrong number of input parameters" , &(chars[0]) );
	print_msg( cerr , "I need a file similar to the following:" , &(chars[0]) );
	cerr << setfill(chars[0]) << setw(50) << chars[0] << setfill(' ') << endl;
	cerr << endl;

	print_msg( cout , "Flash = file or {src[3],drdt[3],num,fname}"   , &(chars[1]) );
	print_msg( cout , "Flash = {10,1e-3,0, 0,0,1e-4, 400, illum.h5}" , &(chars[3]) );

	print_msg( cout , "spin  = {0.9}" , &(chars[3]) );
	print_msg( cout , "theta = {1.0}" , &(chars[3]) );

	print_msg( cout , "Image = file or {imsize,npix,fname}" , &(chars[1]) );
	print_msg( cout , "Image = {100.,80, image.h5}" , &(chars[3]) );

	print_msg( cout , "Disk = {rl1,rl2,nr,nphi}" , &(chars[1]) );
	print_msg( cout , "Disk = {1.3,100,60,1}" , &(chars[3]) );

	print_msg( cout , "time = {tLim1,tLim2,ntime}" , &(chars[1]) );
	print_msg( cout , "time = {0,200,20}" , &(chars[3]) );

	print_msg( cout , "energy = {enL1,enL2,nenergy}" , &(chars[1]) );
	print_msg( cout , "energy = {0.5,2,20}" , &(chars[3]) );

	print_msg( cout , "mode = {mode} 1->flash, 2->image," , &(chars[1]) );
	print_msg( cout , "     3->emissivity, 4->image_flux" , &(chars[1]) );
	print_msg( cout , "     5->image_flux_time, 6->tf" , &(chars[1]) );
	print_msg( cout , "mode = {1}" , &(chars[3]) );

	cerr << endl;
}
/*===========================*/


/********** Check input *************/
void input::check_input(){
	if( mode == 1 ){
		// flash //
		if( not src_init ) raise_error( "Flash mode but no source position given" );
		if( not a_init ) raise_error( "Flash mode but no spin given" );
	}
	if( mode == 2 ){
		// image //
		if( not img_init ) raise_error( "Image mode but no parameters given" );
		if( not a_init ) raise_error( "Image mode but no spin given" );
		if( not incl_init ) raise_error( "Image mode but no theta_o given" );
	}
	if( mode == 3 ){
		// emissivity //
		if( flash_file.empty() or src_init ){
			raise_error( "emiss(3) requires a flash file. Run mode=1." );
		}
		if( not disk_init ) raise_error( "emiss(3) requires disk parameters." );
	}
	if( mode == 4 ){
		// image_flux //
		if( flash_file.empty() or src_init ){
			raise_error( "image_flux(4) requires a flash file. Run mode=1." );
		}
		if( image_file.empty() or img_init ){
			raise_error( "image_flux(4) requires an image file. Run mode=2." );
		}
		if( not disk_init ) raise_error( "image_flux(4) requires disk params." );
	}
	if( mode == 5 ){
		// image_flux_time //
		if( flash_file.empty() or src_init ){
			raise_error( "mode 5 requires a flash file. Run mode=1." );
		}
		if( image_file.empty() or img_init ){
			raise_error( "mode 5 requires an image file. Run mode=2." );
		}
		if( not disk_init ) raise_error( "image_flux_time (5) requires disk params." );
		if( not time_init ) raise_error( "image_flux_time (5) requires time params." );
	}
	if( mode == 6 ){
		// image_flux_time //
		if( flash_file.empty() or src_init ){
			raise_error( "mode 6 requires a flash file. Run mode=1." );
		}
		if( image_file.empty() or img_init ){
			raise_error( "mode 6 requires an image file. Run mode=2." );
		}
		if( not disk_init ) raise_error( "tf (6) requires disk params." );
		if( not time_init ) raise_error( "tf (6) requires time params." );
		if( not en_init ) raise_error( "tf (6) requires energy params." );
	}
}
/*==================================*/

} /* namespace gr */
