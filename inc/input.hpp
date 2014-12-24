/*
 * input.hpp
 *
 *  Created on: Dec 23, 2014
 *      Author: abzoghbi
 */

#ifndef INC_INPUT_HPP_
#define INC_INPUT_HPP_

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace std;

namespace gr {

class input {
	// ----- Private Variables ------ //
	bool		src_init,img_init,disk_init,time_init,en_init,a_init,incl_init;

	// ----- Private Functions ------ //
	static void print_msg( ostream& out , const string msg , char* sep );
	void raise_error( string err );
	void check_input();
public:

	// ----- Public Variables ------ //
	double		src_pos[4],src_drdt[3],spin,theta_o,imsize,rL[2],tL[2],enL[2];
	int			flash_num,npix,nr,nphi,ntime,nenergy,mode;
	string		flash_file,image_file;

	// ----- Public Functions ------ //

	input( char* infile );
	virtual ~input();
	static void print_usage();
};

} /* namespace gr */

#endif /* INC_INPUT_HPP_ */
