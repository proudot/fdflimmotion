// Purpose : Challenge ISBI 2012

#include <iostream>
#include <fstream>
#include <list>
#include "math.h"
#include <algorithm>
#include <time.h>

//My include
#include "CImg.h"
#include "lifetime.hpp"


using namespace cimg_library;
using namespace std;

int main(int argc,char** argv){

	const char * file_name  = cimg_option("-tiff",(char*)NULL,"input file");
	const char * ref_file_name  = cimg_option("-ref",(char*)NULL,"input reference stack");
	const float ref_lifetime=cimg_option("-ref_lf",4.1,"fluorescence lifetime of the reference sample"); 
	const int ref_frame_nb = cimg_option("-frame_nb",12,"nb of frame to process");
	const int frame_nb = cimg_option("-frame_nb",ref_frame_nb,"nb of frame to process");

	typedef float img_type;

	CImg<img_type> stack; stack.load_tiff(file_name);
	if(frame_nb>0) stack.slices(0,frame_nb-1);
	phase_processor pp; pp.assign(stack.width(),stack.height(),frame_nb);
	cfl ref_phases,ref_stack;
	if((ref_file_name != NULL))  {
		ref_stack.assign(ref_file_name).slices(0,ref_frame_nb-1);	
		ref_phases = pp.process_pixels(ref_stack,1,0);
		ref_phases.save("ref_phases.tiff");
	}

	cfl stack_phases=pp.process_pixels(stack,1);
	lifetime_processor lp;
	CImg<> plain_lf = lp.process_phase(stack_phases,ref_phases,ref_lifetime);
	plain_lf.save("fourier_lifetime.tiff");
}

