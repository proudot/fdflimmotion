#ifndef _SPOT_H
#define _SPOT_H

#include <fstream>
#include "math.h"

#include "CImg.h"
#include <list>
#include "tinyxml2.h"

using namespace tinyxml2;
using namespace cimg_library;
using namespace std;

template<class T>
static void max_idx(CImg<T> & img, T & max, int & xmax, int & ymax){
	max=img(0,0);
	xmax=0;
	ymax=0;
	cimg_forXY(img,x,y){
		if(img(x,y) > max){
			max = img(x,y);
			xmax = x;
			ymax = y;
		}
	}
}


template<class T>
class Spot {

protected:
  int width;
  int height;
  CImg<float> center;
  CImgList<float> move_list;
	int starting_frame_idx;
  CImg<T> intensities;
  float lifetime;
	float phase;
	float amp;
	float C;
  float score; // general benchmark score, set for classification & sorting
	float sig_x;	 // TODO : adding a model class, which define the parameter of a model;
	float sig_y;

public:
  Spot();
  Spot(const Spot<T> &spot );
  Spot(float,float,bool set_frame_move=true);
  Spot(float,float,int size);
  Spot(int,CImgList<float> &);
  void print_all (ostream &) const ;
  void set_intensities_from_move (const CImg<T> &,int mode=1);
	void recenter_moves (const CImg<T> & stack);
  inline void  pop_move_list() {  move_list.pop_back();};
  int restore (ifstream & ifs);
  void save ( ofstream &of) const ;
  float distance (Spot<T> &  sp);
  void put_intensity (T value );
  void intensities_graph();
  template<class t> Spot<T> & min_dist (std::list<t> & sp_l);
  float RMSE (const Spot<T> & sp );
  float RMSE_intensities (const Spot<T> & sp );
  bool is_out(const CImg<T> & img);
	void guess_C_and_amp();
	void get_gradient_no_interp( CImg<T>&img,
								CImgList<float>& grad,
								CImg<float>& kernel,
								CImg<float>& kernel_deriv);
	void gradient( CImg<T>&img,
								CImgList<float>& grad,
								CImg<float>& kernel,
								CImg<float>& kernel_deriv);

  // Getters, setters and what-have-yous
  inline  float get_phase () const{ return phase; } ;
  inline void set_phase (const  float value ) { phase	= value;}	;
	inline void set_sig_y(float val ){sig_y=val;};
	inline float get_sig_y(){return sig_y;};
	inline void set_sig_x(float val ){sig_x=val;};
	inline float get_sig_x(){return sig_x;};
  inline void set_amp ( const float value ) { amp	= value;};
  inline  float get_amp () const{ return amp; } ;
  inline  float get_C () const{ return C; } ;
  inline void set_C (const float value ) { C	= value;};
  inline  float get_lifetime () const{ return lifetime; } ;
  inline void set_lifetime (const float value ) { lifetime	= value;}
  inline CImg<float> get_center() const { return center;};
  inline void set_center(const CImg<float> & val) { center=val;};
  inline int get_width() const {return width;};
  inline int get_height(){return height;};
  inline void set_center(const float x, const float y){center(0)=x;center(1)=y;}
	inline void set_x(float val ){center(0)=val;};
	inline void set_y(float val ){center(1)=val;};
	inline float get_x( ){return center(0);};
	inline float get_y( ){return center(1);};
	inline int get_starting_frame_idx() const {return starting_frame_idx;};
	inline void set_starting_frame_idx(int val ){starting_frame_idx=val;};



  inline void set_new_move(const CImg<float> &val){this->move_list.push_back(val);};
  inline void set_new_move(const float x, const float y ){this->move_list.push_back(CImg<float>(1,2,1,1,x,y));};
  inline CImg<float>  get_move(int i) const{ return move_list[i];} ;
  inline void  set_move(int i,const CImg<float> &val) {  move_list[i]=val;} ;
	inline void set_move(const int i , const float x, const float y ){move_list[i]=CImg<float>(1,2,1,1,x,y);};
  inline CImgList<float>  get_move_list() const{ return move_list;};
	inline void set_move_list(CImgList<float> & val ){move_list=val;};

  inline void erase_move_list() { move_list.assign(CImgList<float>::empty());};
  CImg<T>  get_intensities () const;
	CImg<T> &  get_intensities_ref () ;
  T get_intensity   (int) const ;
  void set_intensities (const CImg<T> & );
  void set_intensity (int idx, T val);
	void interpolate(const CImg<T>&img,CImg<float>& interp_window);
	CImg<float > get_values(const CImg<T> &img);
	CImg<float > get_values_no_interp(const CImg<T> &img);
	void challenge_write_spot_XML(XMLDocument & doc,XMLNode *) const ;
	void write_spot_XML(XMLDocument & doc,XMLNode *,int idx) const ;

};
// END OF HEADER

template<class T>
Spot<T>::Spot(const Spot<T> & spot) {
	height=spot.height;
	width=spot.width;
	center=spot.center;
	move_list=spot.move_list;
	intensities=spot.intensities;
	lifetime=spot.lifetime;
	phase=spot.phase;
	C=spot.C;
	amp=spot.amp;
	sig_x=spot.sig_x;
	sig_y=spot.sig_y;
	starting_frame_idx=spot.starting_frame_idx;
}
template<class T>
Spot<T>::Spot() {
	width=7;
	height=7;
	sig_x=2;
	sig_y=2;
	center.assign(1,2);
	starting_frame_idx=0;
}

template<class T>
Spot<T>::Spot(float x, float y,bool set_frame_move) {
  center=CImg<float>(1,2,1,1,x,y);
	width=7;
	height=7;
	if(set_frame_move) move_list=CImgList<float>(center);
	lifetime=0.;
	phase=0.;
	C=0;
	amp=0;
	sig_x=2;
	sig_y=2;
	starting_frame_idx=0;
}
template<class T>
Spot<T>::Spot(float x, float y, int size) {
  center=CImg<float>(1,2,1,1,x,y);
	width=size;
	height=size;
	move_list=CImgList<float>(center);
	starting_frame_idx=0;
	lifetime=0.;
	phase=0.;
	C=0;
	amp=0;
	sig_x=2;
	sig_y=2;
}
template<class T>
Spot<T>::Spot(int starting_frame_idx,CImgList<float> & move_list) {
	width=7;
	height=7;
	sig_x=2;
	sig_y=2;
	center.assign(1,2);
	set_starting_frame_idx(starting_frame_idx);
	set_move_list(move_list);
}

template < class T >
T Spot<T>::get_intensity ( int idx ) const
{
	return intensities(idx);
}		/* -----  end of method Spot<T>::get_intensities  ----- */

template < class T >
CImg<T>  Spot<T>::get_intensities () const
{
	return intensities;
}		/* -----  end of method Spot<T>::get_intensities  ----- */
template < class T >
  CImg<T> &  Spot<T>::get_intensities_ref ()
{
	return intensities;
}		/* -----  end of method Spot<T>::get_intensities  ----- */

template < class T >
  void Spot<T>::set_intensities (const CImg<T> & value )
{
	intensities	= value;
}		/* -----  end of method Spot<T>::set_intensities  ----- */

template<class T>
bool Spot<T>::is_out(const CImg<T> & img) {
  return (Spot<T>::center(0)<width)
		||(Spot<T>::center(1)<height)
		||(Spot<T>::center(1)>=(img.height()-height))
		||(Spot<T>::center(0)>=(img.width()-width));
  }
template < class T >
void Spot<T>::set_intensity (int idx,  T value )
{
	intensities(idx)	= value;
}		/* -----  end of method Spot<T>::set_intensities  ----- */

template < class T >
float Spot<T>::distance (Spot<T> & sp)
{
	return (sp.move_list[0]-move_list[0]).magnitude();
}		/* -----  end of method Spot<T>::set_intensities  ----- */

template < class T >
void Spot<T>::recenter_moves (const CImg<T> & stack)
{
	int xmax=0;
	int ymax=0;
	T max;
	CImg<T> slice;
	for(int z =0 ; z < move_list.size() ; z++){
		slice=stack.get_slice(z);
		slice.crop(move_list[z](0)-width/2,move_list[z](1)-height/2,
							 move_list[z](0)+width/2,move_list[z](1)+height/2);
		max_idx(slice,max,xmax,ymax);
		move_list[z](0)=move_list[z](0)-width/2+xmax;
		move_list[z](1)=move_list[z](1)-height/2+ymax;
	}
}		/* -----  end of method Spot<T>::set_intensities  ----- */

template < class T >
void Spot<T>::set_intensities_from_move (const CImg<T> & stack, int mode)
{
	CImg<T> slice;
	intensities.assign(move_list.size());
	cimg_forX(intensities,z){
		switch(mode){
		case 0:
			intensities(z)=stack.linear_atXY(move_list[z](0),move_list[z](1),z);
		case 1:
			slice=stack.get_slice(z);
			slice.crop(move_list[z](0)-width/2,move_list[z](1)-height/2,
								 move_list[z](0)+width/2,move_list[z](1)+height/2);
			intensities(z)=slice.max();
			break;
		case 2:
			intensities(z)=stack.cubic_atXY(move_list[z](0),move_list[z](1),z);
			break;
		default:
			cout << "unknown set_intensities_from_move mode : " << mode << endl;
		}

// 		intensities(z)=stack(move_list[z](0),move_list[z](1),z);
	}
}		/* -----  end of method Spot<T>::set_intensities  ----- */


template < class T >
int Spot<T>::restore (ifstream & ifs)
{
	bool v = true;
	if(v) cout << "Restore spot : " ;
	char test[200];
	ifs.getline(test,200); //xy
	if(ifs.eof()) return 0;
	center.assign(1,2);
	ifs >> center(0) >> center(1);

	unsigned int move_list_size;
	ifs.getline(test,200,' ');
	ifs >> move_list_size ;
	CImg<float> move(1,2);
	move_list.clear();
	for ( unsigned int i=0 ; i<move_list_size ; i++){
		ifs >>   move(0);        // x
		move_list.push_back(move);
	}
	for ( unsigned int i=0 ; i<move_list_size ; i++){
		ifs >>   move_list[i](1) ;       // x
	}
	if(v) cout << " first pos : " << move_list[0](0) << " , "
		   << move_list[0](1)  << endl;
	ifs.getline(test,200,' ');
	ifs.getline(test,200,' ');
	ifs >> move_list_size;
	if(v) cout << "intensities : " << move_list_size << endl;
	intensities.assign(move_list_size);
	T pfff=0;
	cimg_forX(intensities,x){
		ifs  >> pfff;
		intensities(x)=pfff;
		if(v) cout << " pff : " << pfff << endl;
	}
	if(move_list_size==0) ifs.getline(test,200);
	ifs.getline(test,200);
	ifs.getline(test,200);
	ifs >> lifetime;
	if(v) cout << lifetime << endl;
	ifs.getline(test,200);
	ifs.getline(test,200);
	ifs >> phase;
	if(v) cout << phase << endl;
	ifs.getline(test,200);
	ifs.getline(test,200);
	ifs >> C;
	if(v) cout << C << endl;
	ifs.getline(test,200);
	ifs.getline(test,200);
	ifs >> amp;
	if(v) cout << amp << endl;
	ifs.getline(test,200);
	ifs.getline(test,200);
	return 1;

}		/* -----  end of method Spot<T>::set_intensities  ----- */

template<class T>
void Spot<T>::interpolate(const CImg<T>&img,CImg<float>& interp_window){
	int w_s = width;
	float pos_x=0.;
	float pos_y=0.;
	//could be done in an other function, let here to keep the simetry with the
	//get_interp_and_gradient function
	interp_window.assign(width,height);
  cimg_forXY(interp_window,x,y){
		pos_x=Spot<T>::center(0)+float(x-w_s/2);
		pos_y=Spot<T>::center(1)+float(y-w_s/2);
		interp_window(x,y)=img.linear_atXY(pos_x,pos_y);
  }
}
template<class T>
CImg<float> Spot<T>::get_values(const CImg<T> &img){
	CImg<float> res;
	interpolate(img,res);
	return res;
}
template<class T>
CImg<float> Spot<T>::get_values_no_interp(const CImg<T> &img){
	CImg<float> res;
	CImg<float> rounded_center=get_center().get_round();
	CImg<float> center_bk=get_center();
	set_center(rounded_center); // "Don't try this at home" 3 lines
	interpolate(img,res);
	set_center(center_bk);
	return res;
}
template < class T >
void Spot<T>::print_all (ostream & of) const
{
	of << "xy" << endl << center(0) <<" " <<center(1) << endl;
	of << "move_list" << " " << move_list.size() << endl ;
	for (unsigned int i=0 ; i<move_list.size() ; i++){
		of <<   (get_move(i))(0)        // x
			<< " ";
	}
	of << endl;
	for (unsigned int i=0 ; i<move_list.size() ; i++){
		of <<   (get_move(i))(1)        // x
			<< " ";
	}
	of << endl;
	of << "intensities" << " "<< intensities.size() << endl;
	cimg_forX(intensities,x){
		of  << intensities(x)
					<< " ";
	}
	of<< endl;
	of << "lifetime" <<endl;
	if (this->get_lifetime() != this->get_lifetime()){
		of << 0 << endl ;         /* NaN detection trick */
	}else{
		of << this->get_lifetime() << endl ;
	}
	of << "phase" <<endl<< this->get_phase() << endl;
	of << "C" <<endl<< this->get_C() << endl;
	of << "amp" <<endl<< this->get_amp() << endl;
	of<<endl;
}		/* -----  end of method Spot<T>::set_intensities  ----- */

template < class T >
float Spot<T>::RMSE_intensities (const Spot<T> & sp )
{
	float res=0;
	cimg_forX(intensities ,x){
		res+=pow(double(intensities(x) - sp.get_intensities()(x)),2.);
	}
	res/=intensities.size();
	res=sqrt(res);
	return res;

}		/* -----  end of method Spot<T>::set_intensities  ----- */
template < class T >
float Spot<T>::RMSE (const Spot<T> & sp )
{
	float res=0;
	for (unsigned int i=0 ; i<move_list.size() ; i++){
		res+=(move_list[i] - sp.get_move(i)).pow(2).sum();
	}
	res/=move_list.size();
	res=sqrt(res);
	return res;

}

template < class T >
void Spot<T>::save (  ofstream &of) const
{

	for (unsigned int i=0 ; i<move_list.size() ; i++){
		of << "Point " << (i+1)
			<<  " " << (move_list[i])(0)        // x
			<<	" " << (move_list[i])(1)				// y
			<<  " " << 1.0                        // z
			<<  " " << (i+1)                      // frame nb
			<<  " " << 1.0                        // channel
			<< endl;
	}
}

template < class T > template<class t>
Spot<T> &  Spot<T>::min_dist (std::list<t> & sp_l)
{
	typedef std::list<t> t_spot_l;
	typename t_spot_l::iterator it=sp_l.begin();
	float min_res=distance(*it);
	float tmp;
	Spot<T> match(*it);
	++it;
	for (; it!=sp_l.end(); ++it) {
		tmp=distance(*it);
		if (min_res > tmp){
			match=*it;
			min_res=tmp;
		}
	}
	score=min_res;
	return match;

}

// Contract : at least 2 intensity in intensities
template< class T>
void Spot<T>::guess_C_and_amp(){
	static ofstream log("guess.log") ;
  static CImg<float> sys(2,3);
  static CImg<float> t_sys;
  static float Pi=3.141592;			// Pi value
  static float wt=2*Pi*0.04; 		// the Li-Flim software is usually set up to 40
  static float w=2*Pi/12;
  static float phase=atan(2.5*wt)+2.42342-atan(4.1*wt);	/* FIXME */
	set_lifetime(2.5);
  sys(0,0) = 1 ;
  sys(0,1) = 1 ;
  sys(0,2) = 1 ;
  sys(1,0) = sin(phase) ;
  sys(1,1) = sin(w+phase) ;
  sys(1,2) = sin(2*w+phase) ;
	t_sys=sys.get_transpose();


  CImg<float> B(1,3);
  B(0)=get_intensity(0);
  B(1)=get_intensity(1);
  B(2)=get_intensity(2);
  CImg<float> res;
  res = (t_sys*B).solve((t_sys*sys));
  set_C(res(0));
  set_amp(res(1));
	set_phase(phase);
	// set_sig_x(2);
	print_all(log);

}

template< class T>
void Spot<T>::intensities_graph(){
	const unsigned char white[]={255};
	CImg<unsigned char> graph(600,500);
	CImgDisplay  draw_disp(graph,"Intensity profile");
	while (!draw_disp.is_closed()) {
		draw_disp.wait();
		graph.fill(0).draw_graph(intensities,white,1,1,1,0,0,0,true);
		graph.display(draw_disp);
	}

}

template<class T>
void Spot<T>::get_gradient_no_interp( CImg<T>&img,
																						CImgList<float>& grad,
																						CImg<float>& kernel,
																						CImg<float>& kernel_deriv){
	CImg<float> rounded_center=get_center().get_round();
	CImg<float> center_bk=get_center();
	set_center(rounded_center); // "Don't try this at home" 3 lines
	gradient(img,grad,kernel,kernel_deriv);
	set_center(center_bk);
}

template<class T>
void Spot<T>::gradient( CImg<T>&img,
												CImgList<float>& grad,
												CImg<float>& kernel,
												CImg<float>& kernel_deriv){

  CImg<float> normal_crop;
	CImgList<float> grad_tmp(2);
	float pos_x=0.;
	float pos_y=0.;

  unsigned int border_size=max(kernel.width(),kernel_deriv.width())/2;
  int w= width+border_size*2;
  int h= height+border_size*2;
	CImg<float> rounded_center=get_center().get_round();
	normal_crop=img.get_crop(rounded_center(0) - w/2, rounded_center(1) - h/2,
													 rounded_center(0) + w/2, rounded_center(1) + h/2);


  //Compute gradient
	//interpolated_crop.display();
  grad_tmp[0]=normal_crop ;
	grad_tmp[0].convolve(kernel_deriv);
  grad_tmp[0].convolve(kernel.get_transpose());

  grad_tmp[1]=normal_crop ;
	grad_tmp[1].convolve(kernel);
  grad_tmp[1].convolve(kernel_deriv.get_transpose());


	float dx = get_center()(0)-rounded_center(0);
	float dy = get_center()(1)-rounded_center(1);
	grad.assign(2,width,height);
  cimg_forXY(grad[0],x,y){
		pos_x=float(x)+dx;
		pos_y=float(y)+dy;
		grad[0](x,y)=grad_tmp[0].linear_atXY(border_size+ pos_x,border_size+ pos_y);
		grad[1](x,y)=grad_tmp[1].linear_atXY(border_size+ pos_x,border_size+ pos_y);
  }

}
template<class T>
void Spot<T>::write_spot_XML(XMLDocument & doc, XMLNode * current_node, int idx )const {
	if(move_list.size()>0){
		XMLElement * particle=doc.NewElement("particle");
		particle->SetAttribute("idx",idx);
		particle->SetAttribute("lifetime",get_lifetime());
		current_node->InsertEndChild(particle);
		cimg_for1(move_list.size(),i){
			CImg<float> m=get_move(i);
			XMLElement *	detection = doc.NewElement("detection");
			detection->SetAttribute("t",starting_frame_idx+i);
			detection->SetAttribute("x",m(0));
			detection->SetAttribute("y",m(1));
			detection->SetAttribute("z",0);
			particle->InsertEndChild(detection);
		}
	}
	// current_node->InsertEndChild(particle);
}


template<class T>
void Spot<T>::challenge_write_spot_XML(XMLDocument & doc,XMLNode * current_node )const {
	if(move_list.size()>0){
		XMLElement * particle=doc.NewElement("particle");
		current_node->InsertEndChild(particle);
		cimg_for1(move_list.size(),i){
			CImg<float> m=get_move(i);
			XMLElement *	detection = doc.NewElement("detection");
			detection->SetAttribute("t",starting_frame_idx+i);
			detection->SetAttribute("x",m(0));
			detection->SetAttribute("y",m(1));
			detection->SetAttribute("z",0);
			particle->InsertEndChild(detection);
		}
	}
	// current_node->InsertEndChild(particle);
}

template< class T>
void challenge_write_spot_XML(const typename  std::list<T> & spot_l,const char * filename,
															int snr = 7, const char * density = "low" , const char * scenario="vesicle"){
	typedef std::list<T> t_spot_l;
	typename t_spot_l::const_iterator s;
	XMLDocument doc;
	doc.InsertFirstChild(doc.NewDeclaration("xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\""));
	XMLElement * root = doc.NewElement("root");
	doc.InsertEndChild(root);
	XMLElement * contestHeader  = doc.NewElement("TrackContestISBI2012");
	contestHeader->SetAttribute("snr",snr);
	contestHeader->SetAttribute("density",density);
	contestHeader->SetAttribute("scenario",scenario);
	root->InsertFirstChild(contestHeader);
	for (s=spot_l.begin(); s!=spot_l.end(); ++s){
		s->challenge_write_spot_XML(doc,contestHeader);
	}

	doc.SaveFile(filename);
}
template< class T>
void write_spot_XML(const typename  std::list<T> & spot_l,const char * filename){
	typedef std::list<T> t_spot_l;
	typename t_spot_l::const_iterator s;
	XMLDocument doc;
	doc.InsertFirstChild(doc.NewDeclaration("xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\""));
	XMLElement * root = doc.NewElement("root");
	doc.InsertEndChild(root);
	int fidx=0;
	for (s=spot_l.begin(); s!=spot_l.end(); ++s){
		s->write_spot_XML(doc,root,fidx++);
	}
	doc.SaveFile(filename);
}

template< class T>
std::list<T> restore_spots( const char * infile){
	ifstream of;
	of.open(infile);
	typedef std::list<T> t_spot_l;
	t_spot_l res;
	T sp;
	for (;!of.eof();){
		if(sp.restore(of)) res.push_back(sp);
	}
	of.close();
	return res;
}

template< class T>
CImg<float> get_array(typename std::list<T> & spot_l,float (T::*getter)() const ){
	typedef std::list<T> t_spot_l;
	typename t_spot_l::iterator s;
	CImg<float> res(spot_l.size());
	int i=0;
	for (s=spot_l.begin(); s!=spot_l.end(); ++s){
		res(i)=((*s).*getter)();
		i++;
	}
	return res;
}

template< class T,class t>
void recenter_spots(typename std::list<T> & spot_l, CImg<t> stack){
	typedef std::list<T> t_spot_l;
	typename t_spot_l::iterator s;
	for (s=spot_l.begin(); s!=spot_l.end(); ++s){
			s->recenter_moves(stack);
	}
}

template< class T>
void save_spots(typename std::list<T> & spot_l, const char * outfile,bool print_lost=true){
	ofstream of;
	of.open(outfile);
	typedef std::list<T> t_spot_l;
	typename t_spot_l::iterator s;
	for (s=spot_l.begin(); s!=spot_l.end(); ++s){
		if((print_lost)||(!s->get_lost())) s->print_all(of);
	}
	of.close();
}


template< class T>
void save_MTrackJ(typename std::list<T> & spot_l, const char * outfile,bool print_lost=false){
  ofstream of;
  of.open(outfile);
  of << "MTrackJ 1.3.1 Data File" << endl;
  of << "Assembly 1 FF0000" << endl;
  of << "Cluster 1 " << endl;
  typedef std::list<T> t_spot_l;
  typename t_spot_l::iterator s;
  int i=1;
  for (s=spot_l.begin(); s!=spot_l.end(); ++s){
//     if(!s->get_lost()||print_lost){
      of << "Track " << i << endl;
      s->save(of);
      i++;
//     }
  }
  of << "End of MTrackJ Data File" << endl;
  of.close();
}

// #include <sin_fit.hpp>
// #include <varEstimator.hpp>
// template< class T>
// void compute_lifetime_irls(typename std::list<T> & spot_l,
// 													 float phase_ref,
// 													 float tau_ref,
// 													 float err_sig,
// 													 float alpha,
// 													 float (*pt2CostFunc)(float, float,float),
// 													 bool display = false,
// 													 CImg<float> * err_sigs = NULL){
// 	typename std::list<T>::iterator s=spot_l.begin();
// 	// fill  tracked spot intensity
// 	// And fit intensities
// 	cout << "fit intensities" << endl;
// 	float Pi=3.141592;	 //  Pi value
// 	float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
// 	// Mhz
// 	ofstream log_fit("log_fit_irls.txt");
// 	float phase;
// 	float C;
// 	float amp;
// 	float lifetime;
// 	CImg<float> Y;
// 	CImg<float> beta;
// 	Sin_fit fitter;
// 	Leclerc cf_leclerc(err_sig,alpha);
// 	fitter.set_cost_function(&cf_leclerc);
// 	float phase_shift=phase_ref-atan(w*tau_ref);
// 	CImg<float> residues;
// 	for ( s = spot_l.begin(); s != spot_l.end(); ++s) {
// 		Y.assign(s->get_intensities());
// 		fitter.assign(Y);
// 		fitter.run(beta);
// 		fitter.irls_no_init(beta,alpha,pt2CostFunc);
// 		residues.append(fitter.get_residue(),'x');
// 		log_fit << "display fit : spot ( "
// 						<< s->get_center()(0) << ","
// 						<< s->get_center()(1)
// 						<< " )" << endl ;
// 		log_fit << "fit param : "
// 						<< beta(0) << " "
// 						<< beta(1) << " "
// 						<< beta(2) << " "
// 						<< beta(3) << " "
// 						<< endl ;
// 		phase=beta(2);
// 		amp=beta(1);
// 		C=beta(0);
// 		s->set_C(C);
// 		s->set_amp(amp);
// 		s->set_phase(phase);
// 		if(display) fitter.display_fit();
// 		lifetime=tan(-phase_shift+phase)/w;
// 		log_fit << "lifetime : " << lifetime << endl;
// 		if(lifetime != lifetime) lifetime=0; 		// Nan detection trick
// 		s->set_lifetime(fabs(lifetime));
// 	}
// 	log_fit.close();
// }

// template< class T,class c>
// void compute_lifetime_multiple(typename std::list<T> & spot_l,
// 															 CImg<c> stack,
// 															 float phase_ref=0,
// 															 float tau_ref=0,
// 															 int scale_parameter_type=0,
// 															 float (*pt2CostFunc)(float, float,float)=NULL,
// 															 bool display = false){
// 	typename std::list<T>::iterator s=spot_l.begin();
// 	// fill  tracked spot intensity
// 	// And fit intensities
// 	float Pi=3.141592;	 //  Pi value
// 	float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
// 	// Mhz
// 	float phase;
// 	float C;
// 	float amp;
// 	float lifetime;
// 	CImg<float> Y;
// 	CImg<float> X;
// 	int idx = 0;
// 	CImg<float> feature_value;
// 	CImg<float> beta(3);
// 	beta.fill(0);
// 	Sin_fit fitter;
// 	Leclerc cf_leclerc(0,3);
// 	fitter.set_cost_function(&cf_leclerc);

// 	float phase_shift=phase_ref-atan(w*tau_ref);

// 	for ( s = spot_l.begin(); s != spot_l.end(); ++s) {
// 		Y.assign(s->get_height()*s->get_width()*s->get_move_list().size());
// 		X.assign(s->get_height()*s->get_width()*s->get_move_list().size());
// 		idx=0;
// 		for(unsigned int t = 0; t < s->get_move_list().size();t++){
// 			s->set_center(s->get_move(t));
// 			feature_value=s->get_values(stack.get_slice(t));
// 			cimg_forXY(feature_value,x,y){
// 				Y(idx)=feature_value(x,y);
// 				X(idx)=t;
// 				idx++;
// 			}
// 		}

// 		fitter.init(X,Y,beta);
// 		fitter.set_starting_value();
// 		fitter.run(beta);

// 		if(scale_parameter_type){
// 			float err_sig=0;
// 			switch(scale_parameter_type){
// 			case 1 : {
// 				varEstimator<c> ve;
// 				err_sig = ve.estimate(stack);
// 				break;
// 			}
// 			case 2: {
// 				varEstimator<c> ve;
// 				CImg<float> r = fitter.get_residue();
// 				// r.get_histogram(10).display_graph(0,1,1,0,r.min(),r.max());
// 				ve.set_residue(r);
// 				cout << "r.median() :" << r.median() << endl;
// 				cout << "r.mean() :" << r.mean() << endl;

// 				err_sig = ve.sigmaMAD();
// 				break;
// 			}
// 			}
// 			fitter.set_scale_parameter(err_sig);
// 			fitter.irls_no_init(beta,3,pt2CostFunc);
// 		}


// 		phase=beta(2);
// 		amp=beta(1);
// 		C=beta(0);
// 		s->set_C(C);
// 		s->set_amp(amp);
// 		s->set_phase(phase);
// 		if(display) fitter.display_fit();
// 		lifetime=tan(-phase_shift+phase)/w;
// 		if(lifetime != lifetime) lifetime=0; 		// Nan detection trick
// 		s->set_lifetime(fabs(lifetime));
// 	}
// }

// template< class T>
// void compute_lifetime(typename std::list<T> & spot_l,
// 											float phase_ref=0,
// 											float tau_ref=0,
// 											bool display = false){
// 	typename std::list<T>::iterator s=spot_l.begin();
// 	// fill  tracked spot intensity
// 	// And fit intensities
// 	cout << "fit intensities" << endl;
// 	float Pi=3.141592;	 //  Pi value
// 	float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
// 	// Mhz
// 	ofstream log_fit("log_fit_normal.txt");
// 	float phase;
// 	float C;
// 	float amp;
// 	float lifetime;
// 	CImg<float> Y;
// 	CImg<float> beta;
// 	Sin_fit fitter;
// 	float phase_shift=phase_ref-atan(w*tau_ref);
// 	for ( s = spot_l.begin(); s != spot_l.end(); ++s) {
// 		Y.assign(s->get_intensities());
// 		fitter.assign(Y);
// 		fitter.run(beta);
// 		log_fit << "display fit : spot ( "
// 						<< s->get_center()(0) << ","
// 						<< s->get_center()(1)
// 						<< " )" << endl ;
// 		log_fit << "fit param : "
// 						<< beta(0) << " "
// 						<< beta(1) << " "
// 						<< beta(2) << " "
// 						<< endl ;
// 		phase=beta(2);
// 		amp=beta(1);
// 		C=beta(0);
// 		s->set_C(C);
// 		s->set_amp(amp);
// 		s->set_phase(phase);
// 		if(display) fitter.display_fit();
// 		lifetime=tan(-phase_shift+phase)/w;
// 		log_fit << "lifetime : " << lifetime << endl;
// 		if(lifetime != lifetime) lifetime=0; 		// Nan detection trick
// 		s->set_lifetime(fabs(lifetime));
// 	}
// 	log_fit.close();
// }


// template< class T>
// void compute_lifetime_fourier(typename std::list<T> & spot_l,
// 															float phase_ref=0,
// 															float tau_ref=0,
// 															bool display = false){
// 	typename std::list<T>::iterator s=spot_l.begin();
// 	// fill  tracked spot intensity
// 	// And fit intensities
// 	cout << "fit intensities" << endl;
// 	float Pi=3.141592;	 //  Pi value
// 	float w=2*Pi*0.04;	 // the Li-Flim software is usually set up to 40
// 	// Mhz
// 	int nb_frame=12;

// 	float phase;
// 	float lifetime;
// 	CImg<float> Cos(nb_frame);
// 	CImg<float> Sin(nb_frame);
// 	cimg_forX(Cos,k){
// 		Cos(k)=cos(2*Pi*(k)/(nb_frame));
// 		Sin(k)=sin(2*Pi*(k)/(nb_frame));
// 	}
// 	CImg<float> beta;
// 	Sin_fit fitter;
// 	float phase_shift=phase_ref-atan(w*tau_ref);

// 	float fcos=0;
// 	float fsin=0;
// 	float val=0;
// 	CImg<float> Y;
// 	ofstream log_fit("log_fourier_fit.txt");
// 	for ( s = spot_l.begin(); s != spot_l.end(); ++s) {
// 		Y.assign(s->get_intensities());
// 		fcos=0;
// 		fsin=0;
// 		cimg_forX(Cos,z){
// 			val=Y(z);
// 			fcos+=val*Cos(z);
// 			fsin+=val*Sin(z);
// 		}
// 		phase=-atan(fsin/fcos);
// 		phase+=Pi/2;
// 		s->set_phase(phase);

// 		lifetime=tan(phase-phase_shift)/w;
// 		if(lifetime != lifetime) lifetime=0; 		// Nan detection trick
// 		log_fit << "display fit : spot ( "
// 						<< s->get_center()(0) << ","
// 						<< s->get_center()(1)
// 						<< " )" << endl ;
// 		log_fit << "phase " << phase  <<endl;
// 		log_fit << "lifetime : " << lifetime << endl;
// 		s->set_lifetime(fabs(lifetime));
// 	}
// 	log_fit.close();
// };
#endif /* _SPOT_H */

