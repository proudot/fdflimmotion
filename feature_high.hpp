#ifndef _FEAT_HIGHLIGHT_
#define _FEAT_HIGHLIGHT_



#include <iostream>
#include "math.h"
#include <algorithm>


//My include
#include "CImg.h"
#include "sin_fit.hpp"

using namespace cimg_library;
using namespace std;
template<class T>
void draw_line(CImg<T> & img, int x0, int y0, int x1, int y1,T * color, int thickness=0){
	for(int t = -thickness; t<=thickness;t++) img.draw_line(x0+t,y0+t,x1+t,y1+t, color,1,~0U);
}

struct feat_highlight{
	feat_highlight(){
		window_size=7;
		frame_idx=0;
		scale=1;
		print_lost=false;
		color=false;
		color_map =CImg<unsigned char>::default_LUT256();
	}

	float window_size;
	unsigned int frame_idx;
	float scale;
	bool print_lost;
	bool color;
	CImg< unsigned char > color_map;

	template<class T,class t>
	CImg<T>  build_tracks(const CImg<t> & stack, typename std::list<Feature<T> > & feature_l ){
		CImg<T> feat_high(stack);
		float sc_width=scale*stack.width();
		float sc_height=scale*stack.height();
		feat_high.resize(sc_width,sc_height);

		return draw_tracks_on_img(feat_high,feature_l);
	}


	template<class T>
	CImg<T>  build_stack(const CImg<T> & stack, const typename std::list<Feature<T> > & feature_l ){
		CImg<T> feat_high(stack);
		float sc_width=scale*stack.width();
		float sc_height=scale*stack.height();
		feat_high.resize(sc_width,sc_height);

		return draw_feat_on_img(feat_high,feature_l);
	}
	template<class T>
	CImg<unsigned char>  build_stack_color(const CImg<T> & stack, const typename std::list<Feature<T> > & feature_l ){
		CImg<unsigned char> feat_high;
		float sc_width=scale*stack.width();
		float sc_height=scale*stack.height();
		CImg<unsigned char> rstack=  stack.get_normalize(0,255).get_resize(sc_width,sc_height);
		feat_high.append(rstack,'c').append(rstack,'c').append(rstack,'c');
		// feat_high=stack.get_normalize(0,255).get_resize(sc_width,sc_height,1,3);


		return draw_feat_on_img(feat_high,feature_l);
	}

	template<class T,class t>
	CImg<t>  draw_tracks_on_img(CImg<t> & feat_high, typename std::list<Feature<T> > & feature_l ){
		typedef  std::list<Feature<T> > t_feature_l;
		typename  t_feature_l::const_iterator s;
			// current.resize_tripleXY();
		for (s=feature_l.begin(); s!=feature_l.end(); ++s){
				if(!s->get_lost() || (print_lost)) {
					build_feature_tracks<T,t>(*s,feat_high);
				}
			}
		return feat_high;
	}

	template<class T,class t>

	CImg<t>  draw_f_idx_on_img(CImg<t> & feat_high, typename std::list<Feature<T> > & feature_l ){
		typedef  std::list<Feature<T> > t_feature_l;
		typename  t_feature_l::const_iterator s;
		uint feature_idx = 0;
		for (s=feature_l.begin(); s!=feature_l.end(); ++s,feature_idx++){
				if(!s->get_lost() || (print_lost)) {
					int color_index=cimg::rand()*255;
					uint drawn_stack_idx=s->get_starting_frame_idx();
					uint move_idx=0;
					for(;move_idx<s->get_move_list().size();move_idx++, drawn_stack_idx++){
						cimg_forC(feat_high,c){
							char feature_idx_str[2];
							sprintf(feature_idx_str,"%d",feature_idx);
							feat_high.get_shared_plane(drawn_stack_idx,c)
								.draw_text((s->get_move(move_idx)(0)+s->get_width())*scale,
													 (s->get_move(move_idx)(1)+s->get_width())*scale,
													 feature_idx_str,
													 color_map.data(255,1,1,c),
													 color_map.data(0,1,1,c));
							}
						}
					}
		}

		return feat_high;
	}

	template<class T,class t>
	CImg<t>  draw_feat_on_img(CImg<t> & feat_high, const typename std::list<Feature<T> > & feature_l ){
		typedef  std::list<Feature<T> > t_feature_l;
		typename  t_feature_l::const_iterator s;
			// current.resize_tripleXY();
		for (s=feature_l.begin(); s!=feature_l.end(); ++s){
				if(!s->get_lost() || (print_lost)) {
					build_feature_highlight<T,t>(*s,feat_high);
				}
			}
		return feat_high;
	}

	template <class T,class t>
	void  build_feature_tracks(const Feature<T> & f,  CImg<t> & stack){
		int color_index=cimg::rand()*255;
		cimg_for1(f.get_move_list().size()-1,fr_idx){
			uint drawn_stack_idx=f.get_starting_frame_idx();
			for(;drawn_stack_idx<f.get_starting_frame_idx() +f.get_move_list().size();
					drawn_stack_idx++){
				cimg_forC(stack,c){
					stack.get_shared_plane(drawn_stack_idx,c)
					.draw_line(f.get_move(fr_idx)(0)*scale,
										 f.get_move(fr_idx)(1)*scale,
										 f.get_move(fr_idx+1)(0)*scale,
										 f.get_move(fr_idx+1)(1)*scale,
										 color_map.data(color_index,1,1,c),1);
				}
			}
		}

	}
	template <class T>
	void  build_feature_tracks(const Feature<T> & f,  CImg<T> & current, int fr_idx){
		CImg<float> f_center=f.get_move(fr_idx);
		static T white[1]={current(f_center(0)*scale,f_center(1)*scale)};

		for(unsigned int i =0 ; i < min(f.get_move_list().size(),(unsigned int)11)-1; i++){
			draw_line(current,
								f.get_move(i)(0)*scale,
								f.get_move(i)(1)*scale,
								f.get_move(i+1)(0)*scale,
								f.get_move(i+1)(1)*scale,
								white,1);
		}

	}
	template <class T,class t>
	void build_feature_highlight(const Feature<T> & f,  CImg<t> & current){

		int color_index= cimg::rand()*255;
		cimg_for1(f.get_move_list().size(),fr_idx){
			CImg<float> f_center=f.get_move(fr_idx);
			cimg_forC(current,c){
				current.get_shared_plane(f.get_starting_frame_idx()+fr_idx,c)
					.draw_rectangle( f_center(0)*scale-window_size*scale/2-1,
													 f_center(1)*scale-window_size*scale/2-1,
													 f_center(0)*scale+window_size*scale/2+1,
													 f_center(1)*scale+window_size*scale/2+1,
													 color_map.data(color_index,1,1,c),1,~0U);
			}
		}
	}

	template <class T>
	void build_feature_highlight(const Feature<T> & f,  CImg<T> & current,unsigned int fr_idx){
		if(fr_idx >= f.get_move_list().size()) return ;
		CImg<float> f_center=f.get_move(fr_idx);
		static T white[1]={current(f_center(0)*scale,f_center(1)*scale)};

		current.draw_rectangle( f_center(0)*scale-window_size*scale/2-1,
														f_center(1)*scale-window_size*scale/2-1,
														f_center(0)*scale+window_size*scale/2+1,
														f_center(1)*scale+window_size*scale/2+1,
														white,1,~0U);
		// current.draw_rectangle( f_center(0)*scale-window_size*scale/2-2,
		// 												f_center(1)*scale-window_size*scale/2-2,
		// 												f_center(0)*scale+window_size*scale/2+2,
		// 												f_center(1)*scale+window_size*scale/2+2,
		// 												white,1,~0U);

	}

	template<class T>
	CImg<T>  build_frame(const CImg<T> & stack, typename std::list<Feature<T> > & feature_l ){
		CImg<T> feat_high(stack.get_slice(frame_idx));
		// feat_high.resize_tripleXY();
		typedef std::list<Feature<T> > t_feature_l;
		typename t_feature_l::iterator s;
		for (s=feature_l.begin(); s!=feature_l.end(); ++s){
			if(!s->get_lost()||(print_lost)) {
				build_feature_tracks<T>(*s,feat_high,frame_idx);
				build_feature_highlight<T>(*s,feat_high,frame_idx);
			}
		}
		// feat_high.crop(scale*261,scale*224,scale*(261+75),scale*(224+45));
		return feat_high;

	}


};
#endif /* _VAR_STAB_H_ */
