#ifndef _EVO_VISUALIZER_
#define _EVO_VISUALIZER_

//My include
#include "CImg.h"
#include "CFigure.h"
#include <vector>

using namespace cimg_library;
using namespace std;

typedef CImg<float> cfl;

struct Evo_visualizer{
	CImg<float> data;
	unsigned int data_idx;
	CFigure fig;
	Evo_visualizer(unsigned int max_data_nb){
		data.assign(max_data_nb);
		data_idx=0;
	}
	void assign(unsigned int max_data_nb){
		data.assign(max_data_nb);
		data_idx=0;
	}

	Evo_visualizer(){
		data.assign(1000);
		data_idx=0;
	}

	void append(Evo_visualizer & evo){
		data.assign(data._data,data_idx);
		data.append(evo.data,'x');
		data_idx+=evo.data_idx;
	}

	void reset(){data_idx=0;}

	template<class T>
	void  add_data(T d){
		data(data_idx++)=d;
	}

	cfl get_data(){
		cfl data_cut(data.data(),data_idx);
		return data_cut;
	}
	CFigure draw_hist(){
		cfl data_cut(data.data(),data_idx);
		CFigure fig;
		int bins=20;
		CImg<float> X(bins); X.sequence(data_cut.min(),data_cut.max());
		float blue[3]={0,0,255};

		fig.bar(X,data_cut.histogram(bins),blue);
		fig.erase();
		fig.set_axis();
		fig.replot();
		return fig;
	}

	Evo_visualizer& display_hist(){
		cfl data_cut(data.data(),data_idx);
		data_cut.histogram(data_idx/5).display_graph(0,3,1,0,data_cut.min(),data_cut.max()); 
		// int bins=100;
		// int bin_min=1;
		// int bin_max=4;
		// CImg<float> raw_hist;
		// raw_hist.assign(bins).get_sequence(bin_min,bin_max).get_transpose()
		// 	.append(lp.non_null.get_histogram(bins,bin_min,bin_max).get_transpose(),'x')
		// 	.save_dlm("lifetime_hist.dat");
		return *this;
	}

	Evo_visualizer& plot(CFigure & afig){
	// Use of an assymetric cost function to suppress the influence of the particle in the bg
		CImg<float> Y(data._data,data_idx);
		CImg<float> X(data_idx); X.sequence(0,data_idx-1);
		// afig.set_axis(X,Y);
		afig.plot(X,Y,"r-");
		return *this;
	};
	Evo_visualizer& plot(){
	// Use of an assymetric cost function to suppress the influence of the particle in the bg
		CImg<float> Y(data._data,data_idx);
		CImg<float> X(data_idx); X.sequence(0,data_idx-1);
		fig.clear();
		fig.set_axis(X,Y);
		fig.plot(X,Y,"r-");
		return *this;
	};
	void save(const char * filename){
		fig.save(filename);
	}
};

struct Mult_evo_visualizer{
	std::vector<Evo_visualizer> visu_list;
	CFigure fig;

	Mult_evo_visualizer(){}
	Mult_evo_visualizer(unsigned int param_nb){
		assign(param_nb);
	}

	void assign(unsigned int param_nb, int max_data_nb ){
		visu_list.clear();
		for (unsigned int i = 0; i < param_nb; ++i)
			{
				Evo_visualizer ev(max_data_nb);
				visu_list.push_back(ev);
			}
	}

	void assign(unsigned int param_nb){
		visu_list.clear();
		for (unsigned int i = 0; i < param_nb; ++i)
			{
				Evo_visualizer ev;
				visu_list.push_back(ev);
			}
	}


	void reset(){
		std::vector<Evo_visualizer>::iterator it;
		for (it = visu_list.begin();it != visu_list.end();it++){
			it->reset();
		};
	}

	CImgList<float>  get_data(){
		CImgList<> ret;
		int s=visu_list.size();
		cimg_for1(s,i) ret.push_back( visu_list[i].get_data());
		return ret;
	}
	
 template<class t>
	void  add_data(CImg<t> d){
	 int s=visu_list.size();
	 cimg_for1(s,i) visu_list[i].add_data(d[i]);
	}
	template<class T>
	void  add_data(T d0,T d1){
		 visu_list[0].add_data(d0);
		 visu_list[1].add_data(d1);
	}

	int nb_step(){
		return visu_list[0].data_idx;
	}
	cfl get_vector(unsigned int idx){
		cfl res(1,visu_list.size());
		cimg_foroff(res,i) res[i]=visu_list[i].data(idx);
		return res;
	}

	void append(Mult_evo_visualizer & mult_evo){
		for (unsigned int i =0  ; i< visu_list.size(); i++ ){
			visu_list[i].append(mult_evo.visu_list[i]);
		};
	};

	void plot(){
		std::vector<Evo_visualizer>::iterator it;
		for (it = visu_list.begin();it != visu_list.end();it++){
			it->plot();
		};
	};

	void save(const char * filename){
		unsigned int i=0;
		char graph_filename[400]; 
		for (i=0;i<visu_list.size();i++){
			cimg::number_filename(filename,i,2,graph_filename);
			visu_list[i].fig.save(graph_filename);
		};
	}
};

#endif /* _VAR_STAB_H_ */
