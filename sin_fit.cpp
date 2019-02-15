
#include "sin_fit.hpp"

Sin_fit::Sin_fit(){
	threshold = 0.01;
	max_iter= 100;
	alpha_line = 1;
	measure_nb=12;
	cost_function=NULL;
}
void Sin_fit::init_residue ( 	CImg<float> & aRESIDUE,
															CImg<float> & input_beta
															)
{
	residue=aRESIDUE.get_vector();
	beta=input_beta.get_vector();


	// beta_evo.assign(beta.size());
	spring.assign(1,beta.size()).fill(0);
	W.assign(residue.size(),residue.size());
	residue_weigth.assign(1,residue.size()).fill(1);
	J.assign(beta.size(),residue.size());
	delta_beta.assign(1,beta.size());
}		/* -----  end of method Sin_fit::Sin_fit  ----- */

void Sin_fit::init ( 	CImg<float> & input_X,
											CImg<float> & input_Y,
											CImg<float> & input_beta
											)
{
	// input parameter
	X=input_X.get_vector();
	Y=input_Y.get_vector();
	init_residue(Y,input_beta);
}		/* -----  end of method Sin_fit::Sin_fit  ----- */

// Init of sin_fit
void Sin_fit::init ( 	CImg<float> & input_X,
											CImg<float> & input_Y
											)
{
	beta.assign(3);
	init(input_X,input_Y,beta);
	set_starting_value();
}		/* -----  end of method Sin_fit::Sin_fit  ----- */

void Sin_fit::print_state ( 	ostream & of){
	of << "print fitter state" << endl;
	of << "Y" << endl;
	print_val(Y,of);
	of << "beta" << endl;
	print_val(beta,of );
	of << "J" << endl;
	print_val(J,of);
	of << "residue" << endl;
	print_val(residue,of );
	of << "delta_beta" << endl;
	print_val(delta_beta,of );
	of << "W" << endl;
	cimg_forX(W,x){ of << W(x,x) << "\t" ;}
	of << endl;
}

// Variance estimation used traditionnaly in covariance matrix
// estimation in the OLS case.  see Econometrix (Hayashi) for the
// proof
float Sin_fit::estimate_s_model (){
	float sigma = 0 ;
	float exponent = 2 ;
	fill_residue();
	cimg_forY(residue,i){
		sigma+=pow(residue(i),exponent);
	}
	sigma/=residue.size()-beta.size();
	return sigma;
}

CImg<float> Sin_fit::estimate_cov_matrix(){
	CImg<float> cov;
	float sig = estimate_s_model();
	fill_J();
	cov = sig*((J.get_transpose()*J).get_invert());
	return cov;
}

CImg<float> Sin_fit::cipra_robust_cov_matrix(){
	CImg<float> cov(beta.size(),beta.size());
	cov.fill(0);
	cimg_forY(residue,i){
		float r = residue(i);
		float w =cost_function->get_cost(r);
		cimg_forXY(cov,x,y){
			cov(x,y)=cov(x,y)+w*J(x,i)*J(y,i);
		}
	}
	cov.invert();
	cov*=cost_function->sig_noise; // MAD
	return cov;
}

CImg<float> Sin_fit::ieng_robust_cov_matrix(){
	CImg<float> O1(beta.size(),beta.size());
	CImg<float> O2(beta.size(),beta.size());
	O1.fill(0);
	O2.fill(0);
	float swrr=0;
	float sw=0;
	cimg_forY(residue,i){
		float r = residue(i);
		float w =cost_function->get_cost(r);
		swrr+=w*r*r;
		sw+=w;
		cimg_forXY(O1,x,y){
			float wjj=w*J(x,i)*J(y,i);
			O1(x,y)+=wjj;
			O2(x,y)+=w*wjj;
		}
	}
	CImg<float> O1_inv=O1.get_invert();
	CImg<float> O2O1_inv=O2*O1_inv;
	CImg<float>	cov=swrr*O1_inv*O2O1_inv/(sw-O2O1_inv.trace());
	return cov;
}


void Sin_fit::assign ( CImg<float> & input_Y
											 )
{
	CImg<float>  tmp_X(input_Y.size());
	cimg_forX(tmp_X,i){
		tmp_X(i)=(float)(i);
	}
	init(tmp_X,input_Y);
}		/* -----  end of method Sin_fit::Sin_fit  ----- */

Sin_fit::Sin_fit ( 	CImg<float> & input_X,
										CImg<float> & input_Y)
{
	init(input_X,input_Y);
}		/* -----  end of method Sin_fit::Sin_fit  ----- */

Sin_fit::Sin_fit ( 	CImg<float> & input_Y)
{
	CImg<float>  tmp_X(input_Y.size());
	cimg_forX(tmp_X,i){
		tmp_X(i)=(float)(i);
	}
	init(tmp_X,input_Y);
}		/* ----  end of method Sin_fit::Sin_fit  ----- */

void Sin_fit::assign ( 	CImg<float> & stack, int option)
{
	CImg<float>  tmp_X(stack.size());
	CImg<float>  tmp_Y(stack.size());
	unsigned int i = 0;
	cimg_forZ(stack,z){
		cimg_forXY(stack,x,y){
			tmp_Y(i)=stack(x,y,z);
			tmp_X(i)=(float)(z);
			i++;
		}
	}
	init(tmp_X,tmp_Y);
}		/* ----  end of method Sin_fit::Sin_fit  ----- */
float Sin_fit::reg_funct (float x)
{
	float Pi=3.141592;
	// float w=2.*Pi/(X.max()+1);
	float w=2.*Pi/measure_nb;
	return get_C()+get_amp()*sin(w*x+get_phase());
}		/* -----  end of method Sin_fit<float>::reg_funct  ----- */
float Sin_fit::cos_func (float x)
{
	return beta(0)+beta(1)*cos(beta(2)*x+beta(3));
}		/* -----  end of method Sin_fit<float>::reg_funct  ----- */


float Sin_fit::reg_funct_deriv ( int var_idx,float x)
{
	float res=0;
	float pi=3.141592;
	float w=2.*pi/measure_nb;

	switch ( var_idx ) {
	case 0:
		res = 1;
		break;

	case 1:
		res = sin(w*x+beta(2));
		break;

		// case 2:
		// 	res = beta(1)*x*cos(beta(2)*x+beta(3));
		// 	break;

	case 2:
		res = beta(1)*cos(w*x+beta(2));
		break;

	default:
		printf("error : not enough variable\n");
		break;
	}				/* -----  end switch  ----- */
	return res;
}		/* -----  end of method sin_fit<float>::reg_funct_de & riv_0  ----- */
float Sin_fit::cos_func_deriv ( int var_idx,float x)
{
	float res=0;

	switch ( var_idx ) {
	case 0:
		res = 1;
		break;

	case 1:
		res = cos(beta(2)*x+beta(3));
		break;

	case 2:
		res = -beta(1)*x*sin(beta(2)*x+beta(3));
		break;

	case 3:
		res = -beta(1)*sin(beta(2)*x+beta(3));
		break;

	default:
		printf("ERROR : not enough variable\n");
		break;
	}				/* -----  end switch  ----- */
	return res;
}		/* -----  end of method Sin_fit<float>::reg_funct_de & riv_0  ----- */



void Sin_fit::fill_J ()
{
	cimg_forXY(J,j,i){
		J(j,i)=residue_weigth(i)*reg_funct_deriv(j,X(i));
	}
}
void Sin_fit::recenter_beta ()
{
	// static float C=Y.mean();
	// if (abs(beta(0)-C)> beta(1) )				beta(0)=C;


	// float Pi=3.141592;
	// float w=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	// float amp=(Y.max()-Y.min())/2;	/* FIXME */
	// float phase=atan(2.5*w)+2.42342-atan(4.1*w);	/* FIXME */
	// if (abs(beta(1)-amp)> amp )		beta(1)=amp;
	// if (abs(beta(2)-phase)> 2 )		beta(2)=phase;
}



void Sin_fit::fill_spring ()
{
	// static float Pi=3.141592;
	// static float w=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	// float phase=atan(2.5*w)+2.42342-atan(4.1*w);	/* FIXME */
	float delta=1;
	static float C=Y.mean();
	spring(0)=delta*(beta(0)-C);
}

void Sin_fit::fill_residue ()
{

	// fill_spring();
	// recenter_beta();
	cimg_forY(residue,i){
		residue(i)=residue_weigth(i)*(Y(i)-reg_funct(X(i)));
	}

}		/* -----  end of method Sin_fit::build_delta_beta  ----- */



void Sin_fit::set_starting_value (  )
{
	float Pi=3.141592;
	float w=2*Pi*0.04;	// the Li-Flim software is usually set up to 40 Mhz
	set_C(Y.mean());
	set_amp((Y.max()-Y.min())/2);	/* FIXME */
	// beta(2)=2.*Pi/12.;                                   /* FIXME */
	set_phase(atan(2.5*w)+2.42342-atan(4.1*w));	/* FIXME */
}		/* -----  end of method Sin_fit::starting_value  ----- */

void Sin_fit::run ( CImg<float> & beta_output )
{
	CImg<float> J_t, J_prod(beta.size(),beta.size()), E(1,beta.size()), L, L1;
	nb_iter = 0;
	// Warning ! the order of execution of fill_residue/fill_J is very important
	// due to possible problem with overriding of those method by subclass ( they
	// can share some data).
	CImg<float> beta_init_value = beta;
	float test_value;
	// beta_evo.reset();
	// beta_evo.add_data(beta);
	float r_old=0,r_new=1;
	do {
		r_old=r_new;
		fill_residue();
		r_new=residue.magnitude();
		test_value=std::abs(r_old-r_new)/r_old;

		fill_J();
		E.fill(0);
		J_prod.fill(0);

		cimg_forY(residue,i){
			cimg_forY(delta_beta,j){
				E(j)=E(j)+residue(i)*J(j,i);
			}
			cimg_forXY(J_prod,x,y){
				J_prod(x,y)=J_prod(x,y)+J(x,i)*J(y,i);
			}
		}
		delta_beta=E.solve(J_prod);

		beta+=delta_beta;
		nb_iter++;
		// beta_evo.add_data(beta);
	}while	( (test_value  >threshold ) && (nb_iter < max_iter) );
	// }while	( ( ((delta_beta.get_div(beta)).get_abs().max())>threshold ) && (nb_iter < max_iter) );

	if (beta(0)!=beta(0)){
		// print_state(cout);
		beta=beta_init_value;
		nb_iter=max_iter;
	}

	beta_output=beta;
}		/* -----  End of method Sin_fit::run  ----- */


void Sin_fit::update_weights (  CImg<float> err_sig, float alpha, float (*pt2CostFunc)(float, float,float),
																CImg<float> & W){
	cimg_forX(W,x){
		W(x,x)=(*pt2CostFunc)(residue(x),err_sig(x),alpha)/residue(x);
	}
}

void Sin_fit::update_weights (  float err_sig, float alpha, float (*pt2CostFunc)(float, float,float),
																CImg<float> & W){
	cimg_forX(W,x){
		W(x,x)=(*pt2CostFunc)(residue(x),err_sig,alpha)/(residue(x)+0.0001);
	}
}

// IRLS stands for iteratively reweighted least square
// this method  suppose the beta parameter and the residue are already well initialized
// - beta value is set
// - scale parameter is set
void Sin_fit::irls_no_init  ( CImg<float> & beta_output )
	{
		CImg<float > L;
		CImg<float > L1;

		beta_output=beta;

		nb_iter = 0 ;
		float w= 0;
		CImg<float > J_prod(beta.size(),beta.size());
		CImg<float > E(1,beta.size());

		float r_old=0,r_new=1,test_value;
		do {
			r_old=r_new;
			fill_residue();
			r_new=residue.magnitude();
			test_value=std::abs(r_old-r_new)/r_old;

			fill_J();
			// update_weights(err_sig,alpha,pt2CostFunc,W);
			E.fill(0);
			J_prod.fill(0);
			cimg_forY(residue,i){
				// w =(*pt2CostFunc)(residue(i),err_sig,alpha)/(residue(i));
				w =cost_function->get_cost(residue(i));
				cimg_forY(delta_beta,j){
					E(j)=E(j)+residue(i)*w*J(j,i);
				}

				cimg_forXY(J_prod,x,y){
					J_prod(x,y)=J_prod(x,y)+w*J(x,i)*J(y,i);
				}
			}
			cholesky_decomp(J_prod,L,L1);
			delta_beta=E.solve(L);
			delta_beta=delta_beta.solve(L1);


			// delta_beta=((J_t*W)*residue).solve(J_prod);
			beta+=delta_beta;
			nb_iter++;


			// beta_evo.add_data(beta);
		} while ( ( test_value>threshold ) && (nb_iter < max_iter) );
		// } while ( ( ((delta_beta.get_div(beta)).get_abs().max())>threshold ) && (nb_iter < max_iter) );
		// if (nb_iter>= max_iter)			{cout << "MAX ITER reach during fitting" <<endl;}
		// else												{cout << " M-estimatin converge normally " << endl;}
		if (beta(0)!=beta(0)) {
			// print_state(cout);
			// cout << " M-estimation diverge" <<endl;
			beta=beta_output;
			nb_iter=max_iter;
		}

		beta_output=beta;
	}		/* -----  End of method Sin_fit::run  ----- */

	// IRLS stands for iteratively reweighted least square
	// this version handle the initialisation of parameters using a standard least square
void Sin_fit::irls ( CImg<float> & beta_output  ,
											 float err_sig,
											 float alpha,
											 float (*pt2CostFunc)(float, float,float),
											 CImg<float> * err_sigs=NULL)
	{
		CImg<float> J_t;
		CImg<float > J_prod;
		CImg<float > L;
		CImg<float > L1;

		err_sig=scale_parameter;
		alpha_line=1;							// idem
		W.identity_matrix();
		run(beta_output);
		beta=beta_output;

		// threshold = 0.05;
		// max_iter = 20;
		nb_iter = 0 ;

		do {
			fill_residue();

			if( err_sigs == NULL) {
				update_weights(err_sig,alpha,pt2CostFunc,W);
			}
			else{
				update_weights(*err_sigs,alpha,pt2CostFunc,W);
			}
			fill_J();
			J_t=J.get_transpose();
			J_prod=J_t*W*J;
			cholesky_decomp(J_prod,L,L1);
			delta_beta=((J_t*W)*residue+spring).solve(L);
			delta_beta=delta_beta.solve(L1);
			// delta_beta=((J_t*W)*residue).solve(J_prod);
			beta+=alpha_line*delta_beta;
			nb_iter++;



		} while ( ( delta_beta.magnitude()>threshold ) && (nb_iter < max_iter) );
		// } while ( ( ((delta_beta.get_div(beta)).get_abs().max())>threshold ) && (nb_iter < max_iter) );
		if (nb_iter>= max_iter) cout << "MAX ITER reach during fitting" <<endl;
		if (beta(0)!=beta(0)) {
			print_state(cout);
			beta=beta_output;
			nb_iter=max_iter;
		}

		beta_output=beta;
	}		/* -----  End of method Sin_fit::run  ----- */

	// void Sin_fit::w_fit ( CImg<float> & beta_output , CImg<float> & W){
	// 	CImg<float> J_t;
	// 	CImg<float > J_prod;
	// 	CImg<float > L;
	// 	CImg<float > L1;

	// 	set_starting_value();
	// 	int i=0;
	// 	int max_iter=200;							// arbitrary number.
	// 	float alpha=0.01;
	// 	do {
	// 		fill_residue();
	// 		fill_J();
	// 		J_t=J.get_transpose();
	// 		J_prod=(J_t*W)*J;

	// 		cholesky_decomp(J_prod,L,L1);
	// 		delta_beta=((J_t*W)*residue).solve(L);
	// 		delta_beta=delta_beta.solve(L1);
	// 		// 		delta_beta=((J_t*W)*residue).solve(J_prod);
	// 		beta+=alpha*delta_beta;
	// 		i++;


	// 	} while ( ( ((delta_beta.get_div(beta)).get_abs().max())>threshold ) && (i < max_iter) );

	// 	beta_output=beta;

	// }

	static	void
	fill_id ( CImg<complex<float> > & M )
	{
		cimg_forXY(M,x,y)	{
			M(x,y)=0;
		}
		cimg_forX(M,x)	{
			M(x,x)=1;
		}
	}		/* -----  end of function fill_id  ----- */

	void Sin_fit::cholesky_decomp ( 	CImg<float>  & input,
																		CImg<float>  & L,
																		CImg<float>  & L_trans_conj)
	{

		L=input;
		L_trans_conj=input;
		CImg<complex<float> > A(input);
		CImg<complex<float> > A_copy(input);
		CImg<complex<float> > id(input);
		fill_id(id);
		CImg<complex<float> > L_tmp(id);
		CImg<complex<float> > L_res(id);

		cholesky_decomp_rec ( 	A, id, A_copy, L_tmp, 0, L_res);

		cimg_forXY(L,x,y){
			L(x,y)=real(L_res(x,y));
			L_trans_conj(x,y)=L(y,x);
		}

		L_trans_conj=L.get_transpose();
	}		/* -----  end of method Sin_fit::cholesky_decomp  ----- */


	void Sin_fit::cholesky_decomp_rec ( 	CImg<complex<float> > &A,
																				CImg<complex<float> > &id,
																				CImg<complex<float> > & A_copy,
																				CImg<complex<float> > & L,
																				int n,
																				CImg<complex<float> > & res)
	{
		if(A!=id){
			L=id;
			L(n,n)=sqrt(A(n,n));
			for (int i=n+1 ; i<A.height() ; i++){
				L(n,i)=A(n,i)/L(n,n);
			}
			A_copy.assign(A);
			A=id;
			for (int j=n+1 ; j<A.width() ; j++){
				for (int i=n+1 ; i<A.height() ; i++){
					A(j,i)=A_copy(j,i)-(A_copy(n,i)*conj(A_copy(j,n)))/A_copy(n,n);
				}
			}
			cholesky_decomp_rec ( A,id,A_copy,L,n+1,res*=L);

		}

	}		/* -----  end of method Sin_fit::cholesky_decomp  ----- */


// static const char red_bar[10]="r-";
void Sin_fit::plot_weighted_fit(CFigure & fig ,const char * color = "r-"){
	CImg<float> Y_reg(X.size());			//
	cimg_foroff(Y_reg,i){
		Y_reg[i]=residue_weigth(i)*reg_funct(X[i]);
	}
	fig.plot(X,Y_reg,color);
	fig.set_axis();
	fig.erase();
	fig.replot();
}

void Sin_fit::plot_fit(CFigure & fig ,const char * color = "r-"){
	int pixel_curb_res = 2;
	CImg<float> X_reg(fig._x(X.max())/pixel_curb_res);			//
	CImg<float> Y_reg(fig._x(X.max())/pixel_curb_res);			//
	cimg_foroff(X_reg,i){
		X_reg[i]=fig._ix(i*pixel_curb_res);
		Y_reg[i]=reg_funct(X_reg[i]);
	}
	fig.plot(X_reg,Y_reg,color);
	// fig.set_axis();
	// fig.erase();
	// fig.replot();
}

CFigure Sin_fit::build_1d_weighted_graph(){
		CFigure fig;
		build_1d_weighted_graph(fig);
		return fig;
}
void Sin_fit::build_1d_weighted_graph(CFigure &fig){
		CImg<float> weighted_Y=Y.get_mul(residue_weigth);
		fig.set_axis(X,weighted_Y);
		fig.plot(X,weighted_Y,"bx");
		plot_weighted_fit(fig);
}
CFigure Sin_fit::build_1d_graph(const char * marker){
	CFigure fig;
		fig.set_axis(X,Y);
		fig.plot(X,Y,marker);
		plot_fit(fig);
		return fig;
}
void Sin_fit::build_1d_graph(CFigure & fig, const char * marker){
		fig.plot(X,Y,marker);
		plot_fit(fig,"r-");
}

void  Sin_fit::save_X_Y_gnuplot(const char * filename){
	X.get_append(Y,'x').save_dlm(filename);
}

CFigure Sin_fit::plot_residue(CFigure & fig){
		fig.set_axis(X,residue);
		fig.plot(X,residue,"bx");
		return fig;
}

CFigure Sin_fit::plot_residue(){
		CFigure fig;
		fig.set_axis(X,residue);
		fig.plot(X,residue,"bx");
		return fig;
}

CImg<unsigned char> Sin_fit::display_fit(){
	CFigure fig = build_1d_graph();
	fig.display();
	return fig.get_image();
}
// CImg<unsigned char> Sin_fit::display_fit()
// 	{
// 		CImg<unsigned char> graph(600,500,1,3,0);
// 		CImg<float> value(110);
// 		float i;
// 		for (i = 0; i < 110.; i++) {
// 			value(i)=reg_funct(i/10.);
// 		}
// 		const unsigned char white[]={255,255,255};
// 		const unsigned char red[]={255,0,0};
// 		// CImgDisplay  draw_disp(graph,"Intensity profile");

// 		// while (!draw_disp.is_closed()) {
// 		// 	draw_disp.wait();
// 		// 	graph.draw_graph(value,white,1,0,2,0,0,40000,true);
// 		// 	// graph.draw_graph(Y.get_transpose(),red,1,0,6,0,0,40000,true);
// 		// }	//
// 		graph.draw_graph(value,white,1,0,2,60000,0,true);
// 		graph.draw_graph(Y.get_transpose(),red,1,0,6,60000,0,true);
// 		return graph;
// 	}

// if beta.size()== 2
// save data to be read by gnuplot displaying a plan aspect of the residuals
void Sin_fit::generates_gnuplot_residual_plan(const char * file_name){
	cfl matrix=J.get_append(residue,'x');
	matrix.save_dlm(file_name);
}

void Sin_fit::test (  )
{
	CImg<float> beta_truth(1,4,1,1,5.,3.,0.5,1.);
	CImg<float> beta_output;
	X.assign(1,12);
	Y.assign(1,12);
	beta=beta_truth;
	cout << "test sample :" ;
	cimg_forY(X,k){
		X(k)=(float)k;
		Y(k)=reg_funct(k);
		cout << Y(k) << " ";
	}
	cout << endl;

	beta.fill(0);
	run(beta_output);
	beta_truth.print("truth",false);
	beta_output.print("res",false);
}		/* -----  end of method Sin_fit::test  ----- */
cfl empirical_cov_matrix(cfl & betas,cfl & mean){
	cfl betas_tmp(betas);
	mean.assign(betas.height());
	cfl line;
	cimg_forY(betas_tmp,i){
		line.assign(betas_tmp.get_shared_line(i),true);
		mean(i)=line.mean();
		line-=mean(i);
	}
	cfl cov_matrix=betas_tmp*betas_tmp.get_transpose()/(betas_tmp.width()-1);
	return cov_matrix;
}
