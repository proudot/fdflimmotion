/**
** \file   noise_analysis.h
** \author Jérôme Boulanger (jerome.boulanger@laposte.net)
** \date   2008-2009
** 
** CImg plugin for variance noise stabilization using Generalized
** Anscombe Transform. 
**
** Started on  Thu Feb 12 10:59:46 2009 Jérôme Boulanger
**
** Last update Tue May  5 10:04:01 2009 Jérôme Boulanger
**/

#pragma once
#ifndef cimg_plugin_noise_analysis
#define cimg_plugin_noise_analysis

// to have isnan
#ifndef isnan
#define isnan(x) ((x==x?false:true))
#endif


//! Compute pseudo-residuals in 2D
CImg<Tfloat> get_pseudo_residuals() const {
  CImg<Tfloat> dest(*this);
  CImg_3x3(I,float);
  const Tfloat c=(Tfloat)(1.0/std::sqrt(20.0));
  cimg_forZC(*this,z,v) cimg_for3x3(*this,x,y,z,v,I,Tfloat) 
    dest(x,y,z,v) = c*(Tfloat)((T)4.0*Icc-(Icn+Icp+Ipc+Inc));
  return dest;
}

//! Estimate the noise variance
/**
  \param variance_method = 0 : Least Median of Square, 
                           1 : Least Trimmed of Square, 
	                   2 : Least Mean of Square.
   Robustly estimatate the variance of a the noise using the pseudo-residuals.
   \see variance()
**/
double noise_variance(const int variance_method=2) const{
  return (*this).get_pseudo_residuals().variance(variance_method);
}


//! Robust mean estimation
double robust_mean_estimation(int mean_method=2) const{
  if (is_empty())
    throw CImgInstanceException("CImg<%s>::robust_mean_estimation() : Instance image (%u,%u,%u,%u,%p) is empty.",
				pixel_type(),width(),height(),depth(),spectrum(),data());
  if (mean_method==1) return median();  
  double m=0,s=0,n=0;
  cimg_foroff((*this),i) {
    const double val = (*this)(i);
    m+=val;
    s+=val*val;
    n++;
  }
  m/=n;s/=n;s-=m*m;s=1./(3.*std::sqrt(s));
  double m1=m,m2=m,k=0;
  if (mean_method>0){
    do{
      m2 = m1;
      double sw=0;
      m1=0;
      cimg_foroff((*this),i){
	const double val = (*this)(i), x = (val-m2)*s;
	double w;
	switch(mean_method){
	case 2: w = std::exp(-.5*x*x); break; // Leclerc/Welsh
	case 3: w = x<1.0?1:1./x; break; // huber
	case 4: w = 1/((1+x*x)*(1+x*x)); break; // Geman McClure
	case 5: w = 1/(1+x*x); break; // cauchy
	default: throw CImgException("Undefined M-estimator in robust_mean_estimation()");
	}
	m1 += w*val;
	sw += w;
      }
      m1/=sw;
      k++;
    }while(fabs(m2-m1)/(m2+m1)<1.0e-6 && k<20);
  }
  return m1;
}


//! Print statistics of the mean and noise variance estimated.
template<typename t>
void gat_print_stat(const CImg<t> & x,const CImg<t> & y) const {
  CImg<Tfloat> sX=x.get_stats(),sY=y.get_stats();
  printf("   Stats : E[Y]={%.3f, %.3f}, V[Y]={%.3f, %.3f}, avg=%.3f \n", 
	 sX(0),sX(1), sY(0), sY(1),sY(3));
}

// robust linear fitting y = a+ b x
void linear_fit(const CImg<> &x, const CImg<> & y, 
								double & a, double & b, const int robust=0, CImg<double> & cov = CImg<double>::empty()) const {
  if (x.is_empty()) throw CImgException("linear_fit() in file %s at line %d: x is empty",__FILE__,__LINE__);
  if (y.is_empty()) throw CImgException("linear_fit() in file %s at line %d: y is empty",__FILE__,__LINE__);
  if (x.size()!=y.size()) throw CImgException("linear_fit() in file %s at line %d: CImg<>x(%d,%d,%d,%d) and"
																							"CImg<>y(%d,%d,%d,%d) have not the same size",
																							x.width(),x.height(),x.depth(),x.depth(),
																							y.width(),y.height(),y.depth(),y.depth(),
																							__FILE__,__LINE__);
  double sw=0,swxx=0,swx=0,swy=0,swxy=0, sww = 0 , swrr=0;
  // init least square
  const float *ptr_xi=x.data(), *ptr_yi=y.data();
  for (unsigned int i=0;i<x.size();++i) {
    const float xi = *ptr_xi++, yi = *ptr_yi++;    
    sw++;
    swx += xi;     swy += yi;
    swxx += xi*xi; swxy += xi*yi;
  }
  double d = sw*swxx-swx*swx;
  if (cimg::abs(d)>1.0e-10) {
    a = (swxx*swy-swx*swxy)/d;
    b = (sw*swxy-swx*swy)/d;
  }else{
    a=0; b=1;
#if cimg_debug>=3
    cimg::warn(" %s line %d in linear_fit() : Dividing by %g.\n",__FILE__,__LINE__,d);
#endif   
  }
#if cimg_debug>=3
  printf("   Fitting OLS :  Var[Y]= %.3f E[Y] +  %.3f \n",b,a);
#endif 
  if (robust==1) {
    // compute variance of residuals
		cov.assign(2,2);
    double sigma=0,mean=0,nsamples=0,r_new=0,r_old=0;
    ptr_xi=x.data(); ptr_yi=y.data();
    for (unsigned int i=0;i<x.size();++i) {
      const double xi = *ptr_xi++, yi=*ptr_yi++, delta = yi-(a+b*xi);
      mean += delta; sigma += delta*delta; nsamples++;
    }
    mean /= nsamples; sigma /= nsamples; r_new=sigma; sigma -= mean*mean; 
#if cimg_debug>=3
    printf("   Residuals : %.3g, sigma : %.3g, mean : %.3g\n",r_new,std::sqrt(sigma),mean);
#endif
    // IRLS
    unsigned int k=0;
    while ( k<1000 && (k<5 || cimg::abs((r_old-r_new)/r_old)>1e-12) )  {
      ptr_xi=x.data(); ptr_yi=y.data();
      r_old=r_new; r_new=0; sw=0; swxx=0; swx=0; swy=0; swxy=0;sww = 0 ; swrr=0;
      for (unsigned int i=0;i<x.size();++i){    
				const double xi = *ptr_xi++, yi=*ptr_yi++, 
					delta = yi-(a+b*xi), w = std::exp(-0.5*delta*delta/(9.0*sigma));
				if (!isnan(w)){
					sw += w;
					sww+=w*w;
					swrr+=w*delta*delta;
					swx += w*xi;     swy += w*yi;
					swxx += w*xi*xi; swxy += w*xi*yi;
					r_new += delta*delta;

				}
      }
			/* implementation of p.99 J. Boulanger thesis */
			cov(0,0)=sw; cov(0,1)=swx; cov(1,0)=swx; cov(1,1)=swxx;
			cov.invert();
			cov/=sw*sw;
			cov*=swrr*sww;

			/* Implementation of classic variance estimation */
			/* cov.fill(r_new/(y.size()-2)); */

      d = sw*swxx-swx*swx; r_new /= nsamples;
      if (cimg::abs(d)>1/(cimg::type<double>::max()-1)) {
				a = (swxx*swy-swx*swxy)/d;
				b = (sw*swxy-swx*swy)/d;
      }else
				cimg::warn("Warning in %s line %d in linear_fit() : Dividing by %g.",__FILE__,__LINE__,d);
      ++k;
    }
#if cimg_debug>=3
    printf("   Fitting IRLS(%d) : Var[Y]= %.3f E[Y] +  %.3f , R=%.3g\n",k,b,a,r_new);
#endif    
  }
}

//! Save a plot of the regression using gnuplot
/**
   \note gnuplot should be in the path.
 **/
void gat_save_graph(const CImg<> & x, const CImg<> & y, 
		    const double g0, const double edc, 
		    const char * file_o) const {
  if (!file_o) throw CImgException("gat_save_graph() in file %d at line %d: file is null",__FILE__,__LINE__);
  if (x.is_empty()) throw CImgException("gat_save_graph() in file %d at line %d: x is empty",__FILE__,__LINE__);
  if (y.is_empty()) throw CImgException("gat_save_graph() in file %d at line %d: y is empty",__FILE__,__LINE__);
  if (x.size()!=y.size()) throw CImgException("gat_save_graph() in file %d at line %d: CImg<>x(%d,%d,%d,%d) and"
					      "CImg<>y(%d,%d,%d,%d) have not the same size",__FILE__,__LINE__,
					      x.width(),x.height(),x.depth(),x.depth(),
					      y.width(),y.height(),y.depth(),y.depth());
  char buf[1024];
  const char * file_path="./";
  sprintf(buf,"%s/noise_analysis.txt",file_path);
  // Save results
  FILE * f = fopen(buf,"w");
  for (unsigned int k=0;k<x.size();++k) fprintf(f,"%f\t%f\n",x(k),y(k));  
  fclose(f);
  // Display results
  sprintf(buf,"%s/script_gnuplot.plt",file_path);
  f = fopen(buf,"w");
  const char *ext = cimg::split_filename(file_o);
  if (!cimg::strcasecmp(ext,"tex")) fprintf(f,"set terminal epslatex 'default' 8\n");
  if (!cimg::strcasecmp(ext,"svg")) fprintf(f,"set terminal svg 'default' 8\n");
  if (!cimg::strcasecmp(ext,"eps")) fprintf(f,"set terminal ps 'default' 8\n");
  if (!cimg::strcasecmp(ext,"pdf")) fprintf(f,"set terminal pdf \n");
  if (!cimg::strcasecmp(ext,"png")) fprintf(f,"set terminal png\n");   
  /* if (!cimg::strcasecmp(ext,"png")) fprintf(f,"set terminal png 'default' 8\n"); */
  if (!cimg::strcasecmp(ext,"emf")) fprintf(f,"set terminal emf 'Time Roman' 8\n");
  fprintf(f,"set output '%s' \n",file_o);
  fprintf(f,"set size .6, .7\n");
  fprintf(f,"set xlabel '$\\mathbb{E}[Y]$'\n");
  fprintf(f,"set ylabel '$Var[Y]$'\n");  
  /* fprintf(f,"set ylabel '$Var[Y]$' .5,.0 \n");   */
  fprintf(f,"set grid\n");
  fprintf(f,"set xrange [%d to %d]\n",(int)x.min(),(int)x.max());
  /* fprintf(f,"set yrange [%d to %d]\n",(int)y.min(),(int)y.max()/50); */
  fprintf(f,"plot '%s/noise_analysis.txt' with dots notitle, %f*x+%f notitle\n",
	  file_path,g0,edc);
  fclose(f);
  sprintf(buf,"gnuplot %sscript_gnuplot.plt",file_path);
  if (system(buf)){}  
}

// Remove values for which the variance is negative
void gat_cleanup(CImg<> & x, CImg<> & y,CImg<> & size) const {
  if (x.size()!=y.size()) 
    throw CImgException("gat_cleanup() in file %d at line %d: CImg<>x(%d,%d,%d,%d) and CImg<>y(%d,%d,%d,%d)"
			" have not the same size",__FILE__,__LINE__,x.width(),x.height(),x.depth(),x.depth(),
			y.width(),y.height(),y.depth(),y.depth());
  // count the number of positive values
  unsigned int n0=0;
  cimg_foroff(x,off) if(y(off)>=0 && !isnan(y(off)) && !isnan(x(off))) n0++;
  if (n0<x.size()){
    // fill temporary vectors with positive values
    CImg<double> tx(n0),ty(n0),ts(n0);
    int i=0;
    cimg_foroff(x,off) if(y(off)>=0 && !isnan(y(off)) && !isnan(x(off))) {tx(i)=x(off);ty(i)=y(off);ts(off)=i++;}
    // replace the old by the new one
    x=tx;
    y=ty;
		size=ts;
  }
}

void gat_cleanup(CImg<> & x, CImg<> & y) const {
  if (x.size()!=y.size()) 
    throw CImgException("gat_cleanup() in file %d at line %d: CImg<>x(%d,%d,%d,%d) and CImg<>y(%d,%d,%d,%d)"
			" have not the same size",__FILE__,__LINE__,x.width(),x.height(),x.depth(),x.depth(),
			y.width(),y.height(),y.depth(),y.depth());
  // count the number of positive values
  unsigned int n0=0;
  cimg_foroff(x,off) if(y(off)>=0 && !isnan(y(off)) && !isnan(x(off))) n0++;
  if (n0<x.size()){
    // fill temporary vectors with positive values
    CImg<double> tx(n0),ty(n0);
    int i=0;
    cimg_foroff(x,off) if(y(off)>=0 && !isnan(y(off)) && !isnan(x(off))) {tx(i)=x(off);ty(i)=y(off);i++;}
    // replace the old by the new one
    x=tx;
    y=ty;
  }
}

void gat_compute_levellines_stats(const CImg<> & img, CImg<> &mean, CImg<> &var,const int mean_method=1,const int variance_method=2)const {
  CImg<> imgb=img.get_blur(2.f);
  float mini,maxi=imgb.max_min(mini);
  mini=floor(mini); maxi=ceil(maxi);
  const int nlevels=100;
#if cimg_debug>=3
  printf("   Samples : %d in [%d,%d]\n",nlevels,(int)mini,(int)maxi);
#endif
  mean.assign(nlevels);
  var.assign(nlevels);
  for (int i=0;i<nlevels;++i) {
    int n=0; // number of values for level i
    CImg<> tmp(img.size()); // intensity at level i
    float level=mini+(float)i*(maxi-mini)/(float)nlevels,
      leveln=mini+(float)(i+1)*(maxi-mini)/(float)nlevels;
    cimg_foroff(imgb,j) {
      if (imgb(j)>level && imgb(j)<leveln) {
	tmp(n)=img(j);
	n++;
      }
    }
    if (n>10) {
      tmp.crop(0,n-1);
      mean(i)=tmp.robust_mean_estimation(mean_method);
      var(i)=tmp.variance(variance_method);
    }
  }
  gat_cleanup(mean,var);
}

void gat_quadtree_divide(const int x0, const int y0, const int x1, const int y1,
												 const CImg<> &img, const CImg<> & residuals,
												 CImg<> & mean, CImg<> & variance, CImg<> & size,
												 int & n, CImgDisplay & disp, CImg<> &visu,const int mean_method=1, 
												 const int variance_method=2,const float threshold=1.3) const { 
  const int bs=1; // minimal size of the blocks
  float mmm=0;  if (!disp.is_closed()) mmm=img.max();
  // compute local statistics
  double vi=0,mi=0,ni=0,vn=0,mn=0; 
  for (int y=y0;y<=y1;++y) for (int x=x0;x<=x1;++x) {
      const double vali=img(x,y), valn=residuals(x,y);
      mi+=vali; vi+=vali*vali;
      mn+=valn; vn+=valn*valn;
      ni++;
    }
  vn=(vn-mn*mn/ni)/ni;					/* residual variance */
	vi=(vi-mi*mi/ni)/ni;					/* image variance */
  mi/=ni;												/* image mean */
  float col1[3]={mi,vn,0},col2[3]={mmm,mmm/2,mmm/4};
  // divide the block if the the variance of the image > variance of the noise

  if ( (x1-x0)>2*bs && (y1-y0)>2*bs && cimg::max(vn,vi)/cimg::min(vn,vi)>threshold){
    const int xc=(int)((float)(x1+x0+1)/2.),yc=(int)((float)(y1+y0+1)/2.);
    gat_quadtree_divide(x0,y0,xc,yc,img,residuals,mean,variance,size,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(xc,y0,x1,yc,img,residuals,mean,variance,size,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(x0,yc,xc,y1,img,residuals,mean,variance,size,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(xc,yc,x1,y1,img,residuals,mean,variance,size,n,disp,visu,mean_method,variance_method,threshold);

#if cimg_display==1 || cimg_display==2
    if (!disp.is_closed()) visu.draw_rectangle(x0,y0,x1,y1,col2,.5,~0L).display(disp);
#endif
  }
	else
		if (x1-x0>3*bs && y1-y0>3*bs)
	{
    mean(n)=img.get_crop(x0,y0,x1,y1).robust_mean_estimation(mean_method);
		variance(n)=residuals.get_crop(x0,y0,x1,y1).variance(variance_method);
		size(n)=((x1-x0+1)*(y1-y0+1));
		n++;

#if cimg_display==1 || cimg_display==2
      if (!disp.is_closed())
	visu.draw_rectangle(x0,y0,x1,y1,col1,.5)
	  .draw_rectangle(x0,y0,x1,y1,col2,.7,~0L)
	  .display(disp);
#endif
  }
}

void gat_quadtree_divide(const int x0, const int y0, const int x1, const int y1,
			 const CImg<> &img, const CImg<> & residuals,
			 CImg<> & mean, CImg<> & variance, int & n,
			 CImgDisplay & disp, CImg<> &visu,const int mean_method=1, 
			 const int variance_method=2,const float threshold=1.3) const { 
  const int bs=5; // minimal size of the block
  float mmm=0;  if (!disp.is_closed()) mmm=img.max();
  // compute local statistics
  double vi=0,mi=0,ni=0,vn=0,mn=0; 
  for (int y=y0;y<=y1;++y) for (int x=x0;x<=x1;++x) {
      const double vali=img(x,y), valn=residuals(x,y);
      mi+=vali; vi+=vali*vali;
      mn+=valn; vn+=valn*valn;
      ni++;
    }
  vn=(vn-mn*mn/ni)/ni;  vi=(vi-mi*mi/ni)/ni;  mi/=ni;
  float col1[3]={mi,vn,0},col2[3]={mmm,mmm/2,mmm/4};
  // divide the block if the the variance of the image > variance of the noise
  if (x1-x0>bs && y1-y0>bs && cimg::max(vn,vi)/cimg::min(vn,vi)>threshold){
    const int xc=(int)((float)(x1+x0+1)/2.),yc=(int)((float)(y1+y0+1)/2.);
    gat_quadtree_divide(x0,y0,xc,yc,img,residuals,mean,variance,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(xc,y0,x1,yc,img,residuals,mean,variance,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(x0,yc,xc,y1,img,residuals,mean,variance,n,disp,visu,mean_method,variance_method,threshold);
    gat_quadtree_divide(xc,yc,x1,y1,img,residuals,mean,variance,n,disp,visu,mean_method,variance_method,threshold);
#if cimg_display==1 || cimg_display==2
    if (!disp.is_closed()) visu.draw_rectangle(x0,y0,x1,y1,col2,.5,~0L).display(disp);
#endif
  }else if (x1-x0>bs && y1-y0>bs){
    mean(n)=img.get_crop(x0,y0,x1,y1).robust_mean_estimation(mean_method);
      variance(n)=residuals.get_crop(x0,y0,x1,y1).variance(variance_method);
#if cimg_display==1 || cimg_display==2
      if (!disp.is_closed())
	visu.draw_rectangle(x0,y0,x1,y1,col1,.5)
	  .draw_rectangle(x0,y0,x1,y1,col2,.7,~0L)
	  .display(disp);
#endif
      n++;
  }
}

void gat_compute_quadtree_stats(const CImg<> & img, CImg<> &mean, CImg<> &var, CImg<> & size, const int mean_method=1,
				const int variance_method=2,const bool visualize=false) const {
  CImg<> visu;
  CImgDisplay disp;
  if (visualize) { visu=img.get_resize(-100,-100,-100,3); disp.assign(visu); }
  const CImg<> residuals=img.get_pseudo_residuals();
  float threshold=2.f;
  while (mean.size()<20 && threshold >=1.05f) {
    if (visualize) visu=img.get_resize(-100,-100,-100,3);
    mean.assign(img.size(),1,1,1,-1);
    var.assign(img.size(),1,1,1,-1);
    size.assign(img.size(),1,1,1,-1);
    int n=0;
    cimg_forZC(img,z,v)
      gat_quadtree_divide(0,0,img.width()-1,img.height()-1,
													img.get_shared_plane(z,v),
													residuals.get_shared_plane(z,v),
													mean,var,size,n,disp,visu,mean_method,variance_method,threshold);
    mean.crop(0,n-1);
    var.crop(0,n-1);
		size.crop(0,n-1);
    gat_cleanup(mean,var,size);
#if cimg_debug>=3
    printf("   Samples = (%d,%d) at threshold %.2f\n",(int)mean.size(),(int)var.size(),threshold+.05f);
#endif
    threshold-=0.01f;    
  }
#if cimg_debug>=3
  if (visualize) visu.save_inr("visu.inr");
#endif
  if (mean.size() < 20) {
#if cimg_debug>=3
    cimg::warn("Not enough samples in gat_compute_quadtree_stats(...)");
#endif
    gat_compute_blocks_stats(img,mean,var,mean_method,variance_method,visualize);
  }
}

void gat_compute_quadtree_stats(const CImg<> & img, CImg<> &mean, CImg<> &var,const int mean_method=1,
				const int variance_method=2,const bool visualize=false) const {
  CImg<> visu;
  CImgDisplay disp;
  if (visualize) { visu=img.get_resize(-100,-100,-100,3); disp.assign(visu); }
  const CImg<> residuals=img.get_pseudo_residuals();
  float threshold=1.5f;
  while (mean.size()<20 && threshold >=1.05f) {
    if (visualize) visu=img.get_resize(-100,-100,-100,3);
    mean.assign(img.size(),1,1,1,-1);
    var.assign(img.size(),1,1,1,-1);
    int n=0;
    cimg_forZC(img,z,v)
      gat_quadtree_divide(0,0,img.width()-1,img.height()-1,
			  img.get_shared_plane(z,v),
			  residuals.get_shared_plane(z,v),
			  mean,var,n,disp,visu,mean_method,variance_method,threshold);
    mean.crop(0,n-1);
    var.crop(0,n-1);
    gat_cleanup(mean,var);
#if cimg_debug>=3
    printf("   Samples = (%d,%d) at threshold %.2f\n",(int)mean.size(),(int)var.size(),threshold+.05f);
#endif
    threshold-=0.01f;    
  }
#if cimg_debug>=3
  if (visualize) visu.save_inr("visu.inr");
#endif
  if (mean.size() < 20) {
#if cimg_debug>=3
    cimg::warn("Not enough samples in gat_compute_quadtree_stats(...)");
#endif
    gat_compute_blocks_stats(img,mean,var,mean_method,variance_method,visualize);
  }
}

void gat_compute_blocks_stats(const CImg<> & img, CImg<> &mean, CImg<> &var,const int mean_method=1,
															const int variance_method=2,bool visualize=false) const {
  const int j0=5,j1=5,nj=j1-j0+1;
  CImg<Tfloat> residuals=img.get_pseudo_residuals();
  CImgList<Tfloat> Mean(nj),Variance(nj);
  CImg<> ctrl;
  CImgDisplay disp;
  float col1[3],col2[3];
  if (visualize) { const float mmm=img.max();
    col1[0]=mmm; col1[1]=mmm/2; col1[2]=mmm/4;
    col2[0]=255-mmm; col2[1]=255-mmm/2; col2[2]=255-mmm/4;
    disp.assign(ctrl); ctrl=img.get_resize(-100,-100,-100,3);}
  for (int j=0;j<nj;j++) {
    int bx=j0+j,by=j0+j, bz=1;  //bz=j0+j;
    int nbx=(int)(img.width()/bx), nby=(int)(img.height()/by),nbz=(int)(img.depth()>1?img.depth()/bz:1);
    int nblocks=nbx*nby*nbz*img.spectrum();
#if cimg_debug>=3
    printf("\r j=%d (%d*%d=%d) blocks (%dx%d)",j,nbx,nby,nblocks,bx,by);fflush(stdout);
#endif
    Mean[j].assign(nblocks);
    Variance[j].assign(nblocks);
    int i=0;
    ctrl=img.get_resize(-100,-100,-100,3);
    
    for (int v=0;v<img.spectrum();v++)    
      for (int zi=0;zi<img.depth();zi+=bz)
				for (int yi=0;yi<img.height()-by;yi+=by)
					for (int xi=0;xi<img.width()-bx;xi+=bx) {
						if (i<nblocks) {
							CImg<> tmp=residuals.get_crop(xi,yi,zi,v,xi+bx-1,yi+by-1,(img.depth()>1?zi+bz-1:zi),v);
							Variance[j](i)=tmp.variance(variance_method);	      
							tmp=img.get_crop(xi,yi,zi,v,xi+bx-1,yi+by-1,(img.depth()>1?zi+bz-1:zi),v);
							Mean[j](i)=tmp.median();
							/* Mean[j](i)=tmp.robust_mean_estimation(mean_method); */

#if cimg_display==1 || cimg_display==2
							if (visualize)
								ctrl.draw_rectangle(xi,yi,xi+bx-1,yi+by-1,col1,.7,~0L)
									.draw_rectangle(xi,yi,xi+bx-1,yi+by-1,col2,.25)
									.display(disp);
#endif
							i++;
						}
					}
    gat_cleanup(Mean[j],Variance[j]);
  }
  mean=Mean.get_append('x');
  var=Variance.get_append('x');
#if cimg_debug>=3
  printf("   Samples : %d\n",(int)mean.size());
#endif 
}

//! Estimate the parameters of the generalized Anscombe Transform.
/**

   The parametrization is given by : @f$Y_i=g_0 N_i+m+\sigma_0 W_i@f$
   where @f$N_i@f$ is a Poisson distributed random variable and
   @f$W_i@f$ the centered normally distributed random variable of
   variance 1.  The parameters are estimated using a regression in the
   plane @f$(E[Y],\textrm{Var}[Y])@f$ and we have @f$\textrm{Var}[Y] = g_0
   E[Y] + \sigma^2_0 - g_0 m = g_0 E[Y] + e_{DC} @f$.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m @f$
   @param m : the value of the mean value of the dark current @f$m @f$.
   @param sigma0 : a pointer to the value of the standard deviation
   @f$\sigma_0@f$ of the dark current.
   @param file_o : the filename where to plot the analysis
   @param scheme  : select one of the 3 implemented schemes for computing local mean and noise variance. 0: blocks, 1:
   level_lines, 2: quadtree.
   @param variance_method : select the method for estimating the variance.0: default, 1: unbiased 2: Least median of squares, 3:  Least Trimmed squares.
   @param robust_fit : select the fitting method: 0: least square, 1: iterative robust least square.
   @param visu : Display a visualization of the computation (debug)
**/
const CImg<T> & gat_estimate_parameters(double & g0, double & edc, double & m, double & sigma0, 
																				const char * file_o=0, const int scheme=2, 
																				const int mean_method=1, const int variance_method=3,
																				const bool robust_fit=true,const bool visu=false,
																				CImg<double> & cov=CImg<double>::empty() ) const {
#if cimg_debug>=3
  printf("   Image size : %dx%dx%d\n",width(),height(),depth());
#endif
  CImg<> X,Y;
  switch (scheme){
  case 0:  gat_compute_blocks_stats((*this),X,Y,mean_method,variance_method,visu);   break;
  case 1:  gat_compute_levellines_stats((*this),X,Y,mean_method,variance_method);    break;
  case 2:  gat_compute_quadtree_stats((*this),X,Y,mean_method,variance_method,visu); break;
  default: gat_compute_quadtree_stats((*this),X,Y,mean_method,variance_method,visu); break;
  }
	X.print("X gat",false);

  linear_fit(X,Y,edc,g0,robust_fit,cov);
  if (file_o) gat_save_graph(X,Y,g0,edc,file_o);
  double s0=0, ns=0, thresx=(X.max()-X.min())/X.size()*20.+X.min();
  cimg_foroff(X,i) if (X(i)<thresx) { s0 += Y(i); ns++; }
  if (ns!=0) sigma0 = s0/ns; else  { sigma0=Y.min()*1.5;}  
  m = (sigma0-edc)/g0;
  sigma0=std::sqrt(sigma0);  
#if cimg_debug>=3
  gat_print_stat(X,Y);
  printf("   Anscombe parameters : gain g0 : %.3g, DC : N(%.3g,%.3g)\n",g0,m,sigma0);
#endif
  return (*this);
} 


//! Generate a Poisson+Gaussian noise from an image and a noise model
/**
   The model is given by: @f$ Y_i = g_0 N_i + m + \sigma_0 W_i@f$
   with @f$ W_i@f$ follows a centered Normal distribution of
   variance 1 and @f$X_i@f$ is the Poisson random variable of
   parameter @f$g_1(u_i-m)/g_0@f$ where @f$ u_i @f$ is the input clean
   image.

   @param g0  :  the value of the gain @f$g_0@f$
   @param g1  :  the value of the photon gain @f$g_1@f$
   @param m   :  the value of the mean value of the dark current @f$m@f$.
   @param sigma0 : the value of the standard deviation  @f$\sigma_0@f$ of the dark current.   
   @see gat_estimate_parameters 
**/
CImg<T>& gat_noise(const double g0, const double g1, const double m, const double sigma0) {
  if (cimg::abs(g0)<10.0e-10) throw CImgException("gat_noise() : g0 is null.");
  cimg_for(*this,ptr,T)
    *ptr = (T)((double)cimg::prand((double)(g1*(*ptr)-m)/g0)*g0+m+sigma0*cimg::grand());
  return (*this);
}

//! Get a Poisson+Gaussian noise from an image and a noise model
/**
   The model is given by: @f$ Y_i = g_0 N_i + m + \sigma_0 W_i@f$
   with @f$ W_i@f$ follows a centered Normal distribution of
   variance 1 and @f$X_i@f$ is the Poisson random variable of
   parameter @f$g_1(u_i-m)/g_0@f$ where @f$ u_i @f$ is the input clean
   image.

   @param g0  :  the value of the gain @f$g_0@f$
   @param g1  :  the value of the photon gain @f$g_1@f$
   @param m   :  the value of the mean value of the dark current @f$m@f$.
   @param sigma0 : the value of the standard deviation  @f$\sigma_0@f$ of the dark current.   
   @see gat_estimate_parameters 
**/
CImg<T>& get_gat_noise(const double g0, const double g1, const double m, const double sigma0) const {
  return CImg<Tfloat>(*this,false).gat_noise(g0,g1,m,sigma0);
}


//! Apply the generalized Anscombe variance stabilization transform
/**
   The forward transform is given by:
   @f$ Z_i= \frac{2}{g_0} \sqrt{g_0 Y_i + \frac{3}{8} g^2_0 + e_{DC}}   @f$.

   After the transformation, if the parameters have been correctly
   estimated and that the noise model is correct, the variance of the
   noise is constant on the image and is equal to 1. This property is
   interesting for denoising for example or any further image analysis
   sensitive to the noise level.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
  
   @note A maximum to 0 is taken inside the square root to prevent
   imaginary solutions.
   @
**/
CImg<T>& gat(const double g0, const double edc) {
  if (cimg::abs(g0)<10.0e-10) throw CImgException("gat() : g0 is null.");
  //cimg_for(*this,ptr,T) 
  // (*ptr) = (T)(2.0/g0 * std::sqrt(cimg::max(g0*(double)(*ptr) + 0.375*g0*g0 + edc,0.0)));
  cimg_for(*this,ptr,T) {
    const   double a = g0*(double)(*ptr) + 0.375*g0*g0 + edc, sign=a>0?1:-1;
    (*ptr) = (T)(sign*2.0/g0*std::sqrt(cimg::max(g0*(double)(*ptr) + 0.375*g0*g0 + edc,0.0)));
  }      
  return (*this);
}

//! Get the generalized Anscombe variance stabilization transform
/**
   The forward transform is given by:
   @f$ Z_i= \frac{2}{g_0} \sqrt{g_0 Y_i + \frac{3}{8} g^2_0 + e_{DC}}   @f$.

   After the transformation, if the parameters have been correctly
   estimated and that the noise model is correct, the variance of the
   noise is constant on the image and is equal to 1. This property is
   interesting for denoising for example or any further image analysis
   sensitive to the noise level.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
  
   @note A maximum to 0 is taken inside the square root to prevent
   imaginary solutions.
   @
**/
CImg<Tfloat> & get_gat(const double g0, const double edc) const {
  return CImg<Tfloat>(*this).gat(g0,edc);
}

//! Apply the inverse generalized Anscombe variance stabilization transform
/** The inverse transform is given by:
   @f$ Y_i = g_0 Z_i^2/4 - \frac{3}{8} g_0 - e_{DC}/g_0   @f$.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
   @return The stabilize image.
**/
CImg<T> & igat(const double g0, const double edc) {
  if (is_empty()) return *this;
  if (cimg::abs(g0)<10.0e-10) throw CImgException("gat_psnr : g0 is null.");
  cimg_for(*this,ptr,T) {
    const double val=(double)(*ptr);
    (*ptr) = (T)( g0*(0.25*val*val - 0.375) - edc/g0 ); 
  }
  return (*this);
}

//! Get the inverse generalized Anscombe variance stabilization transform
/** The inverse transform is given by:
   @f$ Y_i = g_0 Z_i^2/4 - \frac{3}{8} g_0 - e_{DC}/g_0   @f$.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
   @return The stabilize image.
**/
CImg<Tfloat> & get_igat(const double g0, const double edc) const {
  return CImg<Tfloat>(*this).igat(g0,edc);
}

//! Apply the unbiased inverse generalized Anscombe variance stabilization transform
/** The inverse transform is given by:
   @f$ Y_i = g_0 Z_i^2/4 - \frac{3}{8} g_0 - e_{DC}/g_0 + \frac{1}{4}(1-\frac{1}{n})  @f$.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
   @param sigma0: the standard devition of the Gaussian noise
   @return The stabilize image.
**/
CImg<T> & uigat(const double g0, const double edc) {
  if (is_empty()) return *this;
  if (cimg::abs(g0)<10.0e-10) throw CImgException("gat_psnr : g0 is null.");
  const double bias =  - 0.25 * (1.0-1.0/(float)size());
  cimg_for(*this,ptr,T) {
    const double val=(double)(*ptr);
    (*ptr) = (T)( g0*(0.25*val*val - 0.375) - edc/g0  - bias ); 
  }
  return (*this);
}

//! Get the unbiased inverse generalized Anscombe variance stabilization transform
/** The inverse transform is given by:
   @f$ Y_i = g_0 Z_i^2/4 - \frac{3}{8} g_0 - e_{DC}/g_0  + \frac{1}{4}(1-\frac{1}{n}) @f$.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
   @param sigma0: the standard devition of the Gaussian noise
   @return The stabilize image.
**/
CImg<Tfloat> & get_uigat(const double g0, const double edc) const {
  return CImg<Tfloat>(*this).uigat(g0,edc);
}

//! Unbiase the generalized Anscombe variance stabilization transform
/** The additive bias is *empirically* estimated by computing the
   average difference with a reference (usually the original noisy
   image).  
   @param ref : reference image @return The unbiased image
**/
CImg<T> gat_unbias(const CImg<> & ref) {
  const T c=(ref-(*this)).mean(); // compute the correction
  printf("Bias correction : %f\n",c);
  return (*this)+c;
}

//! Compute the peak signal-to-noise ratio for Poisson-Gaussian noise
/**
   The peak signal-to-noise ration is defined as: @f$ \textrm{PSNR}=20
   \log_{10}(\mathcal{T}_{GA}(Max(Y_i))-\mathcal{T}_{GA}(Min(Y_i)))
   @f$ The PSNR is a image quality measure expressed in decibels
   (dB). Typical values are 40dB for a good quality and 20dB for very
   noisy images. This definition is adapted for Poisson-Gaussian noise
   and does not apply for Gaussian model.

   @param g0 : the value of the gain @f$g_0@f$
   @param edc : the value of @f$e_{DC} = \sigma^2_0 - g_0 m@f$
   @return the value of the psnr
**/
double gat_psnr(const double g0, const double edc) const {
  if (is_empty()) throw CImgException("gat_psnr : image is empty.");
  if (cimg::abs(g0)<10.0e-10) throw CImgException("gat_psnr : g0 is null.");
  double a,b=max_min(a);
  return 20.0 *  std::log10(2.0/g0*(std::sqrt(g0*b+0.375*g0*g0+edc)-std::sqrt(cimg::max(g0*a+0.375*g0*g0+edc,0.0)))); 
}


#endif 
