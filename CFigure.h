/* -*- mode:c++ -*-
** CFigure.h
** 
** Made by Jerome Boulanger
**
** Todo : legende, image, plot interpoles
**  
** Started      Fri Jun 13 10:16:25 2008 by Jerome Boulanger
** Last update  Thu, 16 Dec 2010 16:12:02 by Jerome Boulanger
*/

#ifndef   	CFIGURE_H_
# define   	CFIGURE_H_


#include <vector>
using namespace std;

#include "CImg.h"
using namespace cimg_library;

template<typename T, typename tc>
CImg<T>& draw_thick_line(CImg<T> &img, 
			 const int x0, const int y0,
			 const int x1, const int y1,
			 const tc *const color, const float opacity=1, const float line_width=0,
			 const unsigned int pattern=~0U, const bool init_hatch=true) {
  if (!color)
    throw CImgArgumentException("draw_line() : Specified color is (null).");
  if (img.is_empty()) return img;
  if (line_width<.5){
    img.draw_line(x0,y0,x1,y1,color,opacity,pattern,init_hatch);
  }
  else{
    static  unsigned int hatch = ~0U - (~0U>>1);
    if (init_hatch) hatch = ~0U - (~0U>>1);
    const float 
      dx=x1-x0,dy=y1-y0,
      D=sqrt(dx*dx+dy*dy);
    for (float i=0;i<D;i+=.5) {
      if (hatch && pattern)
	img.draw_circle((int)(x0+(i*dx)/D),(int)(y0+(i*dy)/D),line_width,
			color,opacity/(float)(line_width));
    }
  }
  return img;
}

class CFigure {
private:
  CImg<unsigned char> _board; // image where the plot are drawn 
  float _xmin,_xmax,_ymin,_ymax; // axis bounds
  float e; // space around the graph
  bool _grid_on;
  int _axis_mode; // auto or fixed
  unsigned char _background_color[3];
  unsigned char _text_color[3];
  char _title_str[1024];
  unsigned int _title_font_height;
  char _xlabel_str[1024];
  unsigned int _xlabel_font_height;
  char _ylabel_str[1024];
  unsigned int _ylabel_font_height;
  unsigned int _axis_font_height;
public:
  // a plot structure
  class plot_t{
    CImg<> _x,_y;
    unsigned int _plot_type, _vertex_type,_pattern;
    unsigned char _color[3];
    float _opacity;
    float _line_width;
  public:
    //! plot_t contructor
    template <typename tc>
    plot_t(const CImg<> x, const CImg<> y, 
					 const unsigned int plot_type, const unsigned int vertex_type, 
					 const tc * const color, const float opacity, const float line_width, const unsigned int pattern){
      _x=x;
      _y=y;
      _plot_type=plot_type;
      _vertex_type=vertex_type;
      if (color!=0) {_color[0]=color[0];_color[1]=color[1];_color[2]=color[2];}
      _opacity=opacity;
      _line_width=line_width;
      _pattern=pattern;
    }
    CImg<> & get_x(){return _x;}
    CImg<> & get_y(){return _y;}
    unsigned int get_plot_type(){return _plot_type;}
    unsigned int get_vertex_type(){return _vertex_type;}
    unsigned char * get_color(){return _color;}
    float get_opacity(){return _opacity;}
    float get_line_width(){return _line_width;}
    unsigned int get_pattern(){return _pattern;}
  };
  
  vector<plot_t> _plot_list; // list of plot (x,y)

  CFigure(){
    assign();
  }

  //! Constructor define by the width and height of the figure in pixels
  CFigure(int width, int height) {
    assign();
    _board.assign(width,height,1,3,255);
  }

  //! reinit
  CFigure & assign(){
    _xmin=0;_xmax=1;_ymin=0;_ymax=1;
    _board.assign(600,400,1,3,255);
    e=.10;
    _grid_on=true;
    _axis_mode=0;
    unsigned char white[3]={255,255,255},black[3]={0,0,0};
    set_background_color(white);
    set_text_color(black);
    sprintf(_title_str," ");
    sprintf(_xlabel_str," ");
    sprintf(_ylabel_str," ");
    _title_font_height=13;
    _xlabel_font_height=13;
    _ylabel_font_height=13;
    _axis_font_height=13;
    return *this;
  }
  
  //! Set the color of the background
  template <typename tc>
  CFigure & set_background_color(const tc * const color) {
    for (int i=0;i<3;i++) _background_color[i]=color[i];
    return *this;
  }

  //! Set the color of the text
  template <typename tc>
  CFigure & set_text_color(const tc * const color) {
    for (int i=0;i<3;i++) _text_color[i]=color[i];
    return *this;
  }

  //! erase the _board but not the list of plots
  CFigure & erase() {
    cimg_forC(_board,c)
      cimg_forXY(_board,x,y) 
      _board(x,y,c) = _background_color[c];
    return *this;
  }

  //! Clear the figure and clear the list of plots
  CFigure & clear() {
    erase();
    _plot_list.clear();
    return *this;
  }
 
  //! Set the axis
  CFigure & set_axis(float xmin, float xmax, float ymin, float ymax) {
    _xmin=xmin;_xmax=xmax;_ymin=ymin;_ymax=ymax;
    return *this;
  }
  //! Set the axis using the values of x and y
  CFigure & set_axis(const CImg<> & x,const CImg<> & y){
    _xmin=x.min();_ymin=y.min();
    _xmax=x.max();_ymax=y.max();
    if (_ymin==_ymax){_ymin--;_ymax++;}
    return *this;
  }

	float get_xmin(){return _xmin;}
	float get_ymin(){return _ymin;}
	float get_xmax(){return _xmax;}
	float get_ymax(){return _ymax;}

  //! Adjust the axis to the plots
  CFigure & set_axis() {
    if (!_plot_list.empty()) {
      _xmin=_plot_list[0].get_x().min();
      _xmax=_plot_list[0].get_x().max();
      _ymin=_plot_list[0].get_y().min()*.9;
      _ymax=_plot_list[0].get_y().max()*1.1;
      for (unsigned int i=1;i<_plot_list.size();i++) {
	_xmin=cimg::min(_plot_list[i].get_x().min(),_xmin);
	_xmax=cimg::max(_plot_list[i].get_x().max(),_xmax);
	_ymin=cimg::min(_plot_list[i].get_y().min()*.9,_ymin);
	_ymax=cimg::max(_plot_list[i].get_y().max()*1.1,_ymax);
      }
      if (_ymin==_ymax){_ymin--;_ymax++;}
    }
    return *this;
  }

  //! set the mode of the axis 1=auto or 0=predefined
  CFigure & set_axis_mode(const int mode){_axis_mode=mode;return *this;}

  CFigure & set_axis_font_size(const unsigned int font_size) {
    _axis_font_height=font_size;
    return *this;
  }

  // Format a number with a precision
  const char * select_precision(const float x){
    if (fabs(x)<0.01f) return "%.f";
    else if (fabs(x)<10.f) return "%.2f";
    return "%.0f";
  }
    
  int number_pixel_length(const float x, const unsigned int font_size=13) {
    char buf[64];
    sprintf(buf,select_precision(x),fabs(x));
    return (font_size*strlen(buf))/2;
  }

  //! Draw the axis in black
  CFigure & draw_axis() {
    
    const int 
      bx0=(int)(e*_board.width()), 
      bx1=(int)((1-e)*_board.width()),
      by0=(int)(e*_board.height()),
      by1=(int)((1-e)*_board.height());
    
    _board.draw_rectangle(bx0,by0,bx1,by1,_text_color,1,~0U);    
    int ts=10;
    for (float x=_xmin;x<=_xmax;x+=(_xmax-_xmin)/10.) {
      const int xx=(int)((x-_xmin)/(_xmax-_xmin)*(1-2*e)*_board.width()+e*_board.width());
      _board.draw_line(xx,by0,xx,by0+ts,_text_color);
      _board.draw_line(xx,by1,xx,by1-ts,_text_color);
      const int tlx=number_pixel_length(x);
      _board.draw_text(xx-tlx/2,by1+ts/2,select_precision(_xmax-_xmin),_text_color,0,1.f,_axis_font_height,x);
    }
    const int tly=number_pixel_length(_ymax);
    for (float y=_ymin;y<=_ymax;y+=(_ymax-_ymin)/10.) {
      const int yy=(int)(_board.height()- ((y-_ymin)/(_ymax-_ymin)*(1-2*e)*_board.height()+e*_board.height()));
      _board.draw_line(bx0,yy,bx0+ts,yy,_text_color);
      _board.draw_line(bx1,yy,bx1-ts,yy,_text_color);      
      _board.draw_text((int)(bx0-tly-(y<0?5:0)),yy-ts/2,select_precision(_ymax-_ymin),_text_color,0,1.f,_axis_font_height,y);
    }
    return *this;

  }

  CFigure & grid(bool v){ _grid_on=v; return *this;}
 
  CFigure & draw_grid(){
    const int 
      bx0=(int)(e*_board.width()), 
      bx1=(int)((1-e)*_board.width()),
      by0=(int)(e*_board.height()),
      by1=(int)((1-e)*_board.height());

    for (float x=_xmin;x<=_xmax;x+=(_xmax-_xmin)/10.) {
      const int xx=(int)((x-_xmin)/(_xmax-_xmin)*(1-2*e)*_board.width()+e*_board.width());
      _board.draw_line(xx,by0,xx,by1,_text_color,.25f,0x55555555);
    }
    for (float y=_ymin;y<=_ymax;y+=(_ymax-_ymin)/10.) {
      const int yy=(int)(_board.height()- ((y-_ymin)/(_ymax-_ymin)*(1-2*e)*_board.height()+e*_board.height()));
      _board.draw_line(bx0,yy,bx1,yy,_text_color,.25f,0x55555555);
    }
    return *this;
  }

  //! Set a title
  CFigure & title(const char * text, const unsigned int font_size=13) {
    if (text!=0) std::strcpy(_title_str,text);
    _title_font_height=font_size;
    return *this;
  }
  
  CFigure & draw_title() {
    if (_title_str!=0) {
      const int l=strlen(_title_str);
      if (l<1024 && l>0){
	const int x0=(int)(_board.width()/2-6*l/2),y0=(int)(e*_board.height())/2;
	_board.draw_text(x0,y0,"%s",_text_color,0,1.f,_title_font_height,_title_str);
      }
    }
    return *this;
  }

  //! Set a label for X-axis
  CFigure & xlabel(const char * text,const unsigned int font_size=13){
    if (text!=0) std::strcpy(_xlabel_str,text);
    _xlabel_font_height=font_size;
    return *this;
  }
  
  CFigure & draw_xlabel() {
    if (_xlabel_str!=0) {
      if (strcmp(_xlabel_str," ")!=0){
	const int l=strlen(_xlabel_str);
	if (l<1024 && l>0){
	  const int x0=(int)(_board.width()/2-6*l/2),y0=(int)((1-e)*_board.height()+_board.height()  )/2;
	  _board.draw_text(x0,y0,"%s",_text_color,0,1.f,_xlabel_font_height,_xlabel_str);
	}
      }
    }    
    return *this;
  }

  //! Set a label for Y-axis
  CFigure & ylabel(const char * text, const unsigned int font_size=13){
    if (text!=0)  std::strcpy(_ylabel_str,text);
    _ylabel_font_height=font_size;
    return *this;
  }
  
  CFigure & draw_ylabel(){
    if (_ylabel_str!=0) {
      const int sl=strlen(_ylabel_str);
      if (sl<1024 && sl>0){
	const int l=number_pixel_length(sl,_ylabel_font_height);
	CImg<unsigned char> tmp(_board.height()-2*e*_board.height(),_ylabel_font_height,1,3);
	cimg_forC(tmp,c)
	  cimg_forXY(tmp,x,y) 
	  tmp(x,y,c) = _background_color[c];
	tmp.draw_text(tmp.width()/2-l/2,0,"%s",_text_color,0,1.f,_ylabel_font_height,_ylabel_str);
	tmp.rotate(-90);
	const int  x0=(int)(0), y0= (int)(e*_board.height());
	_board.draw_image(x0,y0,tmp);
      }
    }    
    return *this;
  }
 
  // Convert from real coordinates to _board coordinate
  int _x(const float x) {
    return (int)((x-_xmin)/(_xmax-_xmin)*(1-2*e)*_board.width()+e*_board.width());
  }
  
  // Convert from real coordinates to _board coordinate
  int _y(const float y) {
    return (int)(_board.height()- ((y-_ymin)/(_ymax-_ymin)*(1-2*e)*_board.height()+e*_board.height()));
  }

  // Convert from _board coordiante to real
  float _ix(const int x){
    return (x-e*_board.width())/((1-2*e)*_board.width())*(_xmax-_xmin)+_xmin;
  }

  float _iy(const int y){
    return (_board.height()-y-e*_board.height())/((1-2*e)*_board.height())*(_ymax-_ymin)+_ymin;
  }

  bool is_inside(const float x, const float y){
    return x>=_xmin && x<=_xmax && y>=_ymin && y<=_ymax;
  }

  // Draw a point
  template<typename tc>
  CFigure & draw_point(const float x, const float y,
		       const tc *const color, const float opacity=1) {
    if (is_inside(x,y))
      _board.draw_point(_x(x),_y(y),color,opacity);
    return *this;
  }
  
  // Draw a cross
  template<typename tc>
  CFigure & draw_cross(const float x, const float y, const int size,
		       const tc *const color=0, const float opacity=1) {
    if (is_inside(x,y)){
      const int X=_x(x),Y=_y(y);
      _board.draw_line(X-size,Y,X+size,Y,color,opacity);
      _board.draw_line(X,Y-size,X,Y+size,color,opacity);
    }
    return *this;
  }

  //! Draw a x-cross
  template<typename tc>
  CFigure & draw_xcross(const float x, const float y, const int size,
			const tc *const color, const float opacity=1) {
    if (is_inside(x,y)){
      const int X=_x(x),Y=_y(y);    
      _board.draw_line(X-size,Y-size,X+size,Y+size,color,opacity);
      _board.draw_line(X-size,Y+size,X+size,Y-size,color,opacity);
    }
    return *this;
  }

  // Draw a circle
  template<typename tc>
  CFigure & draw_circle(const float x, const float y, const int size,
			const tc *const color, const float opacity=1) {
    if (is_inside(x,y)){
      _board.draw_circle(_x(x),_y(y),size,color,opacity);
    }
    return *this;
  }

  // Draw a empty circle 
  template<typename tc>
  CFigure & draw_circle(const float x, const float y, const int size,
			const tc *const color, const float opacity,
			const unsigned int) {
    if (is_inside(x,y)){
      _board.draw_circle(_x(x),_y(y),size,color,opacity,1);
    }
    return *this;
  }

  // Draw a full square
  template<typename tc>
  CFigure & draw_square(const float x, const float y, const int size,
			const tc *const color, const float opacity=1) {
    if (is_inside(x,y)){
      const int X=_x(x),Y=_y(y);
      _board.draw_rectangle(X-size,Y-size,X+size,Y+size,color,opacity);
    }
    return *this;
  }

  // Draw an empty diamond
  template<typename tc>
  CFigure & draw_diamond(const float x, const float y, const int size,
			 const tc *const color, const float opacity=1) {
    if (is_inside(x,y)){
      const int X=_x(x),Y=_y(y);
      _board.draw_line(X,Y+size,X+size,Y,color,opacity);
      _board.draw_line(X-size,Y,X,Y+size,color,opacity);
      _board.draw_line(X-size,Y,X,Y-size,color,opacity);
      _board.draw_line(X,Y-size,X+size,Y,color,opacity);
    }
    return *this;
  }

  // Draw a line
  template<typename tc>
  CFigure & draw_line(const float x0, const float y0, const float x1, const float y1,
		      const tc *const color, const float opacity=1, const float line_width=0,
		      const unsigned int pattern=~0U,bool init_hatch=true) {
    if (is_inside(x0,y0) && is_inside(x1,y1))
      draw_thick_line(_board,_x(x0),_y(y0),_x(x1),_y(y1),color,opacity,line_width,pattern,init_hatch);    
    return *this;
  }

  // Draw an empty rectangle
  template<typename tc>
  CFigure & draw_rectangle(const float x0, const float y0, const float x1, const float y1,
			   const tc *const color, const float opacity, const unsigned int pattern) {
    if (is_inside(x0,y0) && is_inside(x1,y1))
      _board.draw_rectangle(_x(x0),_y(y0),_x(x1),_y(y1),color,opacity,pattern);    
    return *this;
  }
  
  // Draw an filled rectangle
  template<typename tc>
  CFigure & draw_rectangle(const float x0, const float y0, const float x1, const float y1,
			   const tc *const color, const float opacity=1.f){			   
    if (is_inside(x0,y0) && is_inside(x1,y1))
      _board.draw_rectangle(_x(x0),_y(y0),_x(x1),_y(y1),color,opacity);    
    return *this;
  }
  
  // Draw an empty bar
  template<typename tc>
  CFigure & draw_bar(const float x0, const float y0, const float x1, const float y1,
		     const tc *const color, const float opacity,
		     const unsigned int pattern) {
    float y00=cimg::max(_ymin,(cimg::min(_ymax,y0)));
    cimg::unused(y1);
    if (is_inside(x0,y00) && is_inside(x1,y00)){
      _board.draw_rectangle(_x(x0),_y(y00),_x(x1),_y(0),color,opacity,pattern);
    }
    return *this;
  }

  // Draw an filled bar
  template<typename tc>
  CFigure & draw_bar(const float x0, const float y0, const float x1, const float y1,
		     const tc *const color, const float opacity=1.f) {
    float y00=cimg::max(_ymin,(cimg::min(_ymax,y0)));
    cimg::unused(y1);
    if (is_inside(x0,y00) && is_inside(x1,y00)){
      _board.draw_rectangle(_x(x0),_y(y00),_x(x1),_y(0),color,opacity);
    }
    return *this;
  }

  // Draw some text at the given location
  template<typename tc>
  CFigure & draw_text(const float x0, const float y0, const char *const text,
		      const tc *const foreground_color, const int background_color=0,
		      const float opacity=1, const unsigned int font_height=13, ...){
    if (is_inside(x0,y0)) {
      char tmp[2048] = { 0 };  std::va_list ap; 
      va_start(ap,font_height); cimg_vsnprintf(tmp,sizeof(tmp),text,ap); va_end(ap);
      _board.draw_text(_x(x0),_y(y0),tmp,foreground_color, background_color,opacity,font_height);
    }
    return *this;
  }
  
  //! Higher-level interface for plot
  /**
     Supports matlab-like syntax for simple plot type "r", "r-", "r+", "rx" "ro" for colors r,g,b,k
  **/
  template<typename tx, typename ty>
  CFigure & plot(const CImg<tx> &x, const CImg<ty> &y, const char * plot_type){
    if (!cimg::strcasecmp(plot_type,"r")||!cimg::strcasecmp(plot_type,"r-")) {
      plot(x,y,1,0,rgb8("red")); return *this; }
    if (!cimg::strcasecmp(plot_type,"r+")) {plot(x,y,0,2,rgb8("red")); return *this; }
    if (!cimg::strcasecmp(plot_type,"rx")) {plot(x,y,0,3,rgb8("red")); return *this; }
    if (!cimg::strcasecmp(plot_type,"ro")) {plot(x,y,0,5,rgb8("red")); return *this; }
    if (!cimg::strcasecmp(plot_type,"b")||!cimg::strcasecmp(plot_type,"b-")) {
      plot(x,y,1,0,rgb8("blue")); return *this; }
    if (!cimg::strcasecmp(plot_type,"b+")) {plot(x,y,0,2,rgb8("blue")); return *this; }
    if (!cimg::strcasecmp(plot_type,"bx")) {plot(x,y,0,3,rgb8("blue")); return *this; }
    if (!cimg::strcasecmp(plot_type,"bo")) {plot(x,y,0,5,rgb8("blue")); return *this; }
    if (!cimg::strcasecmp(plot_type,"g")||!cimg::strcasecmp(plot_type,"g-")) {
      plot(x,y,1,0,rgb8("green")); return *this; }
    if (!cimg::strcasecmp(plot_type,"g+")) {plot(x,y,0,2,rgb8("green")); return *this; }
    if (!cimg::strcasecmp(plot_type,"gx")) {plot(x,y,0,3,rgb8("green")); return *this; }
    if (!cimg::strcasecmp(plot_type,"go")) {plot(x,y,0,5,rgb8("green")); return *this; }
    if (!cimg::strcasecmp(plot_type,"k")||!cimg::strcasecmp(plot_type,"k-")) {
      plot(x,y,1,0,rgb8("black")); return *this; }
    if (!cimg::strcasecmp(plot_type,"k+")) {plot(x,y,0,2,rgb8("black")); return *this; }
    if (!cimg::strcasecmp(plot_type,"kx")) {plot(x,y,0,3,rgb8("black")); return *this; }
    if (!cimg::strcasecmp(plot_type,"ko")) {plot(x,y,0,5,rgb8("black")); return *this; }
    if (!cimg::strcasecmp(plot_type,"y")||!cimg::strcasecmp(plot_type,"y-")) {
      plot(x,y,1,0,rgb8("yellow")); return *this; }
    if (!cimg::strcasecmp(plot_type,"y+")) {plot(x,y,0,2,rgb8("yellow")); return *this; }
    if (!cimg::strcasecmp(plot_type,"yx")) {plot(x,y,0,3,rgb8("yellow")); return *this; }
    if (!cimg::strcasecmp(plot_type,"yo")) {plot(x,y,0,5,rgb8("yellow")); return *this; }
    printf("Unsupported plot type.\n");
    return *this;
  }

  //! Make a bar plot (histogram like)
  /**
     \param x : X-axis values
     \param y : Y-axis values
     \param color : a pointer to a 3 elements c-array.
     \param opacity : a float [0,1] defining the opacity
  **/
  template<typename tx, typename ty,typename tc>
  CFigure & bar(const CImg<tx> &x, const CImg<ty> &y, const tc *const color, const float opacity=1.0f){
    plot(x,y,4,0,color,opacity);
    return *this;
  }

  //! Plot the data y as a function of x
  /**
     \param x : X-axis values
     \param y : Y-axis values
     \param plot_type : [0,5] type of line (none,segments,bar,vertical lines)
     \param vertex_type : [0,7] type of glyphe
     \param color : a pointer to a 3 elements c-array.
     \param opacity : a float [0,1] defining the opacity
     \param pattern set a pattern using a unsigned long code
  **/
  template<typename tx, typename ty, typename tc>
  CFigure & plot(const CImg<tx> &x, const CImg<ty> &y,
		 const unsigned int plot_type, const int vertex_type,
		 const tc *const color, const float opacity=1.f, 
		 const float line_width=0., const unsigned int pattern=~0U) {
    if (x.size()!=y.size()) throw CImgException("CFigure::plot() x and y have not the same size.");
    if (x.is_empty()) throw CImgException("CFigure::plot() x is empty.");
    if (y.is_empty()) throw CImgException("CFigure::plot() x is empty.");    
    plot_t myplot(x,y,plot_type,vertex_type,color,opacity,line_width,pattern);
    _plot_list.push_back(myplot);
    replot();
    return *this;
  }

  //! Plot a point
  /**
     \param x : X-axis values
     \param y : Y-axis values
     \param plot_type : [0-4] type of line (none,segments,bar,vertical lines)
     \param vertex_type : [0,7] type of glyphe
     \param color : a pointer to a 3 elements c-array.
     \param opacity : a float [0,1] defining the opacity
     \param pattern set a pattern using a unsigned long code
     \note This is use to place isolated points on the graph with differents styles.
  **/
  template<typename tc>
  CFigure & plot(float x, float y,
		 const unsigned int plot_type, const int vertex_type,
		 const tc *const color, const float opacity=1.f){
    CImg<> vx(1,1,1,1,x),vy(1,1,1,1,y);
    plot(vx,vy,plot_type,vertex_type,color,opacity);
    return *this;
  }  

  //! Re-plot all the plots contained in the Figure.
  CFigure & replot() {
    try{
      if (_axis_mode==1) set_axis();
      if (_grid_on) draw_grid();
    
      for (unsigned int i=0;i<_plot_list.size();i++) {
	CImg<> x(_plot_list[i].get_x());
	CImg<> y(_plot_list[i].get_y());
	unsigned int plot_type=_plot_list[i].get_plot_type();
	unsigned int vertex_type=_plot_list[i].get_vertex_type();
	unsigned char * color=_plot_list[i].get_color();
	float opacity = _plot_list[i].get_opacity();
	float line_width = _plot_list[i].get_line_width();
	unsigned int pattern = _plot_list[i].get_pattern();

	switch (plot_type) {
	case 0:{break;}
	case 1:{ // lines
	  bool init_hatch = true;
	  for (int i=0;i<(int)x.size()-1;++i) {
	    draw_line(x(i),y(i),x(i+1),y(i+1),color,opacity,line_width,pattern,init_hatch);
	    init_hatch = false;
	  } 
	} break;
	case 2: { // Spline
	  // not implemented
	} break;
	case 3:{ // bars
	  for (int i=0;i<(int)x.size()-1;++i) 
	    draw_bar(x(i),y(i),x(i+1),y(i+1),color,opacity,pattern);
	} break;
	case 4:{ // filled bars
	  for (int i=0;i<(int)x.size()-1;++i) 
	    draw_bar(x(i),y(i),x(i+1),y(i+1),color,opacity);
	} break;	
	case 5:{ // vertical lines
	  cimg_foroff(x,i)
	    draw_line(x(i),y(i),x(i),cimg::min(_ymax,cimg::max(_ymin,0)),color,opacity,line_width,pattern);
	} break;
	default:{
	  printf("CFigure::plot : plot type not taken into account\n");
	}
	}

	switch (vertex_type) {
	case 0:{break;}
	case 1 : { // a point
	  cimg_foroff(x,i) 
	    draw_point(x(i),y(i),color,opacity); 
	} break;
	case 2 : {  // a + cross
	  cimg_foroff(x,i) 
	    draw_cross(x(i),y(i),2,color,opacity); 
	} break;
	case 3 : { // a x cross
	  cimg_foroff(x,i) 
	    draw_xcross(x(i),y(i),2,color,opacity); 
	} break;
	case 4 : {  // filled circle
	  cimg_foroff(x,i) 
	    draw_circle(x(i),y(i),2,color,opacity); 
	} break;
	case 5 : { // oulined circle
	  cimg_foroff(x,i) 
	    draw_circle(x(i),y(i),2,color,opacity,1); 
	} break;
	case 6 : { // square
	  cimg_foroff(x,i) 
	    draw_square(x(i),y(i),2,color,opacity); 
	} break;
	case 7 : { // diamond
	  cimg_foroff(x,i) 
	    draw_diamond(x(i),y(i),2,color,opacity); 
	} break;
	default:{
	  printf("CFigure::plot : vertex type not taken into account\n");
	}
	}
      }    
      draw_axis();
      draw_ylabel();
      draw_xlabel();
      draw_title();
    }catch(CImgException e){
      printf("%s\n",e.what());
      throw CImgException("CFigure::replot()");
    }
    return *this;
  }

  CFigure & resize(CImgDisplay & disp){
    _board.resize(disp.width(),disp.height(),1,3);
    return *this;
  }

  //! Return a reference to the CImg<uchar> of the plot;
  CImg<unsigned char> & get_image(){
    return _board;
  }

  //! Display the figure on a CImgDislay
  CFigure & display(CImgDisplay & disp){
    if (!_board.is_empty()){
      if (disp.is_resized()) { disp.resize(); _board.resize(disp); replot(); }
      _board.display(disp);      
    }
    else
      printf("CFigure: Empty image.\n");
    return *this;
  }


  //! Display the figure and provide interactions
  CFigure & display_interactive(CImgDisplay & disp) {
    const unsigned char gray[3]={128,128,128},gray2[3]={128,128,192};
    static float xmin0=_xmin,ymin0=_ymin,xmax0=_xmax,ymax0=_ymax;
    static int phase=0;
    static float x0=_xmin,y0=_ymin;   
    int X=disp.mouse_x(),Y=disp.mouse_y();
    float x=_ix(X),y=_iy(Y);
    unsigned int button=disp.button();
    if (button&1){
      if (phase==0){ phase=1; x0=x;y0=y;} 
    }
    if (phase==1&& !(button&1)) { // update the range
      phase=0;
      if (fabs(x0-x)>.05*(_xmax-_xmin) && fabs(y0-y)>.05*(_ymax-_ymin)) { // if not too small
	_xmin=cimg::min(x0,x); _ymin=cimg::min(y0,y);
	_xmax=cimg::max(x0,x); _ymax=cimg::max(y0,y);
      }else { // re-init the range
	_xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;
      }
    }
    switch(disp.key()){
    case cimg::keyHOME:{_xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;} break;
    case cimg::keyS:{
      if (disp.is_keyCTRLLEFT()){
	_board.save("figure.png");
	printf("CFigure: figure.png saved\n");
      }
    } break;
    }
    // Drawing
    _board.assign(disp.width(),disp.height(),1,3);
    erase();
    replot();
    if (is_inside(x,y)) {
      switch(phase) {
      case 0: {
	draw_line(x,_ymin,x,_ymax,gray,1,0xCCCCCCCCU);
	draw_line(_xmin,y,_xmax,y,gray,1,0xCCCCCCCCU);
	_board.draw_text(5,5,"(%g, %g)",gray,0,.8f,13,x,y);
      } break;
      case 1: {
	const float 
	  xa=_x(cimg::min(x0,x)),ya=cimg::min(_y(y0),_y(y)),
	  xb=_x(cimg::max(x0,x)),yb=cimg::max(_y(y0),_y(y));
	cimg_forXYC(_board,xi,yi,ci){
	  if (xi>e*_board.width()&&xi<_board.width()*(1-e) &&
	      yi>e*_board.height()&&yi<_board.height()*(1-e)&&
	    !(xi>xa && xi<xb && yi>ya && yi<yb)) _board(xi,yi,ci)=.25*_board(xi,yi,ci)+.75*_background_color[ci];
      }	  
	draw_rectangle(x0,y0,x,y,gray2,.25);
	draw_rectangle(x0,y0,x,y,gray,1,0xCCCCCCCCU);
	char buf[1024];
	sprintf(buf,"(%.2f, %.2f)-(%.2f, %.2f)",x0,y0,x,y);
	_board.draw_text(.5*(xa+xb)-2*number_pixel_length(strlen(buf),13),.5*(ya+yb),buf,gray,_background_color,1,13);	
      } break;
    }
  }
    display(disp);
		
  
  return *this;
}

  //! Display the figure and provide interactions
	// float xmin0,ymin0,xmax0,ymax0;
	// int phase;
	// float x0,y0;   
  CFigure & display_interactive_debug(CImgDisplay & disp) {
    const unsigned char gray[3]={128,128,128},gray2[3]={128,128,192};
    float xmin0=_xmin,ymin0=_ymin,xmax0=_xmax,ymax0=_ymax;
    int phase=0;
    float x0=_xmin,y0=_ymin;
      int X=disp.mouse_x(),Y=disp.mouse_y();
      float x=_ix(X),y=_iy(Y);
      unsigned int button=disp.button();
      if (button&1){
	if (phase==0){ phase=1; x0=x;y0=y;} 
      }
      if (phase==1&& !(button&1)) { // update the range
	phase=0;
	if (fabs(x0-x)>.05*(_xmax-_xmin) && fabs(y0-y)>.05*(_ymax-_ymin)) { // if not too small
	  _xmin=cimg::min(x0,x); _ymin=cimg::min(y0,y);
	  _xmax=cimg::max(x0,x); _ymax=cimg::max(y0,y);
	}else { // re-init the range
	  _xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;
	}
      }
      switch(disp.key()){
      case cimg::keyHOME:{_xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;} break;
      case cimg::keyS:{
	if (disp.is_keyCTRLLEFT()){
	  _board.save("figure.png");
	  printf("CFigure: figure.png saved\n");
	}
      } break;
      }
      // Drawing
      _board.assign(disp.width(),disp.height(),1,3);
      erase();
      replot();
      if (is_inside(x,y)) {
	switch(phase) {
	case 0: {
	  draw_line(x,_ymin,x,_ymax,gray,1,0xCCCCCCCCU);
	  draw_line(_xmin,y,_xmax,y,gray,1,0xCCCCCCCCU);
	  _board.draw_text(5,5,"(%g, %g)",gray,0,.8f,13,x,y);
	} break;
	case 1: {

	  const float 
	    xa=_x(cimg::min(x0,x)),ya=cimg::min(_y(y0),_y(y)),
	    xb=_x(cimg::max(x0,x)),yb=cimg::max(_y(y0),_y(y));
	  cimg_forXYC(_board,xi,yi,ci){
	    if (xi>e*_board.width()&&xi<_board.width()*(1-e) &&
		yi>e*_board.height()&&yi<_board.height()*(1-e)&&
		!(xi>xa && xi<xb && yi>ya && yi<yb)) _board(xi,yi,ci)=.25*_board(xi,yi,ci)+.75*_background_color[ci];
	  }	  
	  draw_rectangle(x0,y0,x,y,gray2,.25);
	  draw_rectangle(x0,y0,x,y,gray,1,0xCCCCCCCCU);
	  char buf[1024];
	  sprintf(buf,"(%.2f, %.2f)-(%.2f, %.2f)",x0,y0,x,y);
	  _board.draw_text(.5*(xa+xb)-2*number_pixel_length(strlen(buf),13),.5*(ya+yb),buf,gray,_background_color,1,13);	
	} break;
	}
      }
      display(disp);
    return *this;
  }


  //! Display the figure in a stand alone display.
  CFigure & display(const char * window_title="") {
    CImgDisplay disp(_board);
    disp.set_title(window_title);
    const unsigned char gray[3]={128,128,128},gray2[3]={128,128,192};
    float xmin0=_xmin,ymin0=_ymin,xmax0=_xmax,ymax0=_ymax;
    int phase=0;
    float x0=_xmin,y0=_ymin;
    while (!disp.is_closed()) {
      int X=disp.mouse_x(),Y=disp.mouse_y();
      float x=_ix(X),y=_iy(Y);
      unsigned int button=disp.button();
      if (button&1){
	if (phase==0){ phase=1; x0=x;y0=y;} 
      }
      if (phase==1&& !(button&1)) { // update the range
	phase=0;
	if (fabs(x0-x)>.05*(_xmax-_xmin) && fabs(y0-y)>.05*(_ymax-_ymin)) { // if not too small
	  _xmin=cimg::min(x0,x); _ymin=cimg::min(y0,y);
	  _xmax=cimg::max(x0,x); _ymax=cimg::max(y0,y);
	}else { // re-init the range
	  _xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;
	}
      }
      switch(disp.key()){
      case cimg::keyHOME:{_xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;} break;
      case cimg::keyS:{
	if (disp.is_keyS()){
	  _board.save("figure.png");
	  printf("CFigure: figure.png saved\n");
	}
      } break;
      }
      // Drawing
      _board.assign(disp.width(),disp.height(),1,3);
      erase();
      replot();
      if (is_inside(x,y)) {
	switch(phase) {
	case 0: {
	  draw_line(x,_ymin,x,_ymax,gray,1,0xCCCCCCCCU);
	  draw_line(_xmin,y,_xmax,y,gray,1,0xCCCCCCCCU);
	  _board.draw_text(5,5,"(%g, %g)",gray,0,.8f,13,x,y);
	} break;
	case 1: {

	  const float 
	    xa=_x(cimg::min(x0,x)),ya=cimg::min(_y(y0),_y(y)),
	    xb=_x(cimg::max(x0,x)),yb=cimg::max(_y(y0),_y(y));
	  cimg_forXYC(_board,xi,yi,ci){
	    if (xi>e*_board.width()&&xi<_board.width()*(1-e) &&
		yi>e*_board.height()&&yi<_board.height()*(1-e)&&
		!(xi>xa && xi<xb && yi>ya && yi<yb)) _board(xi,yi,ci)=.25*_board(xi,yi,ci)+.75*_background_color[ci];
	  }	  
	  draw_rectangle(x0,y0,x,y,gray2,.25);
	  draw_rectangle(x0,y0,x,y,gray,1,0xCCCCCCCCU);
	  char buf[1024];
	  sprintf(buf,"(%.2f, %.2f)-(%.2f, %.2f)",x0,y0,x,y);
	  _board.draw_text(.5*(xa+xb)-2*number_pixel_length(strlen(buf),13),.5*(ya+yb),buf,gray,_background_color,1,13);	
	} break;
	}
      }
      display(disp);
      disp.wait();      
    }
    // Re-init the range of the plot for the next plots.
    _xmin=xmin0; _ymin=ymin0; _xmax=xmax0; _ymax=ymax0;
    return *this;
  }
  
CFigure & save(const char *const filename) {
  _board.save(filename);
  return *this;
}
 
  // Provide color codes
  static const unsigned char * rgb8(const char* name) {
    static const unsigned char red[3]={255,0,0};
    static const unsigned char green[3]={0,255,0};
    static const unsigned char blue[3]={0,0,255};
    static const unsigned char orange[3]={255,175,0};
    static const unsigned char pink[3]={255,20,147};
    static const unsigned char purple[3]={147,112,219};
    static const unsigned char black[3]={0,0,0};
    static const unsigned char gray128[3]={128,128,128};
    static const unsigned char white[3]={255,255,255};
    static const unsigned char yellow[3]={255,255,0};
    static const unsigned char navy[3]={0,0,128};
    static const unsigned char darkgreen[3]={0,0,128};
    static const unsigned char maroon[3]={128,0,0};
    if (!cimg::strcasecmp(name,"red")) return red;
    if (!cimg::strcasecmp(name,"green")) return green;
    if (!cimg::strcasecmp(name,"blue")) return blue;
    if (!cimg::strcasecmp(name,"orange")) return orange;
    if (!cimg::strcasecmp(name,"pink")) return pink;
    if (!cimg::strcasecmp(name,"purple")) return purple;
    if (!cimg::strcasecmp(name,"yellow")) return yellow;
    if (!cimg::strcasecmp(name,"navy")) return navy;
    if (!cimg::strcasecmp(name,"darkgreen")) return darkgreen;
    if (!cimg::strcasecmp(name,"maroon")) return maroon;
    if (!cimg::strcasecmp(name,"black")) return black;
    if (!cimg::strcasecmp(name,"gray128")) return gray128;
    if (!cimg::strcasecmp(name,"white")) return white;  
    else return 0;
  }

};

#endif 	    /* !CFIGURE_H_ */
