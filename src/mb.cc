#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <cmath>
#include <complex>
#include <assert.h>
using namespace std;


float floatMod(float a, float b)
{
    return (a - b * floor(a / b));
}

float dumbRound(float x){

  float low = floor(x);
  if ((x - low) > 0.5) {
    return ceil(x);
  } else {
    return low;
  }

}

/*takes an array of row major, contiguous rbg values as unsigned char (0-255) and writes them
 * to a ppm file*/
void write_to_p6(char* filename,int dim_x, int dim_y, unsigned char data []) {
	int length = 3*dim_x*dim_y; //colors
	int max_value = 255;
  	ofstream f;
  	f.open(filename);
  	f << "P6\n";
	f << dim_x << " " << dim_y << endl; //sprintf(buff,"%d %d\n",dim_x,dim_y);
	f << max_value << endl; 
	f.close();
	f.open(filename,ios::app | ios::binary);
	f.write((char *) data, length);
}


struct hsv_color {
	double h;
	double s;
	double v;
};

struct rgb_color {
	unsigned char r;
	unsigned char g;
	unsigned char b;
};

rgb_color hsv_to_rgb(double h, double s, double v) {
	h = 360*h;
	
	double c = v*s;
	double hp = h/60;
	double x = c*(1 - fabs(floatMod(hp,2.0)-1));
	double r,g,b;
	rgb_color rgb;
	if (0 <= hp && hp < 1) {r = c; g = x; b = 0;}
	if (1 <= hp && hp < 2) {r = x; g = c; b = 0;}
	if (2 <= hp && hp < 3) {r = 0; g = c; b = x;}
	if (3 <= hp && hp < 4) {r = 0; g = x; b = c;}
	if (4 <= hp && hp < 5) {r = x; g = 0; b = c;}
	if (5 <= hp && hp <= 6) {r = c; g = 0; b = x;}
	double m = v-c;
	r += m;
	g += m;
	b += m;
	rgb.r = (unsigned char) dumbRound(255*r);
	rgb.g = (unsigned char) dumbRound(255*g);
	rgb.b = (unsigned char) dumbRound(255*b);
	return rgb;
}

rgb_color get_color(double distance,int iter,int iter_max, double pixel_size, double radius, double radius_max) {
	rgb_color ret;
	double hue,val;
	double sat = 0.7;
	if (iter >= iter_max) {
		ret.r = 255; ret.g = 255; ret.b = 255	;
		return ret;
	} 
	if (distance < 0.5*pixel_size) {
		val = pow(distance/(0.5*pixel_size),1.0/3);
	}	else {
		val = 1.0;
	}	
	double iter_cont = iter - log2(log(radius)/log(radius_max));
	hue = log(iter_cont)/log(iter_max);
	hue = hue*10;
	hue = hue - floor(hue);
	return hsv_to_rgb(hue,sat,val);
}


double mult_complex_re(double x1,double y1,double x2, double y2) { return x1*x2 - y1*y2;}
double mult_complex_im(double x1,double y1,double x2, double y2) { return x1*y2 + y1*x2;}

double square_complex_re(double x,double y) { return x*x - y*y;}
double square_complex_im(double x,double y) { return 2*x*y;}
double abs_complex(double x,double y) {return pow(x*x + y*y,0.5);}


rgb_color Mandelbrot(double x0, double y0,double pixel_size) {

  int iter = 0;
  int iter_max = 10000;
  double radius_max = 1.0* ( 1 << 18);

  double radius = 0.0;
  double x = 0.0;
  double x_tmp = 0.0;
  double dx_tmp = 0.0;
  double y = 0.0;
	//complex<double> z (0.0,0.0);
	//complex<double> c  (x0,y0);//x0 + std::complex::complex_literals::i*y0;
  double dx = 0.0;
  double dy = 0.0;

  
  while (radius < radius_max && iter < iter_max) {
		dx_tmp = dx;
		dx = 2*mult_complex_re(x,y,dx,dy) + 1;
		dy = 2*mult_complex_im(x,y,dx_tmp,dy) + 1;
		//z = z*z + c;
		//radius = abs(z);

		x_tmp = x;
		x = square_complex_re(x,y) + x0;
		y = square_complex_im(x_tmp,y) + y0;


		radius = abs_complex(x,y);
	  iter++;
  }
	//rgb_color ret;
	//if (iter < iter_max) {ret.r = 255;ret.g = 255;ret.b = 255;}
	//else {ret.r = 0;ret.g = 0;ret.b = 0;}
	
	double distance = 2*log(radius)*radius/abs_complex(dx,dy);
	return get_color(distance,iter,iter_max,pixel_size,radius,radius_max);
}

void downsample_pixels(rgb_color* pixels_hr,rgb_color* pixels_lr,int fac,int pixel_count_x,int pixel_count_y) {
	int accum_r = 0;
	int accum_g = 0;
	int accum_b = 0;
	int pixel_hr_x;
	int pixel_hr_y;
	rgb_color col;
	rgb_color ret;
		
  for (int pixel_y=0; pixel_y<pixel_count_y/fac; pixel_y++) {
		//#pragma acc parallel
    for (int pixel_x=0; pixel_x<pixel_count_x/fac; pixel_x++) {
			accum_r= 0;
			accum_g= 0;
			accum_b= 0;
			for (int j =0; j < fac; j++) {
				for (int i=0; i < fac; i++) {
					pixel_hr_y = fac*pixel_y+j;
					pixel_hr_x = fac*pixel_x+j;
					col = pixels_hr[pixel_hr_y*pixel_count_x+pixel_hr_x];
					accum_r += col.r;
					accum_g += col.g;
					accum_b += col.b;
				}
			}	
			accum_r /= fac*fac;
			accum_g /= fac*fac;
			accum_b /= fac*fac;
			assert(accum_r <= 255);
			assert(accum_r <= 255);
			assert(accum_r <= 255);
			ret.r = (unsigned char) accum_r;
			ret.g = (unsigned char) accum_g;
			ret.b = (unsigned char) accum_b;

      pixels_lr[pixel_y*pixel_count_x/fac+pixel_x] = ret;
		}
	}

}


int main () {

	
	int super_res = 5;
	int pixel_count_raw = 1000;
  int pixel_count_x = super_res*pixel_count_raw;

  double center_x = -.74364386269; //-0.75
  double center_y = 0.13182590271;//0.00
  double length_x = 0.00000013526; //2.75
  double length_y = 0.00000013526; //2.0


  double pixel_size = length_x / pixel_count_x;
	int pixel_count_y = floor(length_y/pixel_size);
 	length_y = pixel_count_y*pixel_size;

  double minx = center_x - length_x/2.0;
  double maxy = center_y + length_y/2.0; 
  

  rgb_color * pixels = (rgb_color *) malloc( sizeof(rgb_color)*pixel_count_x*pixel_count_y );
  rgb_color * out_pixels = (rgb_color *) malloc( sizeof(rgb_color)*pixel_count_x*pixel_count_y/(super_res*super_res));

	int buffer_size = pixel_count_x*pixel_count_y;
	#pragma acc data copyout(pixels[0:buffer_size])
	#pragma acc parallel
	{
	#pragma acc loop independent
  for (int pixel_y=0; pixel_y<pixel_count_y; pixel_y++) {
		//#pragma acc parallel
		#pragma acc loop independent
    for (int pixel_x=0; pixel_x<pixel_count_x; pixel_x++) {

      double x = minx + pixel_x*pixel_size;
      double y = maxy - pixel_y*pixel_size;

			//double x = minx + i*pixel_size;
      //double y = maxy - j*pixel_size;

      pixels[pixel_y*pixel_count_x+pixel_x] = Mandelbrot(x,y,pixel_size);
    }
  }	
	}
	downsample_pixels(pixels,out_pixels,super_res,pixel_count_x,pixel_count_y);
  write_to_p6((char *) "out.ppm",pixel_count_x/super_res,pixel_count_y/super_res, (unsigned char *) out_pixels);
	free(pixels);
	free(out_pixels);
  //unsigned char test_array [] = {10, 10, 10, 20, 20, 20};
  //write_to_p6((char *) "test.ppm",2,1,test_array);

  return 0;

}
