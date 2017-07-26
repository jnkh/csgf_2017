#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
#include <cmath>
using namespace std;


float floatMod(float a, float b)
{
    return (a - b * floor(a / b));
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
	float h;
	float s;
	float v;
};

struct rgb_color {
	unsigned char r;
	unsigned char g;
	unsigned char b;
};

rgb_color hsv_to_rgb(float h, float s, float v) {
	h = 360*h;
	
	float c = v*s;
	float hp = h/60;
	float x = c*(1 - abs(floatMod(hp,2.0)-1));
	float r,g,b;
	rgb_color rgb;
	if (0 <= hp && hp < 1) {r = c; g = x; b = 0;}
	if (1 <= hp && hp < 2) {r = x; g = c; b = 0;}
	if (2 <= hp && hp < 3) {r = 0; g = c; b = x;}
	if (3 <= hp && hp < 4) {r = 0; g = x; b = c;}
	if (4 <= hp && hp < 5) {r = x; g = 0; b = c;}
	if (5 <= hp && hp <= 6) {r = c; g = 0; b = x;}
	float m = v-c;
	r += m;
	g += m;
	b += m;
	rgb.r = (unsigned char) round(255*r);
	rgb.g = (unsigned char) round(255*g);
	rgb.b = (unsigned char) round(255*b);
	return rgb;
}

rgb_color get_color(float distance,int iter,int iter_max, float pixel_size, float radius, float radius_max) {
	rgb_color ret;
	float hue,val;
	float sat = 0.7;
	if (iter >= iter_max) {
		ret.r = 255; ret.g = 255; ret.b = 255	;
		return ret;
	}
	if (distance < 0.5*pixel_size) {
		val = pow(distance/(0.5*pixel_size),1.0/3);
	}	else {
		val = 1.0;
	}	
	float iter_cont = iter - log2(log(radius)/log(radius_max));
	hue = log(iter_cont)/log(iter_max);
	hue = hue*10;
	hue = hue - floor(hue);
	return hsv_to_rgb(hue,sat,val);
}

float abs_complex(float x,float y) {return sqrt(x*x + y*y);}

float mult_complex_re(float x1,float y1,float x2, float y2) { return x1*x2 - y1*y2;}
float mult_complex_im(float x1,float y1,float x2, float y2) { return x1*y2 + y1*x2;}

float square_complex_re(float x,float y) { return x*x - y*y;}
float square_complex_im(float x,float y) { return 2*x*y;}


rgb_color Mandelbrot(float x0, float y0,float pixel_size) {

  int iter = 0;
  int iter_max = 10000;
  float radius_max = 1 << 18;

  float radius = 0.0;
  float x = 0.0;
  float y = 0.0;
  float dx = 0.0;
  float dy = 0.0;

  
  while (radius < radius_max && iter < iter_max) {

		dx = 2*mult_complex_re(x,y,dx,dy) + 1;
		dy = 2*mult_complex_im(x,y,dx,dy);

		x = square_complex_re(x,y) + x0;
		y = square_complex_im(x,y) + y0;


		radius = abs_complex(x,y);
	  iter++;
  }
	
	float distance = 2*log(radius)*radius/abs_complex(dx,dy);
	return get_color(distance,iter,iter_max,pixel_size,radius,radius_max);
}


int main () {

  int pixel_count_x = 40;

  float center_x = -0.75;
  float center_y = 0.00;
  float length_x = 2.75;
  float length_y = 2.0;


  float pixel_size = length_x / pixel_count_x;
	int pixel_count_y = floor(length_y/pixel_size);
 	length_y = pixel_count_y*pixel_size;

  float minx = center_x - length_x/2.0;
  float maxy = center_y + length_y/2.0; 
  

  rgb_color * pixels = (rgb_color *) malloc( sizeof(rgb_color)*pixel_count_x*pixel_count_y );

  for (int pixel_y=0; pixel_y<pixel_count_y; pixel_y++) {
		#pragma acc parallel
    for (int pixel_x=0; pixel_x<pixel_count_x; pixel_x++) {

      float x = minx + pixel_x*pixel_size;
      float y = maxy - pixel_y*pixel_size;

			//float x = minx + i*pixel_size;
      //float y = maxy - j*pixel_size;

      pixels[pixel_y*pixel_count_x+pixel_x] = Mandelbrot(x,y,pixel_size);
    }
  }	

  write_to_p6((char *) "out.ppm",pixel_count_x,pixel_count_y, (unsigned char *) pixels);
	free(pixels);
  //unsigned char test_array [] = {10, 10, 10, 20, 20, 20};
  //write_to_p6((char *) "test.ppm",2,1,test_array);
  return 0;

}
