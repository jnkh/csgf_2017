#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <math.h>
using namespace std;

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

float abs_complex(float x,float y) {return sqrt(x*x + y*y);}

float mult_complex_re(float x1,float y1,float x2, float y2) { return x1*x2 - y1*y2;}
float mult_complex_im(float x1,float y1,float x2, float y2) { return x1*y2 + y1*x2;}

float square_complex_re(float x,float y) { return x*x - y*y;}
float square_complex_im(float x,float y) { return 2*x*y;}


int Mandelbrot(float x0, float y0) {

  int iter = 0;
  int iter_max = 10000;
  float radius_max = 1 << 18;

  float radius = 0.0;
  float x = 0.0;
  float y = 0.0;
  float dx = 0.0;
  float dy = 0.0;

  
  while (radius < radius_max && iter < iter_max) {
		x = square_complex_re(x,y) + x0;
		y = square_complex_im(x,y) + y0;

		dx = 2*mult_complex_re(x,y,dx,dy) + 1;
		dy = 2*mult_complex_im(x,y,dx,dy);

		radius = abs_complex(x,y);
	  iter++;
  }

  //cout << x0 << " " << y0 << endl;
  //cout << "final iter " << iter << endl; 
  return (float)iter;
}


int main () {

  int pixel_count_x = 100;
  //int pixel_count_x = 8192;

  float center_x = -0.75;
  float center_y = 0.00;
  float length_x = 2.75;
  float length_y = 2.0;

  float pixel_size = length_x / pixel_count_x;

  int pixel_count_y = floor(length_y/pixel_size);
  length_y = pixel_count_y*pixel_size;


  float minx = center_x - length_x/2.0;
  float maxy = center_y + length_y/2.0; 
  

  float *pixels = (float *) malloc( sizeof(float)*pixel_count_x*pixel_count_y );

  for (int i=0; i<pixel_count_x; i++) {
    for (int j=0; j<pixel_count_y; j++) {

      float x = minx + i*pixel_size;
      float y = maxy - j*pixel_size;

      pixels[i+pixel_count_x*j] = Mandelbrot(x,y);

    }
  }	

  ofstream myfile;
  unsigned char test_array [] = {10, 10, 10, 20, 20, 20};
  write_to_p6((char *) "test.ppm",2,1,test_array);
  return 0;

}
