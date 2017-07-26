#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
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


int Mandelbrot(float x0, float y0) {

  int iter = 0;
  int iter_max = 1000;
  float radius_max = 2.0;

  float x = 0.0;
  float y = 0.0;

  float radius = 0.0;
  
  while (radius < radius_max && iter < iter_max) {

    float xtemp = x*x - y*y + x0;
    y = 2*x*y + y0;
    x = xtemp;
    radius = x*x + y*y;
    iter++;

  }

  cout << "final iter " << iter << endl; 

  return (float)iter;

}

int main () {

  int pixel_count_x = 8192;

  float center_x = -0.75;
  float center_y = 0.00;
  float length_x = 2.75;
  float length_y = 2.0;

  float pixel_size = length_x / pixel_count_x;
  float pixel_count_y = length_y / pixel_size;

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
