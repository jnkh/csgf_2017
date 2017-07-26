#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

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
  ofstream myfile;
  myfile.open ("example.txt");
  myfile << "Writing this to a file.\n";
  myfile.close();

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

  return 0;
}
