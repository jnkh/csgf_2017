#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void write_to_p6(char *,int,int,unsigned char []);



int main () {
  	ofstream myfile;
	unsigned char test_array [] = {10, 10, 10, 20, 20, 20};
	write_to_p6((char *) "test.ppm",2,1,test_array);
	return 0;
}

/*takes an array of row major, contiguous rbg values as unsigned char (0-255) and writes them
 * to a ppm file*/
void write_to_p6(char* filename,int dim_x, int dim_y, unsigned char data []) {
	int length = 3*dim_x*dim_y; //colors
	int max_value = 255;
	char buff [20];
  	ofstream f;
  	f.open(filename);
  	f << "P6\n";
	f << dim_x << " " << dim_y << endl; //sprintf(buff,"%d %d\n",dim_x,dim_y);
	f << max_value << endl; 
	f.close();
	f.open(filename,ios::app | ios::binary);
	f.write((char *) data, length);
}

	
