#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cuComplex.h>
#include "kernel.h"
#define _USE_MATH_DEFINES
#include <math.h>

const int SIZE = 1024;
const int MINFREQ = (1 * pow(10, -6));
using namespace std;

//__global__ void VectorAdd(int *a, int *b, int *c, int n)
//{
//	int i = threadIdx.x;
//	if(i < n)
//		c[i] = a[i] + b[i];
//}
//
//void readFile(char *fileName, int *data)
//{
//	int i = threadIdx.x;
//	if (i < n)
//		c[i] = a[i] + b[i];
//}

int main()
{
	string dataFileName = "radial_stub^S.txt";

	//Ali: var to store freq points in data (upto 1024 data points)
	double *freq;
	cudaMallocManaged(&freq, SIZE * sizeof(double));

	//Ali: store collected data in complex form (upto 1024 data points)
	cuComplex *data;
	cudaMallocManaged(&data, SIZE * sizeof(cuComplex));

	//Ali: store info about the stroed data
	//		-lowest freq
	//		-highest freq
	//		-number of freq points
	frequency *frequencyInfo;
	cudaMallocManaged(&frequencyInfo, sizeof(frequency));

	std::string::size_type sz;
	

	ifstream infile(dataFileName);
	
	//Ali: Attempt to open data file 
	if (!infile) {
		std::cout << "While opening data file an error was encountered" << std::endl;
	} else {
		string line;
		int fileColumn, freqCount = 0;
		(*frequencyInfo).fpointCount = 0;
		bool skipline;
		string word;

		//Ali: iterate through each line
		while (getline(infile, line)) {
			istringstream stringOfLine(line);
			fileColumn = 0;

			//Ali: Iterrate through each element of line which is refered to as a column
			while (stringOfLine) {
				//load new element into a string variable called "word"
				stringOfLine >> word;
				skipline = false;

				//Ali: if comment ("%#") is detected then like is skipped
				if (word.compare("%#") == 0) {
					skipline = true;
					break;
				}

				//Ali: if it is the first column then data is stored as the frequency
				if (fileColumn == 0) {
					//Ali: translate string to int 
					freq[freqCount] = stod(word, &sz);

					//Ali: if it is the forst freq reading then initialize max and min
					if (freqCount == 0) {
						(*frequencyInfo).high = freq[freqCount];
						(*frequencyInfo).low = freq[freqCount];
					}else {
						//Ali: if current freq is greater than highest then update 
						if ((*frequencyInfo).high < freq[freqCount]) {
							(*frequencyInfo).high = freq[freqCount];
						}
						//Ali: if current freq is smaller than lowest then update 
						if ((*frequencyInfo).low > freq[freqCount]) {
							(*frequencyInfo).low = freq[freqCount];
						}	
					}

				}
				//Ali: The second column is the real part of the response
				if (fileColumn == 1) {
					data[freqCount].x = stod(word, &sz);
				}
				//Ali: The third column is the imag part of the response
				if (fileColumn == 2) {
					data[freqCount].y = stod(word, &sz);
				}
				fileColumn++;
			}
			//Ali: If line is skipped then dont add to freq count
			if (!skipline) {
				freqCount++;
			}		
		}

		//Ali: Update total freq count
		(*frequencyInfo).fpointCount= freqCount;

		//Ali: make sure the min is atleast 1e-6
		if ((*frequencyInfo).low < MINFREQ) {
			(*frequencyInfo).low = MINFREQ;
		}
	}

	for (int z = 0; z < (*frequencyInfo).fpointCount; z++) {
		printf("Z: %d FREQ: %f RESP: %f(%f) \n", z, freq[z], data[z].x, data[z].y);
	}
	printf("********************************************************\n");
	printf("HiegestFREQ: %f GHz \n LowestFREQ: %f GHz \n FreqPoints %d\n", (*frequencyInfo).high, (*frequencyInfo).low, (*frequencyInfo).fpointCount);

	//##################################################################################

	int Real_part_Divisor = 100;
	int NRealPoles = 1;
	int NComplexPoles = 16;

	int NumberOfPoles = NRealPoles + (NComplexPoles / 2);

	double *Poles_imag_part;
	cudaMallocManaged(&Poles_imag_part, NumberOfPoles * sizeof(double));

	double *Poles_real_part;
	cudaMallocManaged(&Poles_real_part, NumberOfPoles * sizeof(double));

	double poleSpacing=((*frequencyInfo).high - (*frequencyInfo).low) / (NumberOfPoles-1);


	for (int z = 0; z < NumberOfPoles; z++) {
		Poles_imag_part[z] = (*frequencyInfo).low + poleSpacing*z;
		Poles_real_part[z] = Poles_imag_part[z] / Real_part_Divisor;
	}

	double *Real_Poles;
	cudaMallocManaged(&Real_Poles, NRealPoles * sizeof(double));

	for (int z = 0; z < NRealPoles; z++) {
		Real_Poles[z] = Poles_real_part[z];
	}

	cuComplex *Complex_Poles;
	cudaMallocManaged(&Complex_Poles, NComplexPoles * sizeof(cuComplex));

	int B[2] = { 1, 1 };
	int C[2] = { -1, 1 };

	int o = 0;
	for (int z = 0; z < NComplexPoles/2; z++) {
		for (int p = 0; p< 2; p++) {
			Complex_Poles[o].x = 2*M_PI*(Poles_imag_part[NRealPoles + z] * B[p]);
			Complex_Poles[o].y = 2*M_PI*(Poles_imag_part[NRealPoles + z] * C[p]);
			o++;
		}
	}


	for (int z = 0; z < NumberOfPoles; z++) {
		printf("Pole[%d]: %f (%f)  \n", z, Poles_imag_part[z], Poles_real_part[z]);
	}
	printf("********************************************************\n");
	printf("NumberOfPoles: %d GHz \n poleSpacing: %f \n", NumberOfPoles, poleSpacing);

	for (int z = 0; z <  NComplexPoles; z++) {
		printf("Complex Pole[%d]: %f(%f) \n", z, Complex_Poles[z].x, Complex_Poles[z].y);
	}
	printf("********************************************************\n");
	printf("NumberOfPoles: %d GHz \n poleSpacing: %f \n", NumberOfPoles, poleSpacing);


/*	Poles_imag_part = linspace(f.L, f.H, IP.Nreal + IP.Ncomplex / 2);
	Poles_real_part = -Poles_imag_part / Real_part_Divisor;

	Real_Poles = Poles_real_part(1 : IP.Nreal);

	Complex_Poles = ...
		kron(Poles_real_part(IP.Nreal + 1:end), [1, 1]) + ...
		kron(Poles_imag_part(IP.Nreal + 1:end), [-j, j]);

	initial_Poles = 2 * pi*transpose(cat(2, Real_Poles, Complex_Poles)); */

	

	cudaFree(freq);
	cudaFree(data);
	cudaFree(frequencyInfo);

	return 0;
}