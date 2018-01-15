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
const double MINFREQ = 1e-6;
using namespace std;
bool debug;

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
	debug = true;
	string dataFileName = "radial_stub^S.txt";
	int NRealPoles = 1;
	int NComplexPoles = 16;

	//###########################Reading File########################################

	

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

	if (debug) {
		printf("\n\n********************************************************\n");
		printf("extracted data\n");
		printf("********************************************************\n");
		for (int z = 0; z < (*frequencyInfo).fpointCount; z++) {
			printf("Z: %d FREQ: %f RESP: %f(%f) \n", z, freq[z], data[z].x, data[z].y);
		}
		printf("^^^^^^^^^^^^^^^^^^^\n");
		printf("HiegestFREQ: %f GHz \n LowestFREQ: %f GHz \n FreqPoints %d\n", (*frequencyInfo).high, (*frequencyInfo).low, (*frequencyInfo).fpointCount);
	}



	//###########################Initial Pole Guess########################################

	//This constant is predetermined in a paper
	int Real_part_Divisor = 100;

	int NumberOfPoles = NRealPoles + (NComplexPoles / 2);

	double *Poles_imag_part;
	cudaMallocManaged(&Poles_imag_part, NumberOfPoles * sizeof(double));

	double *Poles_real_part;
	cudaMallocManaged(&Poles_real_part, NumberOfPoles * sizeof(double));

	double *Real_Poles;
	cudaMallocManaged(&Real_Poles, NRealPoles * sizeof(double));

	cuComplex *Complex_Poles;
	cudaMallocManaged(&Complex_Poles, NComplexPoles * sizeof(cuComplex));

	cuComplex *Poles;
	cudaMallocManaged(&Poles, (NRealPoles+ NComplexPoles) * sizeof(cuComplex));

	int B[2] = { 1, 1 };
	int C[2] = { 1, -1 };

	double poleSpacing = ((*frequencyInfo).high - (*frequencyInfo).low) / (NumberOfPoles-1);

	for (int z = 0; z < NumberOfPoles; z++) {
		Poles_imag_part[z] = (*frequencyInfo).low + poleSpacing*z;
		Poles_real_part[z] = -Poles_imag_part[z] / Real_part_Divisor;
	}

	//Set Real Poles
	for (int z = 0; z < NRealPoles; z++) {
		Real_Poles[z] = 2*M_PI*Poles_real_part[z];
	}
	
	//Set Complex Poles
	int poleIndex = 0;
	for (int z = 0; z < NComplexPoles/2; z++) {
		for (int i = 0; i< 2; i++) {
			Complex_Poles[poleIndex].x = 2*M_PI*(Poles_real_part[NRealPoles + z] * B[i]);
			Complex_Poles[poleIndex].y = 2*M_PI*(Poles_imag_part[NRealPoles + z] * C[i]);
			poleIndex++;
		}
	}

	//merge Poles into one matrix
	printf("444444\n");
	for (int z = 0; z < (NComplexPoles+ NRealPoles); z++) {
		if (z < NRealPoles) {
			Poles[z].x = Real_Poles[z];
		}else {
			Poles[z].x = Complex_Poles[z - NRealPoles].x;
			Poles[z].y = Complex_Poles[z - NRealPoles].y;
		}
	}

	printf("444444\n");
	printf("\n\n********************************************************\n");
	printf("Initial Poles\n");
	printf("********************************************************\n");
	for (int z = 0; z < NRealPoles; z++) {
		printf("Real Pole[%d]: %f(%f) \n", z, Real_Poles[z]);
	}
	for (int z = 0; z <  NComplexPoles; z++) {
		printf("Complex Pole[%d]: %f(%f) \n", z, Complex_Poles[z].x, Complex_Poles[z].y);
	}
	printf("^^^^^^^^^^^^^\n");
	for (int z = 0; z < NComplexPoles + NRealPoles; z++) {
		printf("merged Pole[%d]: %f(%f) \n", z, Poles[z].x, Poles[z].y);
	}
	printf("^^^^^^^^^^^^^\n");
	printf("NumberOfPoles: %d  \n poleSpacing: %f \n", NumberOfPoles, poleSpacing);

	//###########################Ahat set up########################################
	cuComplex *Ahat;
	//size_t pitch;
	//cudaMallocPitch(&devPtr, &devPitch, Ncols * sizeof(float), Nrows);
	//cudaMallocPitch(&Ahat, &pitch, Ncols * sizeof(float), Nrows));
	int NPorts = 1;
	int* Apattern;
	cudaMallocManaged(&Apattern, (NRealPoles + NComplexPoles+ NPorts) * sizeof(int));

	int isReal = 1;
	for (int i = 0; i < NRealPoles*2 + NComplexPoles * 2 + NPorts; i++) {
		if (i < NRealPoles) {
			Apattern[i]=1;	
		}else if (i < NRealPoles+NComplexPoles) {
			Apattern[i] = (isReal) ? 2 : 3;
			isReal ^= 1;	
		}else if (i < NRealPoles + NComplexPoles+NPorts) {
			Apattern[i] = 4;
		}else if ((i- NComplexPoles - NPorts- NRealPoles) < NRealPoles) {
			Apattern[i] = -1;
		}
		else if ((i - NComplexPoles - NPorts - NRealPoles)  < NRealPoles + NComplexPoles) {
			Apattern[i] = (isReal) ? -2 : -3;
			isReal ^= 1;
		}
	};



	printf(" A Pattern: ");
	for (int i = 0; i < NRealPoles*2 + NComplexPoles * 2 + NPorts; i++) {
			printf("%d ", Apattern[i]);
	};
	printf("\n");

	int NCol = NRealPoles * 2 + NComplexPoles * 2;
	int NRow = (*frequencyInfo).fpointCount;
	cudaMallocManaged(&Ahat, NRow*((2*NRealPoles + 2*NComplexPoles+NPorts) * sizeof(double)));
	double s;
	printf("here 1\n");
	int g = 0;
/*	for (int col = 0; col < NRealPoles; col++) {
		for (int row = 0; row < (*frequencyInfo).fpointCount; row++) {
			s = 2 * M_PI*freq[row];
			//Ahat[col*NRow + row].x = (1 / (freq[row] - Real_Poles[col]);
			Ahat[col*NRow + row].x = -Real_Poles[col]/(pow(Real_Poles[col],2) + pow(s,2));
			Ahat[col*NRow + row].y = -s / (pow(Real_Poles[col], 2) + pow(s, 2));
			g++;
		};
	}; */

	printf("here 2\n");
	double real;
	double imag;
	double denum;
	int poleNumb = 0;
	int test = 0;
	for (int col = 0; col < NComplexPoles*2+NRealPoles*2+NPorts; col++) {
		real = Poles[poleNumb].x;
		imag = Poles[poleNumb].y;

		printf("col: %d poleNumb: %d \n", col, poleNumb);
		for (int row = 0; row < (*frequencyInfo).fpointCount; row++) {
			s = 2 * M_PI*freq[row];
			
			//real		
			if(Apattern[col] == 1){
				Ahat[col*NRow + row].x = -real / (pow(real, 2) + pow(s, 2));
				Ahat[col*NRow + row].y = -s / (pow(real, 2) + pow(s, 2));
			}else if (Apattern[col] == 2) {
				denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
				Ahat[col*NRow + row].x = -2 * (real*(pow(real, 2) + pow(s, 2) + pow(imag, 2))) / denum;
				Ahat[col*NRow + row].y = -2 * (s *(pow(real, 2) + pow(s, 2) - pow(imag, 2))) / denum;
			}else if (Apattern[col] == 3) {
				denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
				Ahat[(col)*NRow + row].x = (-2 * imag*(pow(real, 2) - pow(s, 2) + pow(imag, 2))) / denum;
				Ahat[(col)*NRow + row].y = (-4 * real*imag*s) / denum;
			}else if (Apattern[col] == 4) {
				Ahat[col*NRow + row].x = 1;
				Ahat[col*NRow + row].y = 0;
			}
			else if (Apattern[col] == -1) {
				denum = pow(real, 2) + pow(imag, 2) - 2 * imag*s + pow(s, 2);
				Ahat[col*NRow + row].x = (real*data[row].x- data[row].y*s+ data[row].y*imag) / denum;
				Ahat[col*NRow + row].y = (s*data[row].x - data[row].x*imag + data[row].y*real) / denum;
			}
			else if (Apattern[col] == -2) {
				denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
				Ahat[col*NRow + row].x = 2 * (pow(real, 3)*data[row].x - pow(real, 2)*data[row].y*s + real* data[row].x*pow(s, 2) + real* data[row].x*pow(imag, 2) + s*data[row].y*pow(imag, 2) - data[row].y*pow(s, 3)) / denum;
				Ahat[col*NRow + row].y = 2 * (pow(real, 3)*data[row].y + pow(real, 2)*data[row].x*s + real* data[row].y*pow(s, 2) + real* data[row].y*pow(imag, 2) - s*data[row].x*pow(imag, 2) + data[row].x*pow(s, 3)) / denum;
			}
			else if (Apattern[col] == -3) {
				denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
				Ahat[col*NRow + row].x = 2*imag*(pow(real, 2)*data[row].x - data[row].x* pow(s, 2) + data[row].x*pow(imag, 2) -2*real*data[row].y*s) / denum;
				Ahat[col*NRow + row].y = 2 *imag*(2*data[row].x*real*s + pow(real, 2)*data[row].y - data[row].y*pow(s, 2) + data[row].y*pow(imag,2))  / denum;
			}

			//imag			
			g++;		
		};
		if (poleNumb < (NComplexPoles + NRealPoles)) {
			poleNumb++;
		}
		else {
			poleNumb = 0;
		}
	};

	FILE * fp;

	fp = fopen("file.txt", "w+");
	for (int row = 0; row <(*frequencyInfo).fpointCount; row++) {
		for (int col = 0; col < NComplexPoles*2 + NRealPoles*2 + NPorts; col++) {
			fprintf(fp," %.4e(%.4e)", Ahat[col*NRow + row].x, Ahat[col*NRow + row].y);
		};
		fprintf(fp,"\n");
	};

	fclose(fp);

	return(0);
	


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
	cudaFree(Poles_imag_part);
	cudaFree(Poles_real_part);
	cudaFree(Real_Poles);
	cudaFree(Complex_Poles);

	return 0;
}