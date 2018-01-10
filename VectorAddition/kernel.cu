#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cuComplex.h>
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

		double *freq;
		cudaMallocManaged(&freq, SIZE * sizeof(double));

		struct response {
			double real;
			double imag;
		};


		response *data;
		cudaMallocManaged(&data, SIZE * sizeof(response));

		struct frequency {
			double high;
			double low;
			int fpointCount;
		};

		frequency *frequencyInfo;
		cudaMallocManaged(&frequencyInfo, sizeof(frequency));

		std::string::size_type sz;
		int skipline = 0;

	ifstream infile("radial_stub^S.txt");
	
	if (!infile)
	{
		std::cout << "While opening data file an error was encountered" << std::endl;
	} else {
		string line;
		int h = 0;
		(*frequencyInfo).fpointCount = 0;
		int w = 0;
		while (getline(infile, line)) {
			istringstream iss(line);
			int i = 1;			
			while (iss) {
				string word;
				iss >> word;
				skipline = 0;
				if (word.compare("%#") == 0) {
					skipline = 1;
					break;
				}
				if (i == 1) {
					freq[h] = stod(word, &sz);
					(*frequencyInfo).fpointCount++;

					if (w == 0) {
						(*frequencyInfo).high = freq[h];
						(*frequencyInfo).low = freq[h];
						w++;
					}
					else {
						if ((*frequencyInfo).high < freq[h]) {
							(*frequencyInfo).high = freq[h];
						}
							

						if ((*frequencyInfo).low > freq[h]) {
							(*frequencyInfo).low = freq[h];
						}	
					}

				}
				if (i == 2) {
					data[h].real = stod(word, &sz);
				}
				if (i == 3) {
					data[h].imag = stod(word, &sz);
				}
				i++;
			}
			if (skipline == 0) {
				h++;
			}		
		}

		if ((*frequencyInfo).low < MINFREQ) {
			(*frequencyInfo).low = MINFREQ;
		}
	}

	for (int z = 0; z < (*frequencyInfo).fpointCount; z++) {
		printf("Z: %d FREQ: %f RESP: %f(%f) \n", z, freq[z], data[z].real, data[z].imag);
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