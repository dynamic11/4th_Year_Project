#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cuComplex.h>
#include "kernel.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <time.h>
#include "device_launch_parameters.h"

const int SIZE = 1024;
const double MINFREQ = 1e-6;
using namespace std;
bool debug;

struct dataInfo {
	double freqHigh;
	double freqLow;
	int NFreq;
	int NPorts;
};

__global__ void VectorAdd(double *freq, cuComplex *Ahat, cuComplex *data, cuComplex *Poles, dataInfo *frequencyInfo, int *Apattern, int Ahat_size, int NComplexPoles, int NRealPoles)
{
	double real=0;
	double imag=0;
	double denum=0;
	int poleNumb = 0;
	int NRow = (*frequencyInfo).NFreq;
	double s;
	int test = 0;
		int col = blockIdx.x;
		int row = threadIdx.x;
		if (col > 17) {
			poleNumb = col - 17;
		}
		else {
			poleNumb = col;
		}

		real = Poles[poleNumb].x;
		imag = Poles[poleNumb].y;

			s = 2 * M_PI*freq[row];

		//real		
		if (Apattern[col] == 1) {
			Ahat[col*NRow + row].x = -real / (pow(real, 2) + pow(s, 2));
			Ahat[col*NRow + row].y = -s / (pow(real, 2) + pow(s, 2));
		}
		else if (Apattern[col] == 2) {
			denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
			Ahat[col*NRow + row].x = -2 * (real*(pow(real, 2) + pow(s, 2) + pow(imag, 2))) / denum;
			Ahat[col*NRow + row].y = -2 * (s *(pow(real, 2) + pow(s, 2) - pow(imag, 2))) / denum;
		}
		else if (Apattern[col] == 3) {
			denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
			Ahat[(col)*NRow + row].x = (-2 * imag*(pow(real, 2) - pow(s, 2) + pow(imag, 2))) / denum;
			Ahat[(col)*NRow + row].y = (-4 * real*imag*s) / denum;
		}
		else if (Apattern[col] == 4) {
			Ahat[col*NRow + row].x = 1;
			Ahat[col*NRow + row].y = 0;
		}
		else if (Apattern[col] == -1) {
			denum = pow(real, 2) + pow(imag, 2) - 2 * imag*s + pow(s, 2);
			Ahat[col*NRow + row].x = (real*data[row].x - data[row].y*s + data[row].y*imag) / denum;
			Ahat[col*NRow + row].y = (s*data[row].x - data[row].x*imag + data[row].y*real) / denum;
		}
		else if (Apattern[col] == -2) {
			denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
			Ahat[col*NRow + row].x = 2 * (pow(real, 3)*data[row].x - pow(real, 2)*data[row].y*s + real* data[row].x*pow(s, 2) + real* data[row].x*pow(imag, 2) + s*data[row].y*pow(imag, 2) - data[row].y*pow(s, 3)) / denum;
			Ahat[col*NRow + row].y = 2 * (pow(real, 3)*data[row].y + pow(real, 2)*data[row].x*s + real* data[row].y*pow(s, 2) + real* data[row].y*pow(imag, 2) - s*data[row].x*pow(imag, 2) + data[row].x*pow(s, 3)) / denum;
		}
		else if (Apattern[col] == -3) {
			denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
			Ahat[col*NRow + row].x = 2 * imag*(pow(real, 2)*data[row].x - data[row].x* pow(s, 2) + data[row].x*pow(imag, 2) - 2 * real*data[row].y*s) / denum;
			Ahat[col*NRow + row].y = 2 * imag*(2 * data[row].x*real*s + pow(real, 2)*data[row].y - data[row].y*pow(s, 2) + data[row].y*pow(imag, 2)) / denum;
		}

}



//__global__ void AMatrixFormulationGPU(double *freq, cuComplex *Ahat, cuComplex *data, cuComplex *Poles, dataInfo *frequencyInfo, int *Apattern, int Ahat_size, int NComplexPoles, int NRealPoles) {
//	for (int blocky = 0; blocky <= NColToTake; blocky++) {
//		for (int blockx = 0; blockx < NColToTake; blockx++) {
//			if (blocky == NColToTake) {
//				for (int col = 0; col < baseMatrix_NCol; col++) {
//					for (int row = 0; row < baseMatrix_NRow; row++) {
//
//						realBase = baseMatrix[col*NRow + row].x;
//						imagBase = baseMatrix[col*NRow + row].y;
//
//						realData = data[row + NFreq*blockx].x;
//						imagData = data[row + NFreq*blockx].y;
//
//						int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColToTake)+col*NFreq*NColToTake + blockx*NFreq + row;
//
//						Ahat[memPosition].x = (realBase*realData - imagBase*imagData);
//						Ahat[memPosition].y = (realBase*imagData - imagBase*realData);
//					}
//				}//end for col
//			}
//			else if (blocky == blockx) {
//				for (int col = 0; col < baseMatrix_NCol; col++) {
//					for (int row = 0; row < baseMatrix_NRow; row++) {
//						int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColToTake)+col*NFreq*NColToTake + blockx*NFreq + row;
//						Ahat[memPosition].x = baseMatrix[col*NRow + row].x;
//						Ahat[memPosition].y = baseMatrix[col*NRow + row].y;
//					}//end for row
//				}//end for col
//			}
//			else {
//				for (int col = 0; col < baseMatrix_NCol; col++) {
//
//					for (int row = 0; row < baseMatrix_NCol; row++) {
//
//						int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColToTake)+col*NFreq*NColToTake + blockx*NFreq + row;
//						Ahat[memPosition].x = 0;
//						Ahat[memPosition].y = 0;
//					}//end for row
//				}//end for col
//			}//endif
//		}
//	}
//
//}


//
void readFile(string fileName, double *freq, cuComplex *data, dataInfo *dataInfo) {
	std::string::size_type sz;

	ifstream infile(fileName);

	//Ali: Attempt to open data file 
	if (!infile) {
		std::cout << "While opening data file an error was encountered" << std::endl;
	}
	else {
		string line;
		int fileColumn, freqCount = 0;
		bool skipline;
		string word;
		int currentFileDataCol = 0, currentPole = 0, startofPoleCol = 0, endofPoleCol = 0, dataType = 0;

		//Ali: iterate through each line
		while (getline(infile, line)) {
			istringstream stringOfLine(line);
			fileColumn = 0;
			int dataColCount = 0;

			//Ali: Iterrate through each element of line which is refered to as a column
			while (stringOfLine) {

				//Ali: escape after last pole
				if (fileColumn == 2 * pow((*dataInfo).NPorts, 2) + 1) {
					break;
				}

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
						(*dataInfo).freqHigh = freq[freqCount];
						(*dataInfo).freqLow = freq[freqCount];
					}
					else {
						//Ali: if current freq is greater than highest then update 
						if ((*dataInfo).freqHigh < freq[freqCount]) {
							(*dataInfo).freqHigh = freq[freqCount];
						}
						//Ali: if current freq is smaller than lowest then update 
						if ((*dataInfo).freqLow > freq[freqCount]) {
							(*dataInfo).freqLow = freq[freqCount];
						}
					}
					fileColumn++;

				}
				else {
					//Ali: stores the data coloumn number without freq coloumn
					currentFileDataCol = fileColumn - 1;
					//Ali: calculate the current port we are reading
					currentPole = (int)(currentFileDataCol / ((*dataInfo).NPorts * 2));
					//Ali: calculate the start coloumn and end coloumn to start and stop storing data for port
					startofPoleCol = currentPole * 2 * (*dataInfo).NPorts;
					endofPoleCol = startofPoleCol + (currentPole + 1) * 2 - 1;

					//Ali: check if it should store the current coloumn
					if (currentFileDataCol >= startofPoleCol && currentFileDataCol <= endofPoleCol) {
						//Ali: The second column is the real part of the response
						if (dataType == 0) {
							data[freqCount + (*dataInfo).NFreq*(dataColCount)].x = stod(word, &sz);
						}//endif
						//Ali: The third column is the imag part of the response
						if (dataType == 1) {
							data[freqCount + (*dataInfo).NFreq*(dataColCount)].y = stod(word, &sz);
							dataColCount++;
						}//endif
						dataType ^= 1;
					} //endif
					fileColumn++;
				}//endif
			}//endwhile

			//Ali: If line is skipped then dont add to freq count
			if (!skipline) {
				freqCount++;
			}
		}

		//Ali: make sure the min is atleast 1e-6
		if ((*dataInfo).freqLow < MINFREQ) {
			(*dataInfo).freqLow = MINFREQ;
		}
		printf("freqcount %d\n", (*dataInfo).NFreq);
	}
}//enfunction

int main()
{
	debug = true;
	string dataFileName = "radial_stub^S.txt";
	int NRealPoles = 1;
	int NComplexPoles = 2;
	int NPorts = 2;
	int NFreq = 1001;
	int NColOfData = 0;

	//###########################Reading File########################################

	//Ali: Find out how many col we wil have to store based on number of ports
	for (int i = 1; i <= NPorts; i++) {
		NColOfData += i;
	}

	//Ali: var to store freq points in data (upto 1024 data points)
	double *freq;
	cudaMallocManaged(&freq, NFreq * sizeof(double));

	//Ali: store collected data in complex form (upto 1024 data points)
	cuComplex *data;
	cudaMallocManaged(&data, NColOfData * NFreq * sizeof(cuComplex));

	//Ali: store info about the stroed data
	//		-lowest freq
	//		-highest freq
	//		-number of freq points
	dataInfo *dataInfo;
	cudaMallocManaged(&dataInfo, sizeof(dataInfo));

	(*dataInfo).NFreq = NFreq;
	(*dataInfo).NPorts = NPorts;

	//Ali: extract data form file
	readFile(dataFileName, freq, data, dataInfo);


	if (debug) {

		FILE * fp;
		fp = fopen("1_extractedData.txt", "w+");

		fprintf(fp, "********************************************************\n");
		fprintf(fp, "extracted data\n");
		fprintf(fp, "********************************************************\n");
		for (int i = 0; i < NColOfData; i++) {
			fprintf(fp, "\n********************************************************\n");
			fprintf(fp, "col: %d \n", i);
			for (int z = 0; z < (*dataInfo).NFreq; z++) {
				fprintf(fp, "Z: %d FREQ: %f %f(%f) \n", z, freq[z], data[i*NFreq + z].x, data[i*NFreq + z].y);
			}
		}
		fclose(fp);
		printf("^^^^^^^^^^^^^^^^^^^\n");
		printf("HiegestFREQ: %f GHz \n LowestFREQ: %f GHz \n FreqPoints %d\n", (*dataInfo).freqHigh, (*dataInfo).freqLow, (*dataInfo).NFreq);
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
	cudaMallocManaged(&Poles, (NRealPoles + NComplexPoles) * sizeof(cuComplex));

	int B[2] = { 1, 1 };
	int C[2] = { 1, -1 };

	double poleSpacing = ((*dataInfo).freqHigh - (*dataInfo).freqLow) / (NumberOfPoles - 1);

	for (int z = 0; z < NumberOfPoles; z++) {
		Poles_imag_part[z] = (*dataInfo).freqLow + poleSpacing*z;
		Poles_real_part[z] = -Poles_imag_part[z] / Real_part_Divisor;
	}

	//Set Real Poles
	for (int z = 0; z < NRealPoles; z++) {
		Real_Poles[z] = 2 * M_PI*Poles_real_part[z];
	}

	//Set Complex Poles
	int poleIndex = 0;
	for (int z = 0; z < NComplexPoles / 2; z++) {
		for (int i = 0; i < 2; i++) {
			Complex_Poles[poleIndex].x = 2 * M_PI*(Poles_real_part[NRealPoles + z] * B[i]);
			Complex_Poles[poleIndex].y = 2 * M_PI*(Poles_imag_part[NRealPoles + z] * C[i]);
			poleIndex++;
		}
	}

	//merge Poles into one matrix
	for (int z = 0; z < (NComplexPoles + NRealPoles); z++) {
		if (z < NRealPoles) {
			Poles[z].x = Real_Poles[z];
		}
		else {
			Poles[z].x = Complex_Poles[z - NRealPoles].x;
			Poles[z].y = Complex_Poles[z - NRealPoles].y;
		}
	}

	printf("\n\n********************************************************\n");
	printf("Initial Poles\n");
	printf("********************************************************\n");
	for (int z = 0; z < NRealPoles; z++) {
		printf("Real Pole[%d]: %f(%f) \n", z, Real_Poles[z]);
	}
	for (int z = 0; z < NComplexPoles; z++) {
		printf("Complex Pole[%d]: %f(%f) \n", z, Complex_Poles[z].x, Complex_Poles[z].y);
	}
	printf("^^^^^^^^^^^^^\n");
	for (int z = 0; z < NComplexPoles + NRealPoles; z++) {
		printf("merged Pole[%d]: %f(%f) \n", z, Poles[z].x, Poles[z].y);
	}
	printf("^^^^^^^^^^^^^\n");
	printf("NumberOfPoles: %d  \n poleSpacing: %f \n", NumberOfPoles, poleSpacing);





	//########################### Base Matrix Setup ########################################

	// NCol is equal to the number of real poles + number of imaginary poles + 1 (d col)
	int baseMatrix_NCol = NComplexPoles + NRealPoles + 1;
	int baseMatrix_NRow = (*dataInfo).NFreq;

	int NRow = (*dataInfo).NFreq;

	cuComplex *baseMatrix;
	cudaMallocManaged(&baseMatrix, baseMatrix_NRow * baseMatrix_NCol * sizeof(double));

	int CPUenable = 1;
	if (CPUenable == 1 || CPUenable == 2) {
		clock_t tStart = clock();
		/* Do your stuff here */

		double real = 0, imag = 0, denum = 0, s=0;
		int poleNumb = 0, isReal = 1;
	
		//for loop to generate base matrix
		for (int col = 0; col < baseMatrix_NCol; col++) {

			if (poleNumb < NRealPoles + NComplexPoles) {
				real = Poles[poleNumb].x;
				imag = Poles[poleNumb].y;
			}

			for (int row = 0; row < baseMatrix_NRow; row++) {
				s = 2 * M_PI*freq[row];
				//real pole		
				if (col < NRealPoles) {
					baseMatrix[col*NRow + row].x = -real / (pow(real, 2) + pow(s, 2));
					baseMatrix[col*NRow + row].y = -s / (pow(real, 2) + pow(s, 2));

				}
				//imag pole	real part
				else if (col < NRealPoles + NComplexPoles && isReal) {
					denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
					baseMatrix[col*NRow + row].x = -2 * (real*(pow(real, 2) + pow(s, 2) + pow(imag, 2))) / denum;
					baseMatrix[col*NRow + row].y = -2 * (s *(pow(real, 2) + pow(s, 2) - pow(imag, 2))) / denum;
				}
				//imag pole	imag part
				else if (col < NRealPoles + NComplexPoles && !isReal) {
					denum = (pow(real, 2)*(pow(real, 2) + 2 * pow(s, 2) + 2 * pow(imag, 2)) + pow(imag, 4) - 2 * pow(imag, 2)*pow(s, 2) + pow(s, 4));
					baseMatrix[col*NRow + row].x = (-2 * imag*(pow(real, 2) - pow(s, 2) + pow(imag, 2))) / denum;
					baseMatrix[col*NRow + row].y = (-4 * real*imag*s) / denum;
				}
				//d col
				else if (col == NRealPoles + NComplexPoles) {
					baseMatrix[col*NRow + row].x = 1;
					baseMatrix[col*NRow + row].y = 0;
				}//endif
			};//end row for loop
			
			poleNumb++;

			if(col>= NRealPoles && col < NRealPoles + NComplexPoles)
				isReal ^= 1;
		}; //end col for loop

		//Write base b matrix to file
		FILE * fp;
		fp = fopen("2_Basematrix.txt", "w+");
		for (int row = 0; row < baseMatrix_NRow; row++) {
			for (int col = 0; col < baseMatrix_NCol; col++) {
				fprintf(fp, " %.4e(%.4e)", baseMatrix[col*NRow + row].x, baseMatrix[col*NRow + row].y);
			};
			fprintf(fp, "\n");
		};
		fclose(fp);


		double realBase, imagBase, realData, imagData;
		int dim = NColOfData;


		//###########################Ahat set up########################################
		//Ali: generate the AhatMatrix
		cuComplex *Ahat;
		cudaMallocManaged(&Ahat, (dim)*NRow*(baseMatrix_NCol)*(dim+1) * sizeof(double));
		poleNumb = 0;


		for (int blocky = 0; blocky <= NColOfData; blocky++) {
			for (int blockx = 0; blockx < NColOfData; blockx++) {
				if (blocky == NColOfData) {
					for (int col = 0; col < baseMatrix_NCol; col++) {
						for (int row = 0; row < baseMatrix_NRow; row++) {

							realBase = baseMatrix[col*NRow + row].x;
							imagBase = baseMatrix[col*NRow + row].y;

							realData = data[row + NFreq*blockx].x;
							imagData = data[row + NFreq*blockx].y;

							int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColOfData)+col*NFreq*NColOfData + blockx*NFreq + row;

							Ahat[memPosition].x = (realBase*realData - imagBase*imagData);
							Ahat[memPosition].y = (realBase*imagData - imagBase*realData);
						}
					}//end for col
				}
				else if (blocky == blockx) {
					for (int col = 0; col < baseMatrix_NCol; col++) {
						for (int row = 0; row < baseMatrix_NRow; row++) {
							int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColOfData)+col*NFreq*NColOfData + blockx*NFreq + row;
							Ahat[memPosition].x = baseMatrix[col*NRow + row].x;
							Ahat[memPosition].y = baseMatrix[col*NRow + row].y;
						}//end for row
					}//end for col
				}
				else {
					for (int col = 0; col < baseMatrix_NCol; col++) {

						for (int row = 0; row < baseMatrix_NCol; row++) {
							
							int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColOfData)+col*NFreq*NColOfData + blockx*NFreq + row;
							Ahat[memPosition].x = 0;
							Ahat[memPosition].y = 0;
						}//end for row
					}//end for col
				}//endif
			}
		}


		
		clock_t tStop = clock();
		printf("CPU Time taken: %.6fs\n", (double)(tStop - tStart) / CLOCKS_PER_SEC);



	//clock_t start = clock();
	//VectorAdd <<<Ahat_size,(*dataInfo).NFreq  >>> (freq, Ahat, data, Poles, dataInfo, Apattern, Ahat_size, NComplexPoles, NRealPoles);
	//cudaDeviceSynchronize();
	//clock_t stop = clock();
	//printf("GPU Time taken: %.6fs\n", (double)(stop - start) / CLOCKS_PER_SEC);


	//FILE * fp;

	fp = fopen("2_Amatrix.txt", "w+");
	for (int row = 0; row < NColOfData*NFreq; row++) {
		for (int col = 0; col < baseMatrix_NCol*(NColOfData+1); col++) {
			fprintf(fp, " %.4e(%.4e)", Ahat[col*NRow*NColOfData + row].x, Ahat[col*NRow*NColOfData + row].y);
		}
		fprintf(fp, "\n");
	}



	//for (int blocky = 0; blocky <= NColToTake; blocky++) {
	//	for (int blockx = 0; blockx < NColToTake; blockx++) {
	//		for (int row = 0; row < baseMatrix_NRow; row++) {
	//			for (int col = 0; col < baseMatrix_NCol; col++) {
	//				int memPosition = blocky*(baseMatrix_NCol)*(NFreq)*(NColToTake)+col*NFreq*NColToTake + blockx*NFreq + row;
	//				fprintf(fp, " %.4e(%.4e)", Ahat[memPosition].x, Ahat[memPosition].y);
	//			};
	//			fprintf(fp, "\n");
	//		};
	//	}
	//}
	fclose(fp);
}

/*	Poles_imag_part = linspace(f.L, f.H, IP.Nreal + IP.Ncomplex / 2);
	Poles_real_part = -Poles_imag_part / Real_part_Divisor;

	Real_Poles = Poles_real_part(1 : IP.Nreal);

	Complex_Poles = ...
		kron(Poles_real_part(IP.Nreal + 1:end), [1, 1]) + ...
		kron(Poles_imag_part(IP.Nreal + 1:end), [-j, j]);

	initial_Poles = 2 * pi*transpose(cat(2, Real_Poles, Complex_Poles)); */

	

	cudaFree(freq);
	cudaFree(data);
	cudaFree(dataInfo);
	cudaFree(Poles_imag_part);
	cudaFree(Poles_real_part);
	cudaFree(Real_Poles);
	cudaFree(Complex_Poles);

	return 0;
}