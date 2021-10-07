// Copyright 2016 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <array>

#include <cuda_runtime_api.h>
#include <cuda.h>

#include "SamplePoints.h"
#include "AsciiRaster.h"
#include "Utilities.h"

#include "KDtree.h"
#include "CUDA_KDtree.h"

#include "kde_kernel_kdtr.cu"
//#include "CUDA_KDtree.cu" // it seems this is not needed since we already have CUDA_KDtree.h included

using namespace std;

// distance squared between two points
inline  float Distance2(float x0, float y0, float x1, float y1){
	float dx = x1 - x0;
	float dy = y1 - y0;
	return dx*dx + dy*dy;
}

// mean center of points
void MeanCenter(SamplePoints Points, float &mean_x, float &mean_y);

// (squared) standard distance of points
void StandardDistance2(SamplePoints Points, float &d2);

// bandwidth squared
inline float BandWidth2(SamplePoints Points){
	float d2;
	StandardDistance2(Points, d2);
	return sqrtf(2.0f / (3 * Points.numberOfPoints)) * d2;
}

// Gaussian kernel
inline float GaussianKernel(float h2, float d2){
	if(d2 >= CUT_OFF_FACTOR * h2){
		return 0.0f;
	}
	return expf(d2 / (-2.0f * h2)) / (h2*TWO_PI);
}

//Timothy @ 01/21/2020
//EDIT: Changed AllocateDeviceSamplePoints to return void, and instead utilize pointers to an array
//Changed all functions to utilize the array of pointers
void AllocateDeviceSamplePoints(SamplePoints* dPoints, const SamplePoints Points);
void CopyToDeviceSamplePoints(SamplePoints* dPoints, const SamplePoints hPoints);
void CopyFromDeviceSamplePoints(SamplePoints hPoints, const SamplePoints* dPoints);
SamplePoints AllocateSamplePoints(int n); // random points
SamplePoints ReadSamplePoints(const char *csvFile); // points read from a .csv file
// By Guiming @ 2016-09-04
SamplePoints CopySamplePoints(const SamplePoints Points);
void FreeDeviceSamplePoints(SamplePoints* dPoints);
void FreeSamplePoints(SamplePoints* Points);
void WriteSamplePoints(SamplePoints* Points, const char * csvFile);
void WriteSamplePoints(SamplePoints* Points, float* Hs, float* Ws, const char * csvFile);
void ReformPoints(SamplePoints* dPoints); //Timothy @ 08/13/2021
void DividePoints(SamplePoints* dPoints); //Timothy @ 08/13/2021

void AllocateDeviceAsciiRaster(AsciiRaster* dAscii, const AsciiRaster Ascii);
void CopyToDeviceAsciiRaster(AsciiRaster* dAscii, const AsciiRaster hAscii);
void CopyFromDeviceAsciiRaster(AsciiRaster hAscii, const AsciiRaster dAscii);
AsciiRaster AllocateAsciiRaster(int nCols, int nRows, float xLLCorner, float yLLCorner, float cellSize, float noDataValue);
AsciiRaster ReadAsciiRaster(char * asciiFile); // ascii raster read from a .asc file
AsciiRaster CopyAsciiRaster(const AsciiRaster Ascii);
void FreeDeviceAsciiRaster(AsciiRaster* Ascii);
void FreeAsciiRaster(AsciiRaster* Ascii);
void WriteAsciiRaster(AsciiRaster* Ascii, const char * asciiFile);

float* AllocateEdgeCorrectionWeights(SamplePoints Points); 
void FreeEdgeCorrectionWeights(float* weights);
void ReformECWeights(float** dWeights, float* hWeights); //Timothy @ 08/13/2021

void AllocateDeviceEdgeCorrectionWeights(float** dWeights, SamplePoints Points);
void FreeDeviceEdgeCorrectionWeights(float** weights);

///////// Guiming on 2016-03-16 ///////////////
// the array holding bandwidth at each point
float* AllocateBandwidths(int n); // n is number of points
//Allocation on device now done with pointers instead of return
void AllocateDeviceBandwidths(float** dBandwidths, int n); // n is number of points
void CopyToDeviceBandwidths(float** dBandwidth, const float* hBandwidths, const int n);
void CopyFromDeviceBandwidths(float* hBandwidth, const float* dBandwidths, const int n);
void FreeDeviceBandwidths(float* bandwidths);
void FreeBandwidths(float* bandwidths);
void ReformBandwidths(float** dBand, float* hBand); //Timothy @ 08/13/2021 - Reform bandwidth arrays on host and copy back accross devices

// the array holding inclusive/exclusive density at each point
float* AllocateDen(int n); // n is number of points
void AllocateDeviceDen(float** dDen, int n); // n is number of points
void CopyToDeviceDen(float** dDen, const float* hDen, const int n);
void CopyFromDeviceDen(float* hDen, const float* dDen, const int n);
void CopyDeviceDen(float* dDenTo, const float* dDenFrom, const int n);
void FreeDeviceDen(float** den);
void FreeDen(float* den);
void ReformDensities(float** dDen, float* den); //Timothy @ 12/29/21 - Reforms densities from all devices back into one single array

// compute the optimal Maximum Likelihood Estimation fixed bandwidth
// By Guiming @ 2016-02-26
float MLE_FixedBandWidth(AsciiRaster* Ascii, SamplePoints* Points, float **edgeWeights, float h, float* den0 = NULL, float* den1 = NULL, bool useGPU = false, float** dDen0 = NULL, float** dDen1 = NULL);

// compute fixed bandwidth density at sample points
// By Guiming @ 2016-05-21
void ComputeFixedDensityAtPoints(AsciiRaster Ascii, SamplePoints Points, float *edgeWeights, float h, float* den0 = NULL, float* den1 = NULL, float* dDen0 = NULL, float* dDen1 = NULL);

// compute the log likelihood given single bandwidth h
// By Guiming @ 2016-02-26
float LogLikelihood(AsciiRaster* Ascii, SamplePoints* Points, float **edgeWeights, float h, float* den0 = NULL, float* den1 = NULL, bool useGPU = false, float** dDen0 = NULL, float** dDen1 = NULL);

// compute the log likelihood given bandwidths hs
// By Guiming @ 2016-02-26
// float* den0 : density based on all points, including itself
// float* den1 : leave one out density
float LogLikelihood(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float **edgeWeights, float* hs, float* den0 = NULL, float* den1 = NULL, bool useGPU = false, float** dHs = NULL, float** dDen0 = NULL, float** dDen1 = NULL, float h = 1.0f, float alpha = -0.5f, float** dDen0cpy = NULL);

// compute the log likelihood given a center (h0, alpha0) and step (stepH, stepA)
// By Guiming @ 2016-03-06
void hj_likelihood(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float **edgeWeights, float h0, float alpha0, float stepH, float stepA, int lastdmax, float* logLs, float* hs = NULL, float* den0 = NULL, float* den1 = NULL, bool useGPU = false, float** dHs = NULL, float** dDen0 = NULL, float** dDen1 = NULL, float** dDen0cpy = NULL);

// compute the optimal h and alpha (parameters for calculating the optimal adaptive bandwith)
// By Guiming @ 2016-03-06
void hooke_jeeves(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float **edgeWeights, float h0, float alpha0, float stepH, float stepA, float* optParas, float* hs = NULL, float* den0 = NULL, float* den1 = NULL, bool useGPU = false, float** dHs = NULL, float** dDen0 = NULL, float** dDen1 = NULL, float** dDen0cpy = NULL);

float compGML(float* den0, int n);
///////// Guiming on 2016-03-16 ///////////////


// exact edge effects correction (Diggle 1985)
void EdgeCorrectionWeightsExact(SamplePoints Points, float h, AsciiRaster Ascii, float *weights);
void EdgeCorrectionWeightsExact(SamplePoints Points, float *hs, AsciiRaster Ascii, float *weights);

// check whether the result from sequential computation and that from parallel computation agree
void CheckResults(AsciiRaster AsciiSEQ, AsciiRaster AsciiPARA);

// reduction an array on GPU
void ReductionSumGPU(float* dArray, int numberOfElements);

// extract study area boundary from a raster
// By Guiming @ 2016-09-02
void MarkBoundary(AsciiRaster Ascii, bool useGPU = false);

// compute the closest distances from sample points to study area boundary
// By Guiming @ 2016-09-02
void CalcDist2Boundary(SamplePoints Points, AsciiRaster Ascii, bool useGPU = false);

// sort the sample points on their distances to study area boundary
// By Guiming @ 2016-09-04
void SortSamplePoints(SamplePoints Points);

// comparison function for sort
// By Guiming @ 2016-09-04
int compare ( const void *pa, const void *pb );

void BuildCPUKDtree (SamplePoints Points);
void BuildGPUKDtree ();

void EnableP2P(); //Timothy @ 08/13/2020 - Enable P2P Access Across Devices
void nextDev(int numDev, int& curDev); //Timothy @ 08/14/2020 - Determine next Device to be used
void DevProp(); //Timothy @ 08/24/2020 - Check device properties, primarily for troubleshooting purposes

// By Timothy @ 02/26/2020
//This performs the same tasks as ComputeFixedDensityAtPoints function, however it is designed specifically to run in accross multiple
//GPUs asynchronously
void ComputeFixedDensityDevice(cudaStream_t* streams, AsciiRaster* Ascii, SamplePoints* Points, float** edgeWeights, float h, float* den0, float* den1, float** dDen0, float** dDen1);

/* Run in 2 modes
 *
 * Mode 0: Do not read points and mask from files.
 *         User specify # of points and cell size of the estimated intensity surface.
 *         Random points with x, y coordinates in the range [0,100] will be generated.
 *         The cell size (must be less than 100) determines how many cells in the intensity surface raster.
 *
 *         ./kde_cuda [mode] [#points] [cellsize] [skipSEQ] [skipPARA]
 *         e.g., ./kde_cuda 0 100 1.0 0 0
 *
 * Mode 1: Read points and mask from files.
 *
 *         ./kde_cuda [mode] [points_file] [mask_file] [skipSEQ] [skipPARA]
 *         e.g., ./kde_cuda 1 ../Points.csv ../Mask.asc 0 0
 *
*/

/* be very careful with these global variables
 * they are declared in this way to avoid passing additional parameters in functions
*/
KDtree tree; // pointer to the kd tree, can be accessed in any function
CUDA_KDTree GPU_tree[2]; //pointer to the GPU kd tree, can be accessed in any function. EDIT: A copy of the tree 
//is now kept on each GPU with each of these pointers corresponding to a GPU.

vector <Point> dataP; // pointer to the vector to hold data points in kd tree, initilized when building kd tree
float* gpuDen[2]; // this is a global array allocated on gpu to store density values. Used in DensityAtPointsKdtr
//int* gpu_ret_indexes;
//float* gpu_ret_dists;
//float* zeroDen;
int MAX_N_NBRS = 0;

//Timothy @ 08/13/2020
int GPU_N = 1; //Holds number of GPUs on machine
int GPU_C = 0; //Keeps track of our current GPU

cudaStream_t streams[2]; //Streams to be used for parallelism

SamplePoints sPoints; // sample of point events

int main(int argc, char *argv[]){
	int NPNTS = 100;                // default # of points
	float CELLSIZE = 1.0f;          // default cellsize
	char* pntFn = "data/Points.csv";  // default points file
	char* maskFn = "data/Mask.asc";   // default mask file
	bool fromFiles = true;          // by default, read Points and Mask from files

	int SKIPSEQ = 0;                // by default, do not skip sequential execution
	int SKIPPARA = 0;               // by default, do not skip parallel execution

	//Guiming May 1, 2016
	int Hoption = 0; // 0 for rule of thumb
					 // 1 for h optimal
					 // 2 for h adaptive
	char* denSEQfn = "data/den_SEQ.asc";
	char* denCUDAfn = "data/den_CUDA.asc";

	// parse commandline arguments
	if(argc != 9){
		printf("Incorrect arguments provided. Exiting...\n");
		printf("Run in mode 0:\n ./kde_cuda 0 #points cellsize h_option skip_sequential skip_parallel denfn_seq, denfn_cuda\n");
		printf("Run in mode 1:\n ./kde_cuda 1 points_file mask_file h_option skip_sequential skip_parallel denfn_seq, denfn_cuda\n");
        return 1;
	}
	else{
		int mode = atoi(argv[1]);
		if(mode == 0){
			fromFiles = false;
			NPNTS = atoi(argv[2]);
			CELLSIZE = (float)atof(argv[3]);
			Hoption = atoi(argv[4]);
			SKIPSEQ = atoi(argv[5]);
			SKIPPARA = atoi(argv[6]);
			denSEQfn = argv[7];
			denCUDAfn = argv[8];
		}
		else if(mode == 1){
			pntFn = argv[2];
			maskFn = argv[3];
			Hoption = atoi(argv[4]);
			SKIPSEQ = atoi(argv[5]);
			SKIPPARA = atoi(argv[6]);
			denSEQfn = argv[7];
			denCUDAfn = argv[8];
		}
		else{
			printf("Incorrect arguments provided. Exiting...\n");
			printf("Run in mode 0:\n ./kde_cuda 0 #points cellsize h_option skip_sequential skip_parallel denfn_seq, denfn_cuda\n");
			printf("Run in mode 1:\n ./kde_cuda 1 points_file mask_file h_option skip_sequential skip_parallel denfn_seq, denfn_cuda\n");
	        return 1;
		}

	}

	//Timothy @ 08/13/2020
	//Assign and print number of Compute Capable Devices
	cudaGetDeviceCount(&GPU_N);
	printf("Number of Capable Devices: %d\n", GPU_N);
	printf("Current GPU: %d\n", GPU_C);

	/*for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		DevProp();
	}
	cudaSetDevice(0);*/

	//Timothy @ 08/24/2020
	//Enable P2P Access across devices
	//EnableP2P();

	cudaError_t error;

	//Timothy @ 12/29/2020
	//Create streams for each available device
	
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaStreamCreate(&streams[i]);
	}
	if (error != cudaSuccess)
	{
		printf("Failed to create streams (error code %s)!\n", cudaGetErrorString(error));
		exit(EXIT_FAILURE);
	}
	cudaSetDevice(0); //Reset device to first GPU

	//SamplePoints sPoints; // sample of point events
	AsciiRaster Mask;    // a mask indicating the extent of study area
	AsciiRaster DenSurf, DenSurf_CUDA; // the estimated intensity surface
	float *edgeWeights;  // edge effect correct weights (for each point in the sample)
	bool correction = true; // enable edge effect correction
	srand(100); // If not read from files, generate random points

	//Read or generate points
	if (fromFiles){
		sPoints = ReadSamplePoints(pntFn);
		Mask = ReadAsciiRaster(maskFn);
	}
	else{
		sPoints = AllocateSamplePoints(NPNTS);
		Mask = AllocateAsciiRaster(int(100/CELLSIZE), int(100/CELLSIZE), 0.0f, 0.0f, CELLSIZE, -9999.0f);
	}

	DenSurf = CopyAsciiRaster(Mask);

	// parameters
	int numPoints = sPoints.numberOfPoints;
	int nCols = Mask.nCols;
	int nRows = Mask.nRows;
	float xLLCorner = Mask.xLLCorner;
	float yLLCorner = Mask.yLLCorner;
	float noDataValue = Mask.noDataValue;
	float cellSize = Mask.cellSize;

	printf("number of points: %d\n", numPoints);
	printf("cell size: %f\n", cellSize);
	printf("number of cells: %d\n", nCols * nRows);

	printf("skip executing SEQUENTIAL program? %d\n", SKIPSEQ);
	printf("skip executing PARALLEL program? %d\n", SKIPPARA);
	printf("number of threads per block: %d\n", BLOCK_SIZE);

	// do the work
	float cell_x; // x coord of cell
	float cell_y; // y coord of cell
	float p_x;    // x coord of point
	float p_y;    // x coord of point
	float p_w;    // weight of point
	float e_w = 1.0;    // edge effect correction weight

	float h = sqrtf(BandWidth2(sPoints));
	printf("rule of thumb bandwidth h0: %.5f\n", h);

	// timing
	//double start, stop;
	float elaps_seq, elaps_exc, elaps_inc;

	if(SKIPSEQ == 0){

		edgeWeights = NULL;
		edgeWeights = AllocateEdgeCorrectionWeights(sPoints);

	///////////////////////// SEQUENTIAL /////////////////////////////////

		///////////////////////// START CPU TIMING /////////////////////////////
		cudaEvent_t startCPU;
		error = cudaEventCreate(&startCPU);

		if (error != cudaSuccess)
		{
		   printf("Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		cudaEvent_t stopCPU;
		error = cudaEventCreate(&stopCPU);

		if (error != cudaSuccess)
		{
		   printf("Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(startCPU, NULL);
		if (error != cudaSuccess)
		{
		   printf("Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}
		///////////////////////// END OF START CPU TIMING /////////////////////////////

		// By Guiming @ 2016-09-11
		MarkBoundary(Mask); // either on GPU or CPU
		CalcDist2Boundary(sPoints, Mask);
		//WriteAsciiRaster(&Mask, "output/boundary.asc");
		SortSamplePoints(sPoints);

		// By Guiming @ 2016-11-03
		BuildCPUKDtree(sPoints);

		float* hs = AllocateBandwidths(numPoints);
		for (int i = 0; i < numPoints; i++) 
		{
			hs[i] = h;
		}

	    // compute edge effect correction weights
		EdgeCorrectionWeightsExact(sPoints, h, Mask, edgeWeights);

		if(Hoption == 1){
			float hopt = MLE_FixedBandWidth(&Mask, &sPoints, &edgeWeights, h, NULL, NULL, false);
			printf("cross validated optimal fixed bandwidth hopt: %.5f\n", hopt);

			for(int i = 0; i < numPoints; i++){
				hs[i] = hopt;
			}

			// update edge correction weights
			if(UPDATEWEIGHTS){
				EdgeCorrectionWeightsExact(sPoints, hs, Mask, edgeWeights);
			}
		}

		if(Hoption == 2){
			float* den0 = AllocateDen(numPoints);
			float* den1 = AllocateDen(numPoints);
			float h0 = h;
			float alpha0 = -0.5;
			float stepH = h0/10;
			float stepA = 0.1;
			float* optParas = (float*)malloc(3 * sizeof(float));

			hooke_jeeves(&Mask, &sPoints, NULL, &edgeWeights, h0, alpha0, stepH, stepA, optParas, hs, den0, den1, false);
			h0 = optParas[0];
			alpha0 = optParas[1];
			float logL = optParas[2];

			if(DEBUG) printf("h0: %.5f alpha0: %.5f Lmax: %.5f\n", h0, alpha0, logL);

			free(optParas);
			optParas = NULL;

			ComputeFixedDensityAtPoints(Mask, sPoints, edgeWeights, h0, den0, NULL, false);
			float gml = compGML(den0, numPoints);
			for(int i = 0; i < numPoints; i++){
				hs[i] = h0 * powf(den0[i]/gml, alpha0);
			}
			FreeDen(den0);
			FreeDen(den1);

			// update edge correction weights
			if(UPDATEWEIGHTS){
				EdgeCorrectionWeightsExact(sPoints, hs, Mask, edgeWeights);
			}
		}

		// KDE
		for (int row = 0; row < nRows; row++){
			cell_y = ROW_TO_YCOORD(row, nRows, yLLCorner, cellSize);
			for (int col = 0; col < nCols; col++){
				cell_x = COL_TO_XCOORD(col, xLLCorner, cellSize);
				int idx = row * nCols + col;
				if (DenSurf.elements[idx] != noDataValue){

					float den = 0.0;
					float hp;
					for (int p = 0; p < numPoints; p++){
						p_x = sPoints.xCoordinates[p];
						p_y = sPoints.yCoordinates[p];
						p_w = sPoints.weights[p];
						hp = hs[p];
						if (correction){
							e_w = edgeWeights[p];
						}
						float d2 = Distance2(p_x, p_y, cell_x, cell_y);
						den += GaussianKernel(hp * hp, d2) * p_w *e_w;
					}
					DenSurf.elements[idx] = den; // intensity, not probability
				}
			}
		}

		///////////////////////// STOP CPU TIMING /////////////////////////////
	    // Record the stop event
	    error = cudaEventRecord(stopCPU, NULL);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    // Wait for the stop event to complete
	    error = cudaEventSynchronize(stopCPU);
	    if (error != cudaSuccess)
	    {
	        printf("Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    elaps_seq = 0.0f;
	    error = cudaEventElapsedTime(&elaps_seq, startCPU, stopCPU);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }
	    ///////////////////////// END OF STOP CPU TIMING /////////////////////////////
		printf("MAX_N_NBRS=%d\n", MAX_N_NBRS);
		printf("Computation on CPU took %.3f ms\n\n", elaps_seq);

		// write results to file
		WriteAsciiRaster(&DenSurf, denSEQfn);
		WriteSamplePoints(&sPoints, hs, edgeWeights, "pntsSEQ.csv");

		// clean up (only those not needed any more)
		FreeEdgeCorrectionWeights(edgeWeights);
		//FreeAsciiRaster(&DenSurf);
		FreeBandwidths(hs);
	}
////////////////////////// END OF SEQUENTIAL //////////////////////////////

//////////////////////////  CUDA  /////////////////////////////////////////
	if(SKIPPARA == 0){

		//EDIT:Timothy @ 01/29/2021 
		//Changed the way device variables are initialized, allocated, and copied to utilize arrays which identify which GPU
		//the program should be using.
		DenSurf_CUDA = CopyAsciiRaster(Mask);
		SamplePoints dPoints[2]; 
		float* dWeights[2];
		AsciiRaster dAscii[2];

		AllocateDeviceSamplePoints(dPoints, sPoints);
		AllocateDeviceEdgeCorrectionWeights(dWeights, sPoints);
		AllocateDeviceAsciiRaster(dAscii, Mask);

		// Guiming @ 2016-03-17
		float* hs = AllocateBandwidths(sPoints.numberOfPoints);
		float* zeroDen = AllocateDen(sPoints.numberOfPoints);
		for (int i = 0; i < numPoints; i++) {
			hs[i] = h;
			zeroDen[i] = 0.0f;
		}
		
		float* dHs[2];
		AllocateDeviceBandwidths(dHs, sPoints.numberOfPoints);

		float* den0 = AllocateDen(sPoints.numberOfPoints);
		float* dDen0[2]; 
		AllocateDeviceDen(dDen0, sPoints.numberOfPoints);
		float* dDen0cpy[2]; 
		AllocateDeviceDen(dDen0cpy, sPoints.numberOfPoints);

		float* den1 = AllocateDen(sPoints.numberOfPoints);
		float* dDen1[2];
		AllocateDeviceDen(dDen1, sPoints.numberOfPoints);

		AllocateDeviceDen(gpuDen, sPoints.numberOfPoints);
		//gpu_ret_indexes =
		//gpu_ret_dists =

		printf("Allocate DONE...\n"); //DEBUGGING

		///////////////////////// START GPU INCLUSIVE TIMING /////////////////////////////
		cudaEvent_t startInc;
		error = cudaEventCreate(&startInc);

		if (error != cudaSuccess)
		{
		   printf("Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		cudaEvent_t stopInc;
		error = cudaEventCreate(&stopInc);

		if (error != cudaSuccess)
		{
		   printf("Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(startInc, NULL);
		if (error != cudaSuccess)
		{
		   printf("Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		///////////////////////// END OF START GPU INCLUSIVE TIMING /////////////////////////////
		CopyToDeviceBandwidths(dHs, hs, sPoints.numberOfPoints);
		int pNum = sPoints.numberOfPoints;

		int NBLOCK_W = (pNum + BLOCK_SIZE - 1) / BLOCK_SIZE;
		int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
		dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);

		CopyToDeviceSamplePoints(dPoints, sPoints);
		CopyToDeviceAsciiRaster(dAscii, Mask);

		int cells = dAscii[0].nCols * dAscii[0].nRows;
		CopyToDeviceDen(gpuDen, zeroDen, sPoints.numberOfPoints);

		printf("Copied...\n"); //DEBUGGING

		///////////////////////// START GPU EXCLUSIVE TIMING /////////////////////////////
		cudaEvent_t startExc;
		error = cudaEventCreate(&startExc);

		if (error != cudaSuccess)
		{
		   printf("Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		cudaEvent_t stopExc;
		error = cudaEventCreate(&stopExc);

		if (error != cudaSuccess)
		{
		   printf("Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(startExc, NULL);
		if (error != cudaSuccess)
		{
		   printf("Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}
		///////////////////////// END OF START GPU EXLUSIVE TIMING /////////////////////////////

		
		///////////////////////// START SORTING TIMING /////////////////////////////
		cudaEvent_t startSort;
		error = cudaEventCreate(&startSort);

		if (error != cudaSuccess)
		{
		   printf("Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		cudaEvent_t stopSort;
		error = cudaEventCreate(&stopSort);

		if (error != cudaSuccess)
		{
		   printf("Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(startSort, NULL);
		if (error != cudaSuccess)
		{
		   printf("Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}
		///////////////////////// END OF START SORTING TIMING /////////////////////////////
		///*
		// By Guiming @ 2016-09-11
		MarkBoundary(dAscii[0], true); // either on GPU or CPU
		CalcDist2Boundary(dPoints[0], dAscii[0], true);
		CopyFromDeviceSamplePoints(sPoints, dPoints);
		SortSamplePoints(sPoints);

		//EDIT:Timothy @ 12/10/2020
		//When adding back sorted points, divide points as they are copied accross GPUs
		CopyToDeviceSamplePoints(dPoints, sPoints);

		///////////////////////// STOP SORTING TIMING /////////////////////////////
	    // Record the stop event
	    error = cudaEventRecord(stopSort, NULL);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    // Wait for the stop event to complete
	    error = cudaEventSynchronize(stopSort);
	    if (error != cudaSuccess)
	    {
	        printf("Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    float elaps_sort = 0.0f;
	    error = cudaEventElapsedTime(&elaps_sort, startSort, stopSort);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }
	    ///////////////////////// END OF STOP SORTING TIMING /////////////////////////////
		printf("#Sorting took %.3f ms\n", elaps_sort);
	
		printf("Sorted...\n"); //DEBUGGING

		///////////////////////// START KDTREE TIMING /////////////////////////////
		cudaEvent_t startKd;
		error = cudaEventCreate(&startKd);

		if (error != cudaSuccess)
		{
		   printf("Failed to create start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		cudaEvent_t stopKd;
		error = cudaEventCreate(&stopKd);

		if (error != cudaSuccess)
		{
		   printf("Failed to create stop event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}

		// Record the start event
		error = cudaEventRecord(startKd, NULL);
		if (error != cudaSuccess)
		{
		   printf("Failed to record start event (error code %s)!\n", cudaGetErrorString(error));
		   exit(EXIT_FAILURE);
		}
		///////////////////////// END OF START KDTREE TIMING /////////////////////////////
		
		// By Guiming @ 2016-11-03
		if(SKIPSEQ == 1)
		BuildCPUKDtree(sPoints);
		BuildGPUKDtree(); // needs to build the CPUKDtree first

		///////////////////////// STOP KDTREE TIMING /////////////////////////////
	    // Record the stop event
	    error = cudaEventRecord(stopKd, NULL);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    // Wait for the stop event to complete
	    error = cudaEventSynchronize(stopKd);
	    if (error != cudaSuccess)
	    {
	        printf("Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    float elaps_kd = 0.0f;
	    error = cudaEventElapsedTime(&elaps_kd, startKd, stopKd);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }
	    ///////////////////////// END OF STOP KDTREE TIMING /////////////////////////////
		printf("#Building kd tree took %.3f ms\n", elaps_kd);
		//EDIT: Timothy @ 12/29/2020
		//Run Kernel Asynchronously accross GPUs
		for (int i = 0; i < GPU_N; i++)
		{
			cudaSetDevice(i);
			printf("Current Device: %d\n", i);
			//Alg Step: 1
			CalcEdgeCorrectionWeights<<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>> (h * h, dPoints[i], dAscii[i], dWeights[i]);
			cudaStreamSynchronize(streams[i]);
		}
		cudaSetDevice(0); //Reset device to first GPU
		printf("mGPU Kernel...\n"); //DEBUGGING
		cudaSetDevice(0); //Reset device to first GPU
		// Guiming @ 2016-03-17
		/////////////////////////////////////////////////////////////////////////////////////////
		int numPoints = sPoints.numberOfPoints;
		if(Hoption == 1){
			float hopt = MLE_FixedBandWidth(dAscii, dPoints, dWeights, h, NULL, den1, true, NULL, dDen1);
			printf("cross validated optimal fixed bandwidth hopt: %.5f\n", hopt);

			// kind of combusome
			//Timothy @ 02/19/2020
			//Running following kernels accross all GPUs
			for (int i = 0; i < GPU_N; i++)
			{
				cudaSetDevice(i);
				CalcVaryingBandwidths<<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(dPoints[i], hopt, dHs[i]);
				if (UPDATEWEIGHTS)
				{
					CalcEdgeCorrectionWeights<<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(dHs[i], dPoints[i], dAscii[i], dWeights[i]);
				}
				cudaStreamSynchronize(streams[i]);
			}
			cudaSetDevice(0); //Reset device to first GPU
		}

		if(Hoption == 2){
			float h0 = h;
			float alpha0 = -0.5;
			float stepH = h0/10;
			float stepA = 0.1;
			float* optParas = (float*)malloc(3 * sizeof(float));
			hooke_jeeves(dAscii, dPoints, dPoints, dWeights, h0, alpha0, stepH, stepA, optParas, hs, den0, den1, true, dHs, dDen0, dDen1, dDen0cpy);
			h0 = optParas[0];
			alpha0 = optParas[1];
			float logL = optParas[2];
			if(DEBUG) printf("h0: %.5f alpha0: %.5f Lmax: %.5f\n", h0, alpha0, logL);
			free(optParas);
			optParas = NULL;

			ComputeFixedDensityDevice(streams, dAscii, dPoints, dWeights, h0, NULL, NULL, dDen0, dDen1);
			
			for (int i = 0; i < GPU_N; i++)
			{
				cudaSetDevice(i);
				CopyDeviceDen(dDen0cpy[i], dDen0[i], numPoints);
			}
			cudaSetDevice(0);

			//reform points and densities
			ReductionSumGPU(dDen0cpy[0], numPoints);
			
	    	// update bandwidth on GPU
			//Timothy @ 02/19/2020
			//Running following kernels accross all GPUs
			for (int i = 0; i < GPU_N; i++)
			{
				cudaSetDevice(i);
				CalcVaryingBandwidths<<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(sPoints, dDen0[i], h0, alpha0, dHs[i]);

				// update weights
				//CopyToDeviceBandwidths(dHs, hs, numPoints);
				if (UPDATEWEIGHTS) {

					CalcEdgeCorrectionWeights <<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(dHs[i], dPoints[i], dAscii[i], dWeights[i]);
				}
				cudaStreamSynchronize(streams[i]);
			}
			cudaSetDevice(0); //Reset device to first GPU
		}

		//Reform data
		ReformPoints(dPoints);

		cudaStreamSynchronize(streams[0]);
		printf("Done...\n\n");
		
		/////////////////////////////////////////////////////////////////////////////////

		// invoke kernel to do density estimation
		int NBLOCK_K = (dAscii[0].nCols*dAscii[0].nRows + BLOCK_SIZE - 1) / BLOCK_SIZE;
	    int GRID_SIZE_K = (int)(sqrtf(NBLOCK_K)) + 1;
	    dim3 dimGrid_K(GRID_SIZE_K, GRID_SIZE_K);

		KernelDesityEstimation<<<dimGrid_K, BLOCK_SIZE, 0, streams[0]>>>(dHs[0], dPoints[0], dAscii[0], dWeights[0]);

		///////////////////////// STOP GPU EXCLUSIVE TIMING /////////////////////////////
	    // Record the stop event
	    error = cudaEventRecord(stopExc, NULL);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    // Wait for the stop event to complete
	    error = cudaEventSynchronize(stopExc);
	    if (error != cudaSuccess)
	    {
	        printf("Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    elaps_exc = 0.0f;
	    error = cudaEventElapsedTime(&elaps_exc, startExc, stopExc);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }
	    ///////////////////////// END OF STOP GPU EXCLUSIVE TIMING /////////////////////////////

		// copy results back to host
		CopyFromDeviceAsciiRaster(DenSurf_CUDA, dAscii[0]);

		///////////////////////// STOP GPU INCLUSIVE TIMING /////////////////////////////
	    // Record the stop event
	    error = cudaEventRecord(stopInc, NULL);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to record stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    // Wait for the stop event to complete
	    error = cudaEventSynchronize(stopInc);
	    if (error != cudaSuccess)
	    {
	        printf("Failed to synchronize on the stop event (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }

	    elaps_inc = 0.0f;
	    error = cudaEventElapsedTime(&elaps_inc, startInc, stopInc);

	    if (error != cudaSuccess)
	    {
	        printf("Failed to get time elapsed between events (error code %s)!\n", cudaGetErrorString(error));
	        exit(EXIT_FAILURE);
	    }
	    ///////////////////////// END OF STOP GPU INCLUSIVE TIMING /////////////////////////////
	    printf("Computation on GPU took %.3f ms (EXCLUSIVE)\n", elaps_exc);
	    printf("Computation on GPU took %.3f ms (INCLUSIVE)\n", elaps_inc);

	    if(SKIPSEQ == 0){
			printf("SPEEDUP: %.3fx (EXCLUSIVE) %.3fx (INCLUSIVE)\n", elaps_seq / elaps_exc, elaps_seq / elaps_inc);
			// check resutls
			CheckResults(DenSurf, DenSurf_CUDA);
		}
		// write results to file
		WriteAsciiRaster(&DenSurf_CUDA, denCUDAfn);
		WriteSamplePoints(&sPoints, "pntsCUDA.csv");

		printf("Begin cleanup...\n"); //DEBUGGING
		// clean up
		FreeDeviceSamplePoints(dPoints);
		FreeDeviceEdgeCorrectionWeights(dWeights);
		FreeDeviceAsciiRaster(dAscii);
		FreeSamplePoints(&sPoints);
		// By Guiming @ 2016-09-02
		free(sPoints.distances);
		sPoints.distances = NULL;

		FreeAsciiRaster(&DenSurf);
		FreeAsciiRaster(&DenSurf_CUDA);
		FreeAsciiRaster(&Mask);
		FreeAsciiRaster(dAscii);
		FreeBandwidths(hs);
		FreeDeviceBandwidths(dHs[0]);
		FreeDen(den0);
		FreeDeviceDen(dDen0);
		FreeDeviceDen(dDen0cpy);
		FreeDen(den1);
		FreeDeviceDen(dDen1);
		FreeDen(zeroDen);
		FreeDeviceDen(gpuDen);
	}

	printf("MAX_N_NBRS=%d\n", MAX_N_NBRS);
	printf("Done...\n\n");

	return 0;
}

// mean center of points
void MeanCenter(SamplePoints Points, float &mean_x, float& mean_y){
	float sum_x = 0.0;
	float sum_y = 0.0;

	for (int p = 0; p < Points.numberOfPoints; p++){
		sum_x += Points.xCoordinates[p];
		sum_y += Points.yCoordinates[p];
	}

	mean_x = sum_x / Points.numberOfPoints;
	mean_y = sum_y / Points.numberOfPoints;
}

// standard distance squared
void StandardDistance2(SamplePoints Points, float &d2){

	float mean_x, mean_y;
	MeanCenter(Points, mean_x, mean_y);

	float sum2 = 0.0;

	for (int p = 0; p < Points.numberOfPoints; p++){
		sum2 += Distance2(mean_x, mean_y, Points.xCoordinates[p], Points.yCoordinates[p]);
	}

	d2 = sum2 / Points.numberOfPoints;
}

// generate random sample points
SamplePoints AllocateSamplePoints(int n){
	SamplePoints Points;

	Points.numberOfPoints = n;
	Points.start = 0;
	Points.end = n;
	int size = n*sizeof(float);

	Points.xCoordinates = (float*)malloc(size);
	Points.yCoordinates = (float*)malloc(size);
	Points.weights = (float*)malloc(size);
	Points.distances = (float*)malloc(size); // By Guiming @ 2016-09-02

	for (int i = 0; i < n; i++)
	{
		Points.xCoordinates[i] = rand() * 100.0f / RAND_MAX;
		Points.yCoordinates[i] = rand() * 100.0f / RAND_MAX;
		Points.weights[i] = 1.0f;
		Points.distances[i] = 0.0f; // By Guiming @ 2016-09-02
		//printf("x:%.2f y:%.2f w:%.2f\n", Points.xCoordinates[i], Points.yCoordinates[i], Points.weights[i]);
	}
	return Points;
}

// points read from a .csv file
SamplePoints ReadSamplePoints(const char *csvFile){
	FILE *f = fopen(csvFile, "rt");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	const int CSV_LINE_LENGTH = 256;
	SamplePoints Points;
	int n = 0;
	char line[CSV_LINE_LENGTH];
	char ch;

	while (!feof(f))
	{
		ch = fgetc(f);
		if (ch == '\n')
		{
			n++;
		}
	}

	if (n == 1){
		printf("No point in file!\n");
		exit(1);
	}

	n = n - 1; // do not count the header line
	Points.numberOfPoints = n;
	Points.xCoordinates = (float*)malloc(n*sizeof(float));
	Points.yCoordinates = (float*)malloc(n*sizeof(float));
	Points.weights = (float*)malloc(n*sizeof(float));
	Points.distances = (float*)malloc(n*sizeof(float)); // By Guiming @ 2016-09-02

	int counter = 0;
	char * pch;
	float x, y;
	rewind(f); // go back to the beginning of file
	fgets(line, CSV_LINE_LENGTH, f); //skip the header line
	while (fgets(line, CSV_LINE_LENGTH, f) != NULL){
		pch = strtok(line, ",\n");
		x = atof(pch);
		while (pch != NULL)
		{
			pch = strtok(NULL, ",\n");
			y = atof(pch);
			break;
		}
		Points.xCoordinates[counter] = x;
		Points.yCoordinates[counter] = y;
		Points.weights[counter] = 1.0f;
		Points.distances[counter] = 0.0f; // By Guiming @ 2016-09-02

		counter++;
	}

	fclose(f);

	return Points;
}

void AllocateDeviceSamplePoints(SamplePoints* dPoints, const SamplePoints Points){
	//Timothy @ 01/15/2021
	//EDIT: Changing dPoints to be a array of pointers to each set of points on each device.
	for (int i = 0; i < GPU_N; i++)
	{
		dPoints[i] = Points;
	
		dPoints[i].numberOfPoints = Points.numberOfPoints;
		dPoints[i].totalNumPoints = Points.numberOfPoints;
		int size = Points.numberOfPoints * sizeof(float);
		cudaError_t error;

		cudaSetDevice(i);
		error = cudaMalloc((void**)&dPoints[i].xCoordinates, size);
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMalloc((void**)&dPoints[i].yCoordinates, size);
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMalloc((void**)&dPoints[i].weights, size);
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		// By Guiming @ 2016-09-02
		error = cudaMalloc((void**)&dPoints[i].distances, size);
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

//original
//void CopyToDeviceSamplePoints(SamplePoints* dPoints, const SamplePoints hPoints) {
//	int size = hPoints.numberOfPoints * sizeof(float);
//
//	//for(int i = 0; i < hPoints.numberOfPoints; i++)
//	//	printf("x:%.2f y:%.2f w:%.2f\n", hPoints.xCoordinates[i], hPoints.yCoordinates[i], hPoints.weights[i]);
//
//	//printf("copy %d points to device\n", size);
//	cudaError_t error;
//
//	error = cudaMemcpy(dPoints[0].xCoordinates, hPoints.xCoordinates, size, cudaMemcpyHostToDevice);
//	if (error != cudaSuccess)
//	{
//		printf("ERROR in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
//		exit(EXIT_FAILURE);
//	}
//	error = cudaMemcpy(dPoints[0].yCoordinates, hPoints.yCoordinates, size, cudaMemcpyHostToDevice);
//	if (error != cudaSuccess)
//	{
//		printf("ERROR in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
//		exit(EXIT_FAILURE);
//	}
//	error = cudaMemcpy(dPoints[0].weights, hPoints.weights, size, cudaMemcpyHostToDevice);
//	if (error != cudaSuccess)
//	{
//		printf("ERROR in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
//		exit(EXIT_FAILURE);
//	}
//
//	// By Guiming @ 2016-09-02
//	error = cudaMemcpy(dPoints[0].distances, hPoints.distances, size, cudaMemcpyHostToDevice);
//	if (error != cudaSuccess)
//	{
//		printf("ERROR in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
//		exit(EXIT_FAILURE);
//	}
//}


//EDIT: Timothy @ 03/26/2021
//Added additional variable to track division of points across multiple GPUs
void CopyToDeviceSamplePoints(SamplePoints* dPoints, const SamplePoints hPoints) {
	int size = hPoints.numberOfPoints * sizeof(float);
	int n = hPoints.numberOfPoints; //Number of points on GPU
	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
	int div = n / GPU_N; //Division of points to be divided amongst GPUs
	int divNum = 0; //Tracks our place in the original n number of points
	int index = 0; //Tracks indexing for our multiple GPUs
	cudaError_t error;
	dPoints[0].end = div;

	//Timothy @ 01/15/2020
	//Copying the points to each GPU so the data is present across all devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i); //Set device (GPU) being actively copied to
		dPoints[i].start = index; //Begin tracking division of points
		dPoints[i].end = index + div; //Tracking end of division

		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
		//if not add remainder to size on final GPU
		if ((i == GPU_N - 1) && (rem != 0)) 
		{
			div += rem;
		}
		dPoints[i].numberOfPoints = div; //# of points is assigned here to compensate for remainders.
		error = cudaMemcpy(dPoints[i].xCoordinates, hPoints.xCoordinates, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 1 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMemcpy(dPoints[i].yCoordinates, hPoints.yCoordinates, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 2 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMemcpy(dPoints[i].weights, hPoints.weights, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 3 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		error = cudaMemcpy(dPoints[i].distances, hPoints.distances, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 4 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		index = div; //Set starting index of next group of sample points to the end of previous group.
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyFromDeviceSamplePoints(SamplePoints hPoints, const SamplePoints* dPoints){
	int size = dPoints[0].numberOfPoints * sizeof(float);
	cudaError_t error;

	cudaSetDevice(0);

	error = cudaMemcpy(hPoints.xCoordinates, dPoints[0].xCoordinates, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
    {
        printf("ERROR 1 in CopyFromDeviceSamplePoints: %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	error = cudaMemcpy(hPoints.yCoordinates, dPoints[0].yCoordinates, size, cudaMemcpyDeviceToHost);
		if (error != cudaSuccess)
    {
        printf("ERROR 2 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	error = cudaMemcpy(hPoints.weights, dPoints[0].weights, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
    {
        printf("ERROR 3 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	// By Guiming @ 2016-09-02
	error = cudaMemcpy(hPoints.distances, dPoints[0].distances, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
	{
	    printf("ERROR 4 in CopyToDeviceSamplePoints: %s\n", cudaGetErrorString(error));
	    exit(EXIT_FAILURE);
	}
}

// write to .csv file
void WriteSamplePoints(SamplePoints* Points, const char * csvFile){
	FILE *f = fopen(csvFile, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "x, y\n");
	for (int p = 0; p < Points->numberOfPoints; p++){
		fprintf(f, "%f, %f\n", Points->xCoordinates[p], Points->yCoordinates[p]);
	}
	fclose(f);
}

// write to .csv file
void WriteSamplePoints(SamplePoints* Points, float* Hs, float* Ws, const char * csvFile){
	FILE *f = fopen(csvFile, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "x, y, h, w\n");
	for (int p = 0; p < Points->numberOfPoints; p++){
		fprintf(f, "%f, %f, %f, %f\n", Points->xCoordinates[p], Points->yCoordinates[p], Hs[p], Ws[p]);
	}
	fclose(f);
}

void FreeSamplePoints(SamplePoints* Points) {
	free(Points->xCoordinates);
	Points->xCoordinates = NULL;

	free(Points->yCoordinates);
	Points->yCoordinates = NULL;

	free(Points->weights);
	Points->weights = NULL;
	
	// By Guiming @ 2016-09-02
	free(Points->distances);
	Points->distances = NULL;
}

void FreeDeviceSamplePoints(SamplePoints* dPoints){
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Free Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaFree(dPoints[i].xCoordinates);
		if (error != cudaSuccess)
		{
			printf("ERROR 1 in FreeDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		dPoints->xCoordinates = NULL;

		error = cudaFree(dPoints[i].yCoordinates);
		if (error != cudaSuccess)
		{
			printf("ERROR 2 in FreeDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		dPoints->yCoordinates = NULL;

		error = cudaFree(dPoints[i].weights);
		if (error != cudaSuccess)
		{
			printf("ERROR 3 in FreeDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		dPoints->weights = NULL;

		// By Guiming @ 2016-09-02
		error = cudaFree(dPoints[i].distances);
		if (error != cudaSuccess)
		{
			printf("ERROR in FreeDeviceSamplePoints: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		dPoints->distances = NULL;
	}
	cudaSetDevice(0); //Reset device to first GPU
}

// this is a mask
AsciiRaster AllocateAsciiRaster(int nCols, int nRows, float xLLCorner, float yLLCorner, float cellSize, float noDataValue){
	AsciiRaster Ascii;

	Ascii.nCols = nCols;
	Ascii.nRows = nRows;
	Ascii.xLLCorner = xLLCorner;
	Ascii.yLLCorner = yLLCorner;
	Ascii.cellSize = cellSize;
	Ascii.noDataValue = noDataValue;

	int size = Ascii.nCols * Ascii.nRows;
	Ascii.elements = (float*)malloc(size * sizeof(float));

	for (int row = 0; row < Ascii.nRows; row++){
		for (int col = 0; col < Ascii.nCols; col++){
			//if (row < 2 || col < 2)
			//	Ascii.elements[row * nCols + col] = Ascii.noDataValue;
			//else
				Ascii.elements[row * nCols + col] = 0.0f;
		}
	}

	return Ascii;
}

// copy a ascii raster
AsciiRaster CopyAsciiRaster(const AsciiRaster anotherAscii){
	AsciiRaster Ascii;

	Ascii.nCols = anotherAscii.nCols;
	Ascii.nRows = anotherAscii.nRows;
	Ascii.xLLCorner = anotherAscii.xLLCorner;
	Ascii.yLLCorner = anotherAscii.yLLCorner;
	Ascii.cellSize = anotherAscii.cellSize;
	Ascii.noDataValue = anotherAscii.noDataValue;

	int size = Ascii.nCols * Ascii.nRows;
	Ascii.elements = (float*)malloc(size * sizeof(float));

	for (int row = 0; row < Ascii.nRows; row++){
		for (int col = 0; col < Ascii.nCols; col++){
			Ascii.elements[row * Ascii.nCols + col] = anotherAscii.elements[row * Ascii.nCols + col];
		}
	}

	return Ascii;
}

// ascii raster read from a .asc file
AsciiRaster ReadAsciiRaster(char * asciiFile){
	FILE *f = fopen(asciiFile, "rt");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	AsciiRaster Ascii;

	const int HEADER_LINE_LENGTH = 64;
	char hdrLine[HEADER_LINE_LENGTH];
	char* pch;
	float meta[6] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

	// read headers
	for (int i = 0; i < 6; i++){
		fgets(hdrLine, HEADER_LINE_LENGTH, f);
		pch = strtok(hdrLine, " \n");
		while (pch != NULL)
		{
			pch = strtok(NULL, "\n");
			meta[i] = atof(pch);
			break;
		}
	}

	Ascii.nCols = (int)meta[0];
	Ascii.nRows = (int)meta[1];
	Ascii.xLLCorner = meta[2];
	Ascii.yLLCorner = meta[3];
	Ascii.cellSize = meta[4];
	Ascii.noDataValue = meta[5];
	Ascii.elements = (float*)malloc(Ascii.nRows * Ascii.nCols * sizeof(float));

	const int DATA_LINE_LENGTH = Ascii.nCols * 32;
	char* datLine = (char*)malloc(DATA_LINE_LENGTH * sizeof(char));

	int row_counter = 0;
	while (fgets(datLine, DATA_LINE_LENGTH, f) != NULL){
		int col_counter = 0;
		pch = strtok(datLine, " \n");
		Ascii.elements[row_counter*Ascii.nCols+col_counter] = atof(pch);
		while (pch != NULL)
		{
			pch = strtok(NULL, " ");
			if (pch != NULL && col_counter < Ascii.nCols - 1){
				col_counter++;
				Ascii.elements[row_counter*Ascii.nCols + col_counter] = atof(pch);
			}
		}
		row_counter++;
	}
	free(datLine);

	fclose(f);

	return Ascii;
}

void AllocateDeviceAsciiRaster(AsciiRaster* dAscii, const AsciiRaster hAscii){
	//Timothy @ 10/16/2020
	//Allocate Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		dAscii[i].nCols = hAscii.nCols;
		dAscii[i].nRows = hAscii.nRows;
		dAscii[i].xLLCorner = hAscii.xLLCorner;
		dAscii[i].yLLCorner = hAscii.yLLCorner;
		dAscii[i].cellSize = hAscii.cellSize;
		dAscii[i].noDataValue = hAscii.noDataValue;
		int size = hAscii.nCols*hAscii.nRows * sizeof(float);
		cudaError_t error;
	
			cudaSetDevice(i);
			error = cudaMalloc((void**)&dAscii[i].elements, size);
			if (error != cudaSuccess)
			{
				printf("ERROR in AllocateDeviceAsciiRaster: %s\n", cudaGetErrorString(error));
				exit(EXIT_FAILURE);
			}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyToDeviceAsciiRaster(AsciiRaster* dAscii, const AsciiRaster hAscii){
	int size = hAscii.nCols*hAscii.nRows * sizeof(float);
	cudaError_t error;
	//Copy raster to all available devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dAscii[i].elements, hAscii.elements, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR in CopyToDeviceAsciiRaster: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyFromDeviceAsciiRaster(AsciiRaster hAscii, const AsciiRaster dAscii){
	hAscii.nCols = dAscii.nCols;
	hAscii.nRows = dAscii.nRows;
	hAscii.xLLCorner = dAscii.xLLCorner;
	hAscii.yLLCorner = dAscii.yLLCorner;
	hAscii.cellSize = dAscii.cellSize;
	hAscii.noDataValue = dAscii.noDataValue;

	int size = dAscii.nCols*dAscii.nRows * sizeof(float);
	cudaError_t error;
	error = cudaMemcpy(hAscii.elements, dAscii.elements, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
    {
        printf("ERROR in CopyFromDeviceAsciiRaster: %s\n", cudaGetErrorString(error));
		printf("size=%d mode=%d\n", size, cudaMemcpyDeviceToHost);
        exit(EXIT_FAILURE);
    }
}

// write to .asc file
void WriteAsciiRaster(AsciiRaster* Ascii, const char * asciiFile){
	FILE *f = fopen(asciiFile, "w");
	if (f == NULL)
	{
		printf("Error opening file!\n");
		exit(1);
	}

	fprintf(f, "ncols %d\n", Ascii->nCols);
	fprintf(f, "nrows %d\n", Ascii->nRows);
	fprintf(f, "xllcorner %f\n", Ascii->xLLCorner);
	fprintf(f, "yllcorner %f\n", Ascii->yLLCorner);
	fprintf(f, "cellsize %f\n", Ascii->cellSize);
	fprintf(f, "NODATA_value %.0f\n", Ascii->noDataValue);

	for (int row = 0; row < Ascii->nRows; row++){
		for (int col = 0; col < Ascii->nCols; col++){
			fprintf(f, "%.16f ", Ascii->elements[row*Ascii->nCols+col]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

void FreeAsciiRaster(AsciiRaster* Ascii){
	free(Ascii->elements);
	Ascii->elements = NULL;
}

void FreeDeviceAsciiRaster(AsciiRaster* Ascii){
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Free Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaFree(Ascii[i].elements);
		if (error != cudaSuccess)
		{
			printf("ERROR in FreeDeviceAsciiRaster: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		Ascii[i].elements = NULL;
	}
	cudaSetDevice(0); //Reset device to first GPU
}

// edge effects correction weights at each point, weights is allocated somewhere else
void EdgeCorrectionWeightsExact(SamplePoints Points, float h, AsciiRaster Ascii, float *weights){
	float h2 = h * h;
	float cellArea = Ascii.cellSize * Ascii.cellSize;
	float p_x, p_y, cell_x, cell_y;
	float ew;

	for (int p = 0; p < Points.numberOfPoints; p++){
		//printf("%6d / %6d\n", p, Points.numberOfPoints);

		// By Guiming @ 2016-09-03
		if(Points.distances[p] >= CUT_OFF_FACTOR * h2){ // pnts too far away from the study area boundary, skip to save labor!
			weights[p] = 1.0f;
			//printf("bypassed! %f %f %d\n", Points.distances[p], 9.0 * h2, nThreads);
			continue;
		}

		p_x = Points.xCoordinates[p];
		p_y = Points.yCoordinates[p];
		ew = 0.0f;

		// added by Guiming @2016-09-11
		// narrow down the row/col range
		int row_lower = 0;
		int row_upper = Ascii.nRows - 1;
		int col_lower = 0;
		int col_upper = Ascii.nCols - 1;
		if(NARROW){
			int r = YCOORD_TO_ROW(p_y + SQRT_CUT_OFF_FACTOR * h, Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize);
			row_lower = MAX(0, r);
			row_upper = MIN(Ascii.nRows - 1, YCOORD_TO_ROW(p_y - SQRT_CUT_OFF_FACTOR * h, Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize));
			col_lower = MAX(0, XCOORD_TO_COL(p_x - SQRT_CUT_OFF_FACTOR * h, Ascii.xLLCorner, Ascii.cellSize));
			col_upper = MIN(Ascii.nCols - 1, XCOORD_TO_COL(p_x + SQRT_CUT_OFF_FACTOR * h, Ascii.xLLCorner, Ascii.cellSize));
		}

		for (int row = row_lower; row <= row_upper; row++){
			for (int col = col_lower; col <= col_upper; col++){
				if (Ascii.elements[row*Ascii.nCols+col] != Ascii.noDataValue){
					cell_x = COL_TO_XCOORD(col, Ascii.xLLCorner, Ascii.cellSize);
					cell_y = ROW_TO_YCOORD(row, Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize);
					float d2 = Distance2(p_x, p_y, cell_x, cell_y);
					ew += GaussianKernel(h2, d2) * cellArea;
				}
			}
		}
		weights[p] = 1.0 / ew;
	}
}

void EdgeCorrectionWeightsExact(SamplePoints Points, float* hs, AsciiRaster Ascii, float *weights){
	//float h2 = BandWidth2(Points);
	float cellArea = Ascii.cellSize * Ascii.cellSize;
	float p_x, p_y, cell_x, cell_y;
	float ew, h2;

	for (int p = 0; p < Points.numberOfPoints; p++){
		//printf("%6d / %6d\n", p, Points.numberOfPoints);
		p_x = Points.xCoordinates[p];
		p_y = Points.yCoordinates[p];
		ew = 0.0f;
		h2 = hs[p] * hs[p];

		// By Guiming @ 2016-09-03
		if(Points.distances[p] >= CUT_OFF_FACTOR * h2){ // pnts too far away from the study area boundary, skip to save labor!
			weights[p] = 1.0f;
			//printf("bypassed! %f %f %d\n", Points.distances[p], 9.0 * h2, nThreads);
			continue;
		}

		// added by Guiming @2016-09-11
		// narrow down the row/col range
		int row_lower = 0;
		int row_upper = Ascii.nRows - 1;
		int col_lower = 0;
		int col_upper = Ascii.nCols - 1;

		if(NARROW){
			int r = YCOORD_TO_ROW(p_y + SQRT_CUT_OFF_FACTOR * hs[p], Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize);
			row_lower = MAX(0, r);
			row_upper = MIN(Ascii.nRows - 1, YCOORD_TO_ROW(p_y - SQRT_CUT_OFF_FACTOR * hs[p], Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize));
			col_lower = MAX(0, XCOORD_TO_COL(p_x - SQRT_CUT_OFF_FACTOR * hs[p], Ascii.xLLCorner, Ascii.cellSize));
			col_upper = MIN(Ascii.nCols - 1, XCOORD_TO_COL(p_x + SQRT_CUT_OFF_FACTOR * hs[p], Ascii.xLLCorner, Ascii.cellSize));
		}

		for (int row = row_lower; row <= row_upper; row++){
			for (int col = col_lower; col <= col_upper; col++){
				if (Ascii.elements[row*Ascii.nCols+col] != Ascii.noDataValue){
					cell_x = COL_TO_XCOORD(col, Ascii.xLLCorner, Ascii.cellSize);
					cell_y = ROW_TO_YCOORD(row, Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize);
					float d2 = Distance2(p_x, p_y, cell_x, cell_y);
					ew += GaussianKernel(h2, d2) * cellArea;
				}
			}
		}
		weights[p] = 1.0 / ew;
	}
}

float* AllocateEdgeCorrectionWeights(SamplePoints Points){
	return (float*)malloc(Points.numberOfPoints*sizeof(float));
}

void AllocateDeviceEdgeCorrectionWeights(float** dWeights, SamplePoints Points){
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Allocate Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMalloc((void**)&dWeights[i], Points.numberOfPoints * sizeof(float));
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceEdgeCorrectionWeights: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void FreeEdgeCorrectionWeights(float* weights){
	
	free(weights);
	weights = NULL;
}

void FreeDeviceEdgeCorrectionWeights(float** weights){
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Free Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaFree(weights[i]);
		if (error != cudaSuccess)
		{
			printf("ERROR in FreeDeviceEdgeCorrectionWeights: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		weights[i] = NULL;
	}
	cudaSetDevice(0); //Reset device to first GPU
}

///////// Guiming on 2016-03-16 ///////////////
// the array holding bandwidth at each point
float* AllocateBandwidths(int n){ // n is number of points
	return (float*)malloc(n*sizeof(float));
}

void AllocateDeviceBandwidths(float** dBandwidths, int n){ // n is number of points
	cudaError_t error;
	//Allocate bandwidth accross all available devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMalloc((void**)&dBandwidths[i], n * sizeof(float));
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceBandwidths: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyToDeviceBandwidths(float** dBandwidth, const float* hBandwidths, const int n) {
	int size = n * sizeof(float);
	cudaError_t error;
	
	//Copy to each available device
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dBandwidth[i], hBandwidths, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR in CopyToDeviceBandwidths: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0);
}

void CopyFromDeviceBandwidths(float* hBandwidth, const float* dBandwidths, const int n){
	int size = n * sizeof(float);
	cudaError_t error;
	error = cudaMemcpy(hBandwidth, dBandwidths, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
    {
        printf("ERROR in CopyFromDeviceBandwidths: %s\n", cudaGetErrorString(error));
		printf("size=%d mode=%d\n", size, cudaMemcpyDeviceToHost);
        exit(EXIT_FAILURE);
    }
}

void FreeDeviceBandwidths(float* bandwidths){
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Free Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaFree(bandwidths);
		if (error != cudaSuccess)
		{
			printf("ERROR in FreeDeviceBandwidths: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		bandwidths = NULL;
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void FreeBandwidths(float* bandwidths){
	free(bandwidths);
	bandwidths = NULL;
}

// the array holding inclusive density at each point
float* AllocateDen(int n){ // n is number of points
	return (float*)malloc(n*sizeof(float));
}

void AllocateDeviceDen(float** dDen, int n){ // n is number of points
	cudaError_t error;
	//Timothy @ 10/16/2020
	//Allocate Memory Across All Devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMalloc((void**)&dDen[i], n * sizeof(float));
		if (error != cudaSuccess)
		{
			printf("ERROR in AllocateDeviceDen: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyToDeviceDen(float** dDen, const float* hDen, const int n){
	int size = n * sizeof(float);
	cudaError_t error;
	//Copy accross all available devices
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dDen[i], hDen, size, cudaMemcpyHostToDevice);
	}
	if (error != cudaSuccess)
    {
        printf("ERROR in CopyToDeviceDen: %s\n", cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
	cudaSetDevice(0); //Reset device to first GPU
}

void CopyFromDeviceDen(float* hDen, const float* dDen, const int n){
	int size = n * sizeof(float);
	cudaError_t error;
	error = cudaMemcpy(hDen, dDen, size, cudaMemcpyDeviceToHost);
	if (error != cudaSuccess)
    {
        printf("ERROR in CopyFromDeviceDen: %s\n", cudaGetErrorString(error));
		printf("size=%d mode=%d\n", size, cudaMemcpyDeviceToHost);
        exit(EXIT_FAILURE);
    }
}

void CopyDeviceDen(float* dDenTo, const float* dDenFrom, const int n){
	int size = n * sizeof(float);
	cudaError_t error = cudaSuccess;
	error = cudaMemcpy(dDenTo, dDenFrom, size, cudaMemcpyDeviceToDevice);
	if (error != cudaSuccess)
    {
        printf("ERROR in CopyDeviceDen: %s\n", cudaGetErrorString(error));
		printf("size=%d mode=%d\n", size, cudaMemcpyDeviceToDevice);
        exit(EXIT_FAILURE);
    }
}

void FreeDeviceDen(float** den){
	cudaError_t error;
	cudaSetDevice(0);
	error = cudaFree(den[0]);
	if (error != cudaSuccess)
	{
		printf("ERROR in FreeDeviceDen(Elements): %s\n", cudaGetErrorString(error));
		exit(EXIT_FAILURE);
	}
	den = NULL;
	error = cudaFree(den);
	if (error != cudaSuccess)
	{
		printf("ERROR in FreeDeviceDen(Pointer): %s\n", cudaGetErrorString(error));
		exit(EXIT_FAILURE);
	}
}

void FreeDen(float* den){
	free(den);
	den = NULL;
}

// compute the optimal Maximum Likelihood Estimation fixed bandwidth
// By Guiming @ 2016-02-26
float MLE_FixedBandWidth(AsciiRaster* Ascii, SamplePoints* Points, float **edgeWeights, float h, float* den0, float* den1, bool useGPU, float** dDen0, float** dDen1){
	
	float hA = h/10;
	float hD = 4 * h;
	float width = hD - hA;
	float epsilon = width/100;
	float factor = 1 + sqrtf(5.0f);
	int iteration = 0;

	printf("hA: %f hD: %f width: %f, epsilon: %d\n", hA, hD, width, epsilon); //DEBUG
	while(width > epsilon){

		if(DEBUG){
			printf("iteration: %d ", iteration);
			printf("hD: %.6f ", hD);
			printf("hA: %.6f ", hA);
		}

		float hB = hA + width / factor;
		float hC = hD - width / factor;

		//ERROR HERE, ONLY WHEN USING GPU
		float LoghB = LogLikelihood(Ascii, Points, edgeWeights, hB, den0, den1, useGPU, dDen0, dDen1);
		float LoghC = LogLikelihood(Ascii, Points, edgeWeights, hC, den0, den1, useGPU, dDen0, dDen1);

		if(LoghB > LoghC){
			hD = hC;
			if(DEBUG) printf("LoghB: %.6f \n", LoghB);
		}
		else{
			hA = hB;
			if(DEBUG) printf("LoghC: %.6f \n", LoghC);
		}

		width = hD - hA;

		iteration += 1;
	}

	return (hA + hD) / 2;
}

// By Guiming @ 2016-05-21
// computed fixed bandwidth kde
void ComputeFixedDensityAtPoints(AsciiRaster Ascii, SamplePoints Points, float* edgeWeights, float h, float* den0, float* den1, float* dDen0, float* dDen1) {
	
	int numPoints = Points.numberOfPoints;
		// update edge correction weights
		if (UPDATEWEIGHTS) {
			EdgeCorrectionWeightsExact(Points, h, Ascii, edgeWeights);
		}

		for (int i = 0; i < numPoints; i++) {
			float pi_x = Points.xCoordinates[i];
			float pi_y = Points.yCoordinates[i];

			float den = EPSILONDENSITY;
			float den_itself = EPSILONDENSITY;
			for (int j = 0; j < numPoints; j++) {
				float pj_x = Points.xCoordinates[j];
				float pj_y = Points.yCoordinates[j];
				float pj_w = Points.weights[j];
				float pj_ew = edgeWeights[j];

				float d2 = Distance2(pi_x, pi_y, pj_x, pj_y);

				if (j == i) {
					den_itself += GaussianKernel(h * h, d2) * pj_w * pj_ew; // / numPoints;
				}
				else {
					den += GaussianKernel(h * h, d2) * pj_w * pj_ew;
				}
			}

			if (den0 != NULL) {
				den0[i] = den + den_itself;
			}
			if (den1 != NULL) {
				den1[i] = den;
			}
		}
}

// By Timothy @ 04-23-2021
// Separated this function into two separate, rather than using the boolean value as was done previously.
//this was done to enable parallel functionality accross multiple GPUs
void ComputeFixedDensityDevice(cudaStream_t* streams, AsciiRaster* Ascii, SamplePoints* Points, float** edgeWeights, float h, float* den0, float* den1, float** dDen0, float** dDen1){
	int numPoints = Points[0].numberOfPoints;
	// invoke kernels to compute density at each point
	// execution config.
	int NBLOCK_W = (numPoints + BLOCK_SIZE - 1) / BLOCK_SIZE;
	int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
	dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);
	for (int i = 0; i < GPU_N; i++){
		cudaSetDevice(i);
		// update edge correction weights
		if (UPDATEWEIGHTS) {
			CalcEdgeCorrectionWeights <<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(h * h, Points[i], Ascii[i], edgeWeights[i]);
		}
		//// brute force to search for neighbors
		//DensityAtPoints<<<dimGrid_W, BLOCK_SIZE>>>(h*h, Points, edgeWeights, dDen0, dDen1);
		///// use KD Tree to speedup neighor search
		//CopyToDeviceDen(gpuDen, zeroDen, Points.numberOfPoints);
		DensityAtPointsKdtr <<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(GPU_tree[i].m_gpu_nodes, GPU_tree[i].m_gpu_indexes, GPU_tree[0].m_gpu_points, h * h, Points[i], edgeWeights[i], gpuDen[i]);
		
		// have to do this as a separate kernel call due to the need of block synchronization !!!
		// this took me hours to debug!
		dCopyDensityValues <<<dimGrid_W, BLOCK_SIZE, 0, streams[i]>>>(Points[i], edgeWeights[i], h * h, gpuDen[i], dDen0[i], dDen1[i]);
		cudaStreamSynchronize(streams[i]);
	}
	cudaSetDevice(0);
}

// By Guiming @ 2016-02-26
// the log likelihood given single bandwidth h
float LogLikelihood(AsciiRaster* Ascii, SamplePoints* Points, float **edgeWeights, float h, float* den0, float* den1, bool useGPU, float** dDen0, float** dDen1){
	int numPoints = Points[0].numberOfPoints;
	float logL = 0.0f; // log likelihood
	int size = Points[0].totalNumPoints * sizeof(float);
	cudaError error = cudaSuccess;

	if (error != cudaSuccess)
	{
		printf("ERROR 0 in LogLikelihood: %s\n", cudaGetErrorString(error));
		exit(EXIT_FAILURE);
	}
	if (useGPU) { // do it on GPU
		//SamplePoints hostP = AllocateSamplePoints(Points[0].numberOfPoints);
		for (int i = 0; i < GPU_N; i++)
		{
			cudaSetDevice(i);
			///*
			// execution config.
			int NBLOCK_W = (numPoints + BLOCK_SIZE - 1) / BLOCK_SIZE;
			int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
			dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);

			// update edge correction weights
			if (UPDATEWEIGHTS) {
				CalcEdgeCorrectionWeights << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (h * h, Points[i], Ascii[i], edgeWeights[i]);
			}
			cudaStreamSynchronize(streams[i]);
			DensityAtPointsKdtr << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (GPU_tree[i].m_gpu_nodes, GPU_tree[i].m_gpu_indexes, GPU_tree[i].m_gpu_points, h * h, Points[i], edgeWeights[i], gpuDen[i]);
			// have to do this as a separate kernel call due to the need of block synchronization !!!
			// this took me hours to debug!
			cudaStreamSynchronize(streams[i]);
			dCopyDensityValues << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (Points[i], edgeWeights[i], h * h, gpuDen[i], NULL, dDen1[i]);
			cudaStreamSynchronize(streams[i]);
		}
		ReformPoints(Points);
		// compute likelihood on GPU
		ReductionSumGPU(dDen1[0], numPoints);
		cudaMemcpyFromSymbol(&logL, dReductionSum, sizeof(float), 0, cudaMemcpyDeviceToHost);
		//printf("reduction result (likelihood) A: %3.4f \n", logL);
		CopyToDeviceSamplePoints(Points, sPoints);
		//Cleanup
		cudaSetDevice(0);
		/*FreeSamplePoints(&hostP);*/
	}
	else{ // do it on CPU
		// update edge correction weights
		if(UPDATEWEIGHTS){
			EdgeCorrectionWeightsExact(Points[0], h, Ascii[0], edgeWeights[0]);
		}

		// the kd tree appraoch
		float* tmpden = AllocateDen(numPoints);
		float h2 = h * h;
		float range = CUT_OFF_FACTOR * h2;

		for(int i = 0; i < numPoints; i++){
			tmpden[i] = -1.0 * GaussianKernel(h2, 0.0f) *  Points[0].weights[i] * edgeWeights[0][i];
		}

		vector<int> ret_index = vector<int>();
		vector<float> ret_dist = vector<float>(); // squared distance

		for(int i = 0; i < numPoints; i++){
			float pi_x = Points[0].xCoordinates[i];
			float pi_y = Points[0].yCoordinates[i];
			float pj_w = Points[0].weights[i];
			float pj_ew = edgeWeights[0][i];

			// range query
			Point query;
			query.coords[0] = pi_x;
			query.coords[1] = pi_y;
			ret_index.clear();
			ret_dist.clear();
			tree.SearchRange(query, range, ret_index, ret_dist);
			//printf("CPU PNT_%d %d NBRS RANGE=%.1f\n", i, ret_index.size(), range);

			if(ret_index.size() > MAX_N_NBRS) MAX_N_NBRS = ret_index.size();

			float g = 0.0f;
			int idx;
			for(int j = 0; j < ret_index.size(); j++){
					g = GaussianKernel(h2, ret_dist[j]) * pj_w *pj_ew;
					idx = ret_index[j];
					//float t = tmpden[idx];
					tmpden[idx] += g;
					//if(i == 0) printf("CPU PNT_%d g[%d]=%.5f gpuDen[%d]=%.5f gpuDen[%d]=%.5f\n", i, idx, g, idx, t, idx, tmpden[idx]);
			}
		} // END OF COMPUTING DENSITIES AT POINTS



		for(int i = 0; i < numPoints; i++){
			//printf("CPU H2=%.2f DEN[%d]=%.5f\n", h2, i, tmpden[i]);
			logL += logf(tmpden[i] + EPSILONDENSITY);
		}

		if(den0 != NULL){
			for(int i = 0; i < numPoints; i++)
				den0[i] = tmpden[i]  + GaussianKernel(h2, 0.0f) * Points[0].weights[i] * edgeWeights[0][i];
		}
		if(den1 != NULL){
			for(int i = 0; i < numPoints; i++)
				den1[i] = tmpden[i];
		}

		FreeDen(tmpden);
	}
	return logL;
}

// the log likelihood given bandwidths hs
// By Guiming @ 2016-02-26
// float* den0 : density based on all points, including itself
// float* den1 : leave one out density
//EDIT: Timothy @ 03/12/2021
//Added additional variable passed in, so now when useGPU is TRUE, the function will be handled across however many GPUs are available
float LogLikelihood(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float **edgeWeights, float* hs, float* den0, float* den1, bool useGPU, float** dHs, float** dDen0, float** dDen1, float h, float alpha, float** dDen0cpy){
	int numPoints = Points[0].numberOfPoints;
	float logL = 0.0f; // log likelihood
	cudaError_t error = cudaSuccess;
	// execution config.
	int NBLOCK_W = (numPoints + BLOCK_SIZE - 1) / BLOCK_SIZE;
	int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
	dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);
	if(useGPU){ // do it on GPU
		for (int i = 0; i < GPU_N; i++)
		{
			cudaSetDevice(i);
			//CopyToDeviceBandwidths(dHs, hs, numPoints);

			// update bandwidth on GPU
			//CalcVaryingBandwidths<<<dimGrid_W, BLOCK_SIZE>>>(Points, dDen0, h, alpha, dHs);
			// update edge correction weights
			if (UPDATEWEIGHTS) {
				CalcEdgeCorrectionWeights << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (h * h, gpuPoints[i], Ascii[i], edgeWeights[i]);
			}

			// compute (log) density at sample points [h^2, not h! OMG!!! Took me hours for spotting this!]
			//DensityAtPoints<<<dimGrid_W, BLOCK_SIZE>>>(h * h, Points, edgeWeights, dDen0, dDen1);

			//CopyToDeviceDen(gpuDen, zeroDen, Points.numberOfPoints);
			DensityAtPointsKdtr << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (GPU_tree[i].m_gpu_nodes, GPU_tree[i].m_gpu_indexes, GPU_tree[i].m_gpu_points, h * h, gpuPoints[i], edgeWeights[i], gpuDen[i]);
			// have to do this as a separate kernel call due to the need of block synchronization !!!
			// this took me hours to debug!
			dCopyDensityValues << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (gpuPoints[i], edgeWeights[i], h * h, gpuDen[i], dDen0[i], dDen1[i]);
			
			cudaStreamSynchronize(streams[i]);
		}

		ReformDensities(dDen0, den0);
		 

		CopyDeviceDen(dDen0cpy[0], dDen0[0], sPoints.numberOfPoints);
		
		int size = gpuPoints[0].totalNumPoints * sizeof(float);
		//SamplePoints hostP = AllocateSamplePoints(Points[0].numberOfPoints);
		if (error != cudaSuccess)
		{
			printf("ERROR 0 in LogLikelihood (Overload): %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		ReformPoints(gpuPoints);
		
		// compute sum of log densities on GPU
		ReductionSumGPU(dDen0cpy[0], gpuPoints[0].totalNumPoints);
		//float tmp = 0.0f;
		//cudaMemcpyFromSymbol(&tmp, dReductionSum, sizeof(float), 0, cudaMemcpyDeviceToHost);
		//printf("reduction result (geometricmean): %3.4f \n", exp(tmp/numPoints));
		
		CopyToDeviceSamplePoints(gpuPoints, sPoints);
		
		for (int i = 0; i < GPU_N; i++)
		{
			cudaSetDevice(i);
			// update bandwidth on GPU
			CalcVaryingBandwidths << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (gpuPoints[i], dDen0[i], h, alpha, dHs[i]);

			// update edge correction weights
			if (UPDATEWEIGHTS) {
				CalcEdgeCorrectionWeights << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (dHs[i], gpuPoints[i], Ascii[i], edgeWeights[i]);
			}

			//DensityAtPoints<<<dimGrid_W, BLOCK_SIZE>>>(dHs, Points, edgeWeights, dDen0, dDen1);
			//CopyToDeviceDen(gpuDen, zeroDen, Points.numberOfPoints);

			//to-do: BANDWIDTH AND EDGECORRECTIONWEIGHTS NEEDS TO BE REFORMED HERE! Just like densities

			DensityAtPointsKdtr << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (GPU_tree[i].m_gpu_nodes, GPU_tree[i].m_gpu_indexes, GPU_tree[i].m_gpu_points, dHs[i], gpuPoints[i], edgeWeights[i], gpuDen[i]);
			// have to do this as a separate kernel call due to the need of block synchronization !!!
			// this took me hours to debug!
			dCopyDensityValues << <dimGrid_W, BLOCK_SIZE, 0, streams[i] >> > (gpuPoints[i], edgeWeights[i], dHs[i], gpuDen[i], dDen0[i], dDen1[i]);
		
			cudaStreamSynchronize(streams[i]);
		}
		ReformPoints(gpuPoints);
		// compute likelihood on GPU
		ReductionSumGPU(dDen1[0], numPoints);
		cudaMemcpyFromSymbol(&logL, dReductionSum, sizeof(float), 0, cudaMemcpyDeviceToHost);
		CopyToDeviceSamplePoints(gpuPoints, sPoints);
		//printf("reduction result (likelihood): %3.4f \n", logL);
		//Cleanup
		//FreeSamplePoints(&hostP);
		cudaSetDevice(0);
	}
	else{ // do it on CPU
		// update edge correction weights
		if(UPDATEWEIGHTS){
			EdgeCorrectionWeightsExact(Points[0], h, Ascii[0], edgeWeights[0]);
		}

		// kdtree approach
		float h2 = h * h;
		float range = CUT_OFF_FACTOR * h2;
		float* denTmp = AllocateDen(numPoints);
		for(int i = 0; i < numPoints; i++){
			denTmp[i] = 0.0f;
		}

		vector<int> ret_index = vector<int>();
		vector<float> ret_dist = vector<float>(); // squared distance
		for(int i = 0; i < numPoints; i++){
			float pi_x = Points[0].xCoordinates[i];
			float pi_y = Points[0].yCoordinates[i];
			float pj_w = Points[0].weights[i];
			float pj_ew = edgeWeights[0][i];

			// range query
			Point query;
			query.coords[0] = pi_x;
			query.coords[1] = pi_y;
			ret_index.clear();
			ret_dist.clear();
			tree.SearchRange(query, range, ret_index, ret_dist);

			if(ret_index.size() > MAX_N_NBRS) MAX_N_NBRS = ret_index.size();

			int nn = ret_index.size();
			float g = 0.0f;
			int idx;
			for(int j = 0; j < ret_index.size(); j++){
					g = GaussianKernel(h2, ret_dist[j]) * pj_w * pj_ew;
					idx = ret_index[j];
					denTmp[idx] += g;
			}
		} // END OF COMPUTING DENSITIES AT POINTS

		// update bandwidths
		//float gml = compGML(den0, numPoints);
		float gml = compGML(denTmp, numPoints);
	    for(int i = 0; i < numPoints; i++){
				hs[i] = h * powf((denTmp[i] / gml), alpha);
	    }

		// update edge correction weights
		if(UPDATEWEIGHTS){
			EdgeCorrectionWeightsExact(Points[0], hs, Ascii[0], edgeWeights[0]);
		}

		for(int i = 0; i < numPoints; i++){
			float h2 = hs[i] * hs[i];
			denTmp[i] = -1.0 * GaussianKernel(h2, 0.0f) *  Points[0].weights[i] * edgeWeights[0][i];
		}

		for(int i = 0; i < numPoints; i++){
			float pi_x = Points[0].xCoordinates[i];
			float pi_y = Points[0].yCoordinates[i];
			float pj_w = Points[0].weights[i];
			float pj_ew = edgeWeights[0][i];
			float h2 = hs[i] * hs[i];
			float range = CUT_OFF_FACTOR * h2;

			// range query
			Point query;
			query.coords[0] = pi_x;
			query.coords[1] = pi_y;
			ret_index.clear();
			ret_dist.clear();
			tree.SearchRange(query, range, ret_index, ret_dist);

			if(ret_index.size() > MAX_N_NBRS) MAX_N_NBRS = ret_index.size();

			int nn = ret_index.size();
			float g = 0.0f;
			int idx;
			for(int j = 0; j < ret_index.size(); j++){
					g = GaussianKernel(h2, ret_dist[j]) * pj_w * pj_ew;
					idx = ret_index[j];
					denTmp[idx] += g;
			}
		} // END OF COMPUTING DENSITIES AT POINTS

		for(int i = 0; i < numPoints; i++){
			logL += logf(denTmp[i] + EPSILONDENSITY);
		}

		if(den0 != NULL){
			for(int i = 0; i < numPoints; i++){
				float h2 = hs[i] * hs[i];
				den0[i] = denTmp[i] + GaussianKernel(h2, 0.0f) *  Points[0].weights[i] * edgeWeights[0][i];
			}
		}

		if(den1 != NULL){
			for(int i = 0; i < numPoints; i++){
				den1[i] = denTmp[i];
			}
		}

		FreeDen(denTmp);
	}

	return logL;
}

// compute the log likelihood given a center (h0, alpha0) and step (stepH, stepA)
// By Guiming @ 2016-03-06
/*
 return 9 elements log likelihood in float* logLs
**/
void hj_likelihood(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float** edgeWeights, float h0, float alpha0, float stepH, float stepA, int lastdmax, float* logLs, float* hs, float* den0, float* den1, bool useGPU, float** dHs, float** dDen0, float** dDen1, float** dDen0cpy)
{
    //int n = Points.numberOfPoints;

    //float gml;
    // the center (h0, alpha0)
    if(lastdmax == -1){ // avoid unnecessary [expensive] computation
	    //LogLikelihood(Ascii, Points, edgeWeights, h0, den0, den1, useGPU, dDen0, dDen1);
	    float L0 = LogLikelihood(Ascii, Points, gpuPoints, edgeWeights, hs, den0, den1, useGPU, dHs, dDen0, dDen1, h0, alpha0, dDen0cpy);
	    //printf("L0: %.5f\t", L0);
	    logLs[0] = L0;
	}
    // (h0 - stepH, alpha0)
    if(lastdmax != 2){ // avoid unnecessary [expensive] computation
	    //LogLikelihood(Ascii, Points, edgeWeights, h0 - stepH, den0, den1, useGPU, dDen0, dDen1);
	    float L1 = LogLikelihood(Ascii, Points, gpuPoints, edgeWeights, hs, den0, den1, useGPU, dHs, dDen0, dDen1, h0 - stepH, alpha0, dDen0cpy);
	    //printf("L1: %.5f\t", L1);
	    logLs[1] = L1;
	}
    // (h0 + stepH, alpha0)
    if(lastdmax != 1){
	    //LogLikelihood(Ascii, Points, edgeWeights, h0 + stepH, den0, den1, useGPU, dDen0, dDen1);
	    float L2 = LogLikelihood(Ascii, Points, gpuPoints, edgeWeights, hs, den0, den1, useGPU, dHs, dDen0, dDen1, h0 + stepH, alpha0, dDen0cpy);
	    //printf("L2: %.5f\t", L2);
	    logLs[2] = L2;
	}
    // (h0, alpha0 + stepA)
    if(lastdmax != 4){
	    //LogLikelihood(Ascii, Points, edgeWeights, h0, den0, den1, useGPU, dDen0, dDen1);
	    float L3 = LogLikelihood(Ascii, Points, gpuPoints, edgeWeights, hs, den0, den1, useGPU, dHs, dDen0, dDen1, h0, alpha0 + stepA, dDen0cpy);
	    //printf("L3: %.5f\t", L3);
	    logLs[3] = L3;
	}
    // (h0, alpha0 - stepA)
    if(lastdmax != 3){
	    //LogLikelihood(Ascii, Points, edgeWeights, h0, den0, den1, useGPU, dDen0, dDen1);
	    float L4 = LogLikelihood(Ascii, Points, gpuPoints, edgeWeights, hs, den0, den1, useGPU, dHs, dDen0, dDen1, h0, alpha0 - stepA, dDen0cpy);
	    //printf("L4: %.5f\n", L4);
	    logLs[4] = L4;
	}
}

// compute the optimal h and alpha (parameters for calculating the optimal adaptive bandwith)
// By Guiming @ 2016-03-06
/*
 return 3 optmal parameters in float* optParas (optH, optAlpha, LogLmax)
//EDIT: Timothy @ 03/10/2021
//Added aditional variable to this, hj_likelihood and LogLikelihood functions to handle array of SamplePoints whenever multiple GPUs are present
**/
void hooke_jeeves(AsciiRaster* Ascii, SamplePoints* Points, SamplePoints* gpuPoints, float **edgeWeights, float h0, float alpha0, float stepH, float stepA, float* optParas, float* hs, float* den0, float* den1, bool useGPU, float** dHs, float** dDen0, float** dDen1, float** dDen0cpy){
	float* Ls = (float*)malloc(5 * sizeof(float)); // remember to free at the end
	hj_likelihood(Ascii, Points, gpuPoints, edgeWeights, h0, alpha0, stepH, stepA, -1, Ls, hs, den0, den1, useGPU, dHs, dDen0, dDen1, dDen0cpy);

	float Lmax = Ls[0];

	float s = stepH / 20;
	float a = stepA / 20;

	int iteration = 0;
    while ((stepH > s || stepA > a) &&  iteration <= MAX_NUM_ITERATIONS){

        //float Lmax0 = Lmax;
        int dmax = 0;
        for(int i = 0; i < 5; i++){
            if(Ls[i] > Lmax){
            	Lmax = Ls[i];
                dmax = i;
            }
        }
        if(DEBUG)
        	printf ("iteration: %d center: (%.5f %.5f) steps: (%.5f %.5f) dmax: %d Lmax: %.5f \n", iteration, h0, alpha0, stepH, stepA, dmax, Lmax);

        if(dmax == 0){
            stepH = stepH / 2;
            stepA = stepA / 2;
        }

        else{
            if(dmax == 1){
                h0 = h0 - stepH;
                alpha0 = alpha0;
                Ls[2] = Ls[0]; // avoid unnecessary [expensive] computation
                Ls[0] = Ls[1];
            }
            if(dmax == 2){
                h0 = h0 + stepH;
                alpha0 = alpha0;
                Ls[1] = Ls[0];
                Ls[0] = Ls[2];
            }
            if (dmax == 3){
                h0 = h0;
                alpha0 = alpha0 + stepA;
                Ls[3] = Ls[0];
                Ls[0] = Ls[4];
            }
            if(dmax == 4){
                h0 = h0;
                alpha0 = alpha0 - stepA;
                Ls[3] = Ls[0];
                Ls[0] = Ls[4];
            }
        }
	    hj_likelihood(Ascii, Points, gpuPoints, edgeWeights, h0, alpha0, stepH, stepA, dmax, Ls, hs, den0, den1, useGPU, dHs, dDen0, dDen1, dDen0cpy);
	    iteration++;
    }

    optParas[0] = h0;
    optParas[1] = alpha0;
    optParas[2] = Lmax;

	free(Ls);
    Ls = NULL;
}

///////// Guiming on 2016-03-16 ///////////////

// check whether the result from sequential computation and that from parallel computation agree
void CheckResults(AsciiRaster AsciiSEQ, AsciiRaster AsciiPARA){
	float eps = 0.000001f;

	int n = AsciiSEQ.nCols * AsciiSEQ.nRows;

	for(int i = 0; i < n; i++){
		if(abs(AsciiSEQ.elements[i] - AsciiPARA.elements[i]) > eps){
			printf("TEST FAILED. Result from parallel computation does not match that from sequential computation.\n");
			return;
		}
	}
	printf("TEST PASSED. Result from GPU computation does match that from CPU computation.\n");
}

float compGML(float* den0, int n){
	float gml = 0.0f;
	for(int i = 0; i < n; i++){
		gml = gml + log(den0[i]);
	}
	gml = expf(gml / n);
	return gml;
}

// reduction sum on GPU
void ReductionSumGPU(float* dArray, int numberOfElements){

   unsigned int N = numberOfElements;

   int iteration = 0;
   int NUM_ACTIVE_ITEMS = numberOfElements; // # active items need to be reduced

   // approx. # of blocks needed
   int NUM_BLOCKS = (numberOfElements ) / BLOCK_SIZE;

   // decide grid dimension
   int GRID_SIZE = (int)(sqrtf(NUM_BLOCKS)) + 1;
   dim3 dimGrid(GRID_SIZE, GRID_SIZE);

   // call the kernel for the first iteration
   ReductionSum<<<dimGrid, BLOCK_SIZE>>>(dArray, N, iteration, NUM_ACTIVE_ITEMS);

   // update # of items to be reduced in next iteration
   NUM_ACTIVE_ITEMS = (NUM_ACTIVE_ITEMS + BLOCK_SIZE - 1) / BLOCK_SIZE;

   // update numberOfElements (needed for deciding grid dimension)
   numberOfElements = dimGrid.x * dimGrid.y;

   // increment iteraton index
   iteration++;

   // iterate if needed
   while(numberOfElements > 1){
      NUM_BLOCKS = (numberOfElements ) / BLOCK_SIZE;

      GRID_SIZE = (int)(sqrtf(NUM_BLOCKS)) + 1;
      dimGrid.x = GRID_SIZE;
      dimGrid.y = GRID_SIZE;
      ReductionSum<<<dimGrid, BLOCK_SIZE>>>(dArray, N, iteration, NUM_ACTIVE_ITEMS);
      NUM_ACTIVE_ITEMS = (NUM_ACTIVE_ITEMS + BLOCK_SIZE - 1) / BLOCK_SIZE;

      numberOfElements = dimGrid.x * dimGrid.y;

      iteration++;
   }
}

// mark the boundary cells on a raster representing the study area
// By Guiming @ 2016-09-02
void MarkBoundary(AsciiRaster Ascii, bool useGPU){

	if(useGPU){ // do it on GPU
			// invoke kernels to mark the boundary of study area
			// execution config.
			int NBLOCK_W = (Ascii.nRows *  Ascii.nCols + BLOCK_SIZE - 1) / BLOCK_SIZE;
			int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
			dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);
			//printf("In Marking Boundary...\n");
			dMarkBoundary <<<dimGrid_W, BLOCK_SIZE / 2 >> > (Ascii);
	}
	else{ // do it on CPU
		for(int row = 0; row < Ascii.nRows; row++){
			for(int col = 0; col < Ascii.nCols; col++){

				if(Ascii.elements[row * Ascii.nCols + col] == Ascii.noDataValue)
					continue;

				if(row == 0 || (row == Ascii.nRows - 1) || col == 0 || (col == Ascii.nCols - 1)){ // cells in the outmost rows and cols
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}

				if(Ascii.elements[(row - 1) * Ascii.nCols + col - 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				if(Ascii.elements[row * Ascii.nCols + col - 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				if(Ascii.elements[(row + 1) * Ascii.nCols + col - 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}

				if(Ascii.elements[(row - 1) * Ascii.nCols + col] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				if(Ascii.elements[(row + 1) * Ascii.nCols + col] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}

				if(Ascii.elements[(row - 1) * Ascii.nCols + col + 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				if(Ascii.elements[row * Ascii.nCols + col + 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				if(Ascii.elements[(row + 1) * Ascii.nCols + col + 1] == Ascii.noDataValue){
					Ascii.elements[row * Ascii.nCols + col] = 1.0f;
					continue;
				}
				Ascii.elements[row * Ascii.nCols + col] = 0.0f;
			}
		}
	}
}

// compute the closest distances from sample points to study area boundary
// By Guiming @ 2016-09-02
void CalcDist2Boundary(SamplePoints Points, AsciiRaster Ascii, bool useGPU){

	// mark the boundary first
	//MarkBoundary(Ascii, useGPU); // either on GPU or CPU

	//printf("Done Marking Boundary!\n");

	if(useGPU){ // do it on GPU
			// invoke kernels to compute the nearest distance to boundary (squared) at each point
			// execution config.
			int NBLOCK_W = (Points.numberOfPoints + BLOCK_SIZE - 1) / BLOCK_SIZE;
	    int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
	    dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);

			dCalcDist2Boundary<<<dimGrid_W, BLOCK_SIZE>>>(Points, Ascii);
	}
	else{
		float p_x, p_y, cell_x, cell_y;
		for (int p = 0; p < Points.numberOfPoints; p++){
			float minDist = FLOAT_MAX;
			p_x = Points.xCoordinates[p];
			p_y = Points.yCoordinates[p];

			for (int row = 0; row < Ascii.nRows; row++){
				for (int col = 0; col < Ascii.nCols; col++){
					if (Ascii.elements[row*Ascii.nCols+col] == 1.0f){ // cells on boundary
						cell_x = COL_TO_XCOORD(col, Ascii.xLLCorner, Ascii.cellSize);
						cell_y = ROW_TO_YCOORD(row, Ascii.nRows, Ascii.yLLCorner, Ascii.cellSize);
						float d2 = Distance2(p_x, p_y, cell_x, cell_y);

						if(d2 < minDist){
							minDist = d2;
						}
					}
				}
			}

			Points.distances[p] = minDist;
			//printf("p: %d Points.distances[p]: %f minDist: %f\n", p, Points.distances[p]);
		}
	}
}

// By Guiming @ 2016-09-04
SamplePoints CopySamplePoints(const SamplePoints anotherPoints){ // copy points
	int n = anotherPoints.numberOfPoints;
	SamplePoints Points = AllocateSamplePoints(n);
	Points.numberOfPoints = n;
	Points.totalNumPoints = n;
	for(int p = 0; p < n; p++){
		Points.xCoordinates[p] = anotherPoints.xCoordinates[p];
		Points.yCoordinates[p] = anotherPoints.yCoordinates[p];
		Points.weights[p] = anotherPoints.weights[p];
		Points.distances[p] = anotherPoints.distances[p];
	}
	return Points;
}

// comparison function for sort
// By Guiming @ 2016-09-04
int compare ( const void *pa, const void *pb )
{
    const float *a = (const float *)pa;
    const float *b = (const float *)pb;
    if(a[0] == b[0])
        return a[1] - b[1];
    else
        return a[0] > b[0] ? 1 : -1;
}

void SortSamplePoints(SamplePoints Points) {
	const int n = 100;
	SamplePoints temPoints = CopySamplePoints(Points);

	float distances[n][2];
	for (int i = 0; i < n; i++)
	{
		distances[i][0] = Points.distances[i];
		distances[i][1] = i * 1.0f;
	}
	/*
	for(int i = 0; i < n; ++i)
	  printf("%.1f\n", Points.distances[i]);
  printf("\n");
	*/

	qsort(distances, n, sizeof(distances[0]), compare);

	for (int i = 0; i < n; i++)
	{
		int idx = (int)distances[i][1];
		Points.xCoordinates[i] = temPoints.xCoordinates[idx];
		Points.yCoordinates[i] = temPoints.yCoordinates[idx];
		Points.weights[i] = temPoints.weights[idx];
		Points.distances[i] = temPoints.distances[idx];
	}
	/*
	for(int i = 0; i < n; ++i)
	  printf("%.1f\n", Points.distances[i]);
	*/
	FreeSamplePoints(&temPoints);

}

// build a KDtree on sample points
// By Guiming @ 2016-09-07
void BuildCPUKDtree (SamplePoints Points){
	int NPTS = Points.numberOfPoints;
	dataP = vector<Point>(NPTS);
	for(int i = 0; i < NPTS; i++){
		dataP[i].coords[0] = Points.xCoordinates[i];
    dataP[i].coords[1] = Points.yCoordinates[i];
	}
	int max_level = (int)(log(dataP.size())/log(2) / 2) + 1;
	tree.Create(dataP, max_level);
}

void BuildGPUKDtree ()
{
	for(int i=0;i<GPU_N;i++)
	{
		cudaSetDevice(i);
		GPU_tree[i].CreateKDTree(tree.GetRoot(), tree.GetNumNodes(), dataP);	
	}
	cudaSetDevice(0);
}

//Enable P2P Access Across Devices
//Timothy @ 08/13/2020
void EnableP2P()
{
	cudaError_t error = cudaSuccess;
	for (int id = 0; id < GPU_N; ++id)
	{
		cudaSetDevice(id);
		const int top = id > 0 ? id - 1 : (GPU_N - 1); //Int representing first in list of GPUs
		int capable = 1; //(T/F) P2P Access is enabled between devices 
		error = cudaDeviceCanAccessPeer(&capable, id, top);
		if (error != cudaSuccess)
		{
			printf("ERROR 1 in EnableP2P: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		if (capable)
		{
			printf("Enabled P2P for Device %d...\n", id);
			cudaDeviceEnablePeerAccess(top, 0);
		}
		else if (!capable){printf("NOT CAPABLE! P2P for Device %d...\n", id);}
		const int bottom = (id + 1) % GPU_N;
		if (top != bottom)
		{
			error = cudaDeviceCanAccessPeer(&capable, id, bottom);
			if (error != cudaSuccess)
			{
				printf("ERROR 2 in EnableP2P: %s\n", cudaGetErrorString(error));
				exit(EXIT_FAILURE);
			}
			if (capable)
			{
				printf("Enabling P2P for Device %d...\n", id);
				cudaDeviceEnablePeerAccess(bottom, 0);
			}
			else if (!capable){printf("NOT CAPABLE! P2P for Device %d...\n", id);}
		}
	}
	cudaSetDevice(0); //Reset device to first GPU
}

//By Timothy @ 08/14/2020
//Determine next Device to be used based on passed integers assumed to represent their numbers
void nextDev(int numDev, int& curDev)
{
	if (curDev == (numDev - 1))
	{
		curDev = 0;
	}
	else
	{
		curDev++;
	}
}

//Timothy @ 08/24/2020
//Function to check device properties, primarily for troubleshooting purposes
void DevProp()
{
	for (int i = 0; i < GPU_N; i++)
	{
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n", prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n", prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n\n", 2.0 * prop.memoryClockRate * (prop.memoryBusWidth / 8) / 1.0e6);
	}
}

////Timothy @ 12/29/20
////Function which copies each group of points into a temorary place on the host, before copying their values to
////hPoints in order to reform the original group
//void ReformPoints(SamplePoints hPoints, SamplePoints* dPoints)
//{
//	int n = hPoints.numberOfPoints; //Number of TOTAL points
//	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
//	int div = n / GPU_N; //Division of points to be divided amongst GPUs
//	int size; //Size of data chunk being copied to tempPoints
//	int index = 0; //Index for the points we are reforming into
//
//	size = div * sizeof(float);
//
//	cudaError_t error = cudaSuccess;
//
//	SamplePoints tempPoints = AllocateSamplePoints((div + rem));
//	if (error != cudaSuccess)
//	{
//		printf("ERROR 0 in ReformPoints: %s\n", cudaGetErrorString(error));
//		exit(EXIT_FAILURE);
//	}
//
//	for (int device = 0; device < GPU_N; device++)
//	{
//		cudaSetDevice(device);
//		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
//		//if not add remainder to size on final GPU
//		if ((device == GPU_N - 1) && (rem != 0))
//		{
//			div += rem;
//		}
//
//		int NBLOCK_W = (n + BLOCK_SIZE - 1) / BLOCK_SIZE;
//		int GRID_SIZE_W = (int)(sqrtf(NBLOCK_W)) + 1;
//		dim3 dimGrid_W(GRID_SIZE_W, GRID_SIZE_W);
//		cudaSetDevice(device);
//		printf("%d:Points...\n", device); //DEBUGGING
//		PrintPoints << <dimGrid_W, BLOCK_SIZE, 0, streams[device] >> > (dPoints[device], 100);
//		cudaStreamSynchronize(streams[device]);
//	
//		//Copy all data from chunk to tempPoints
//		error = cudaMemcpy(tempPoints.xCoordinates, dPoints[device].xCoordinates, size, cudaMemcpyDeviceToHost);
//		if (error != cudaSuccess)
//		{
//			printf("ERROR 1 in ReformPoints: %s\n", cudaGetErrorString(error));
//			exit(EXIT_FAILURE);
//		}
//		error = cudaMemcpy(tempPoints.xCoordinates, dPoints[device].xCoordinates, size, cudaMemcpyDeviceToHost);
//		if (error != cudaSuccess)
//		{
//			printf("ERROR 2 in ReformPoints: %s\n", cudaGetErrorString(error));
//			exit(EXIT_FAILURE);
//		}
//		error = cudaMemcpy(tempPoints.xCoordinates, dPoints[device].xCoordinates, size, cudaMemcpyDeviceToHost);
//		if (error != cudaSuccess)
//		{
//			printf("ERROR 3 in ReformPoints: %s\n", cudaGetErrorString(error));
//			exit(EXIT_FAILURE);
//		}
//		// By Guiming @ 2016-09-02
//		error = cudaMemcpy(tempPoints.xCoordinates, dPoints[device].xCoordinates, size, cudaMemcpyDeviceToHost);
//		if (error != cudaSuccess)
//		{
//			printf("ERROR 4 in ReformPoints: %s\n", cudaGetErrorString(error));
//			exit(EXIT_FAILURE);
//		}
//
//		//Loop to merge copied chunk of points into hPoints
//		for (int i = 0; i < hPoints.numberOfPoints; i++)
//		{
//			hPoints.xCoordinates[index] = tempPoints.xCoordinates[index];
//			hPoints.yCoordinates[index] = tempPoints.yCoordinates[index];
//			hPoints.weights[index] = tempPoints.weights[index];
//			hPoints.distances[index] = tempPoints.distances[index];
//			index++;
//		}
//	}
//	cudaSetDevice(0); //Reset device to first GPU
//	//Free temp points
//	FreeSamplePoints(&tempPoints);
//}

//Timothy @ 08/13/2021
//Reforms points using indeces rather than actually changing any memory
//We realized when reforming other data values that each GPU already has copies of the full set of point structs
void ReformPoints(SamplePoints* dPoints)
{
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		dPoints[i].numberOfPoints = sPoints.numberOfPoints;
		dPoints[i].start = 0;
		dPoints[i].end = 100;
	}
	cudaSetDevice(0);
}

//Timothy @ 08/13/2021
//Divides points using indeces rather than actually changing any memory
void DividePoints(SamplePoints* dPoints)
{
	int n = sPoints.numberOfPoints; //Number of TOTAL points
	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
	int div = n / GPU_N; //Division of points to be divided amongst GPUs
	int index = 0; //Index to track start of each data division
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
		//if not add remainder to size on final GPU
		if ((i == GPU_N - 1) && (rem != 0))
		{
			div += rem;
		}
		dPoints[i].numberOfPoints = div;
		dPoints[i].start = index; //Begin tracking division of points
		index += div; //Add division size to index
		dPoints[i].end = index; //Tracking end of division
	}
	cudaSetDevice(0);
}

//Timothy @ 08/10/2021
//Reform density arrays on host and copy back accross devices
void ReformDensities(float** dDen, float* hDen)
{
	int n = sPoints.numberOfPoints; //Number of TOTAL points
	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
	int div = n / GPU_N; //Division of points to be divided amongst GPUs
	int size; //Size of data chunk being copied to tempPoints
	int index = 0; //Index for the points we are reforming into
	cudaError_t error = cudaSuccess;
	float* tempDen = (float*)malloc(n * sizeof(float));
	size = n * sizeof(float);
	
	for (int device = 0; device < GPU_N; device++)
	{
		cudaSetDevice(device);
		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
		//if not add remainder to size on final GPU
		if ((device == GPU_N - 1) && (rem != 0))
		{
			div += rem;
		}
		
		//Copy all data from chunk to tempPoints
		error = cudaMemcpy(tempDen, dDen[device], size, cudaMemcpyDeviceToHost);
		if (error != cudaSuccess)
		{
			printf("ERROR 1.%d in ReformDensities: %s\n", device, cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		//Loop to merge copied chunk of points into hPoints
		for (int i = 0; i < div; i++)
		{
			hDen[index] = tempDen[index];
			index++;
		}
	}
	//Copy reformed dDen accross GPUs
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dDen[i], hDen, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 2 in ReformDensities: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	//Cleanup
	cudaSetDevice(0); //Reset device to first GPU
	//Free temp points
	FreeDen(tempDen);
}

//Timothy @ 08/13/2021
//Reform bandwidth arrays on host and copy back accross devices
void ReformBandwidths(float** dBand, float* hBand) 
{
	int n = sPoints.numberOfPoints; //Number of TOTAL points
	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
	int div = n / GPU_N; //Division of points to be divided amongst GPUs
	int size; //Size of data chunk being copied to tempPoints
	int index = 0; //Index for the points we are reforming into
	cudaError_t error = cudaSuccess;
	float* tempBand = (float*)malloc(n * sizeof(float));
	size = n * sizeof(float);

	for (int device = 0; device < GPU_N; device++)
	{
		cudaSetDevice(device);
		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
		//if not add remainder to size on final GPU
		if ((device == GPU_N - 1) && (rem != 0))
		{
			div += rem;
		}

		//Copy all data from chunk to tempPoints
		error = cudaMemcpy(tempBand, dBand[device], size, cudaMemcpyDeviceToHost);
		if (error != cudaSuccess)
		{
			printf("ERROR 1.%d in ReformBandwidths: %s\n", device, cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		//Loop to merge copied chunk of points into hPoints
		for (int i = 0; i < div; i++)
		{
			hBand[index] = tempBand[index];
			index++;
		}
	}
	//Copy reformed dDen accross GPUs
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dBand[i], hBand, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 2 in ReformBandwidths: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	//Cleanup
	cudaSetDevice(0); //Reset device to first GPU
	//Free temp bands
	FreeBandwidths(tempBand);
}

//Timothy @ 08/13/2021
//Reform EC Weight arrays on host and copy back accross devices
void ReformECWeights(float** dWeights, float* hWeights)
{
	int n = sPoints.numberOfPoints; //Number of TOTAL points
	int rem = n % GPU_N; //Remainder to determine if number of GPUs divides Number of Points evenly
	int div = n / GPU_N; //Division of points to be divided amongst GPUs
	int size; //Size of data chunk being copied to tempPoints
	int index = 0; //Index for the points we are reforming into
	cudaError_t error = cudaSuccess;
	float* tempWeights = (float*)malloc(n * sizeof(float));
	size = n * sizeof(float);

	for (int device = 0; device < GPU_N; device++)
	{
		cudaSetDevice(device);
		//If on last GPU, check if GPU_N divided into points evenly (rem==0) 
		//if not add remainder to size on final GPU
		if ((device == GPU_N - 1) && (rem != 0))
		{
			div += rem;
		}

		//Copy all data from chunk to tempPoints
		error = cudaMemcpy(tempWeights, dWeights[device], size, cudaMemcpyDeviceToHost);
		if (error != cudaSuccess)
		{
			printf("ERROR 1.%d in ReformBandwidths: %s\n", device, cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
		//Loop to merge copied chunk of points into hPoints
		for (int i = 0; i < div; i++)
		{
			hWeights[index] = tempWeights[index];
			index++;
		}
	}
	//Copy reformed dDen accross GPUs
	for (int i = 0; i < GPU_N; i++)
	{
		cudaSetDevice(i);
		error = cudaMemcpy(dWeights[i], hWeights, size, cudaMemcpyHostToDevice);
		if (error != cudaSuccess)
		{
			printf("ERROR 2 in ReformBandwidths: %s\n", cudaGetErrorString(error));
			exit(EXIT_FAILURE);
		}
	}
	//Cleanup
	cudaSetDevice(0); //Reset device to first GPU
	//Free temp weights
	FreeEdgeCorrectionWeights(tempWeights);
}