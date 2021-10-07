// Copyright 2016 Guiming Zhang (gzhang45@wisc.edu)
// Distributed under GNU General Public License (GPL) license

#ifndef _SAMPLEPOINTS_H_
#define _SAMPLEPOINTS_H_


// SamplePoints Structure declaration
typedef struct {
	unsigned int totalNumPoints; //Used to keep track of the entire number of points for reforming purposes
	unsigned int numberOfPoints; //Dynamic number of points which helps track number of points being handled by each GPU
	unsigned int start; //Used to identify starting index when points are being worked with accross multiple GPUs
	unsigned int end; //Ending index
	float* xCoordinates;
	float* yCoordinates;
	float* weights;
	float* distances; // closest distances (squared) to study area boundary (by Guiming @ 2016-09-02)
} SamplePoints;


#endif // _SAMPLEPOINTS_H_
