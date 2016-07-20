#include "dataStructures.h"
#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <vector>
#include <numeric>
#include <climits>
#include <thread>

#include <stdlib.h> // DEBUG
using namespace std;

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

namespace SpkDslowFilter {
class InterpDetection {	

public:
	InterpDetection(int cols, int rows, double samplingRate);
	~InterpDetection();
	void detect(unsigned short* vm, int t0, int t1, int tCut);

private:
	inline int interpolateFourChannels(int* V, int ch);
	int* computeFourChInterp(int* V, int start, int tInc);
	inline int interpolateFiveChannels(int* V, int ch);
	int* computeFiveChInterp(int* V, int start, int tInc);
	inline bool isOutlier(unsigned short v);
	int* computeMedian(unsigned short* vm, int tInc);
	inline void updateBaseline(int v, int ch);
	void initialiseVMovingAvg(unsigned short* vm, int* vGlobal);
	void initialiseVGlobalMovingAvg(int* vGlobal);	
	int* preprocessData(unsigned short* vm, int* vGlobal, int start, int tInc);
	void findSpikes(int* fourChInterp, int* fiveChInterp, int start, int t0, int tInc);

	// Parallel functions
	void computeFourChInterpThread(int threadID, int* fourChInterp, int* V, 
		    int start, int tInc);
	void computeFiveChInterpThread(int threadID, int* fiveChInterp, int* V, 
			int start, int tInc);
	void preprocessDataThread(int threadID, unsigned short* vm, int* vGlobal, 
			int start, int tInc, int* Qdiff, int* Qmax, int* vGlobalMovingAvg);

	bool detectionInitialised;
	
	int chCols;
	int chRows;
	int NChannels;

	int movingWindowLength;

	CheckSpace registry;
	int chunkSize;
	int* baseline; 		// Size NChannels
	int* variability; 	// Size: NChannels
	int* vMovingAvg; // Size: NChannels
	int vGlobalMovingAvgInterIt;
	int* outlierWait; // Size: NChannels
	int outlierMark = INT_MIN;
	int outlierWaitingTime = 10; // Frames
	unsigned short scale; // Scale all data to increase the resolution

	int nSpikes;

	// Algorithm parameters
	int initialBaseline;	
	int initialVariability;
	int minVariability;
	unsigned short d_pitch; // Electrode pitch [microV]
	int maxSpikeDelay; // Interval for depolarisation area
	double f_s; // Sampling rate
	double f_v; // Variability update rate
	double f_b; // Baseline update rate
	float theta; // Detection threshold
	float theta_b; // Repolarisation threshold
	float theta_ev; // Minimum depolarisation area
	int tau_event; // Characteristic event length
	int tau_coinc; // Window for coincident events
	int w_cs; // Center/surround weighting (5-channel interpolation) 
	float tau_pre; // Cut-out window before peak
	float tau_post; // Cut-out window after peak

	unsigned short max_out_threshold; // Outlier threshold (max)
	unsigned short min_out_threshold; // 

	int DEBUG_CH = 2310;

	int nthreads;
	std::thread* threads;
};
};