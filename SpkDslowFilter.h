#include "dataStructures.h"
#include <iostream>
#include <numeric>

using namespace std;

namespace SpkDslowFilter {
class InterpDetection {	

public:
	InterpDetection(int cols, int rows, double samplingRate);
	~InterpDetection();
	void detect(unsigned short* vm, long t0, long t1);

private:
	inline int interpolateFourChannels(unsigned short* vm, int ch);
	int* computeFourChInterp(unsigned short* vm, int* vGlobal, long tInc);
	inline int interpolateFiveChannels(unsigned short* vm, int ch);
	int* computeFiveChInterp(unsigned short* vm, int* vGlobal, long tInc);
	int* computeMedian(unsigned short* vm, long tInc);
	inline bool isOutlier(unsigned short v);
	void updateBaseline(int* fiveChInterp, int t);
	void findSpike(int* fourChInterp, int* fiveChInterp, unsigned short* vm,
				   int i, int j, int t);
	
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
	int medianMovingAvg;

	unsigned short scale; // Scale all data to increase the resolution

	// Algorithm parameters
	int initialBaseline;	
	int initialVariability;
	int minVariability;
	unsigned short d_pitch; // Electrode pitch [microV]
	double f_s; // Sampling rate
	double f_v; // Variability update rate
	double f_b; // Baseline update rate
	float theta; // Detection threshold
	float theta_b; // Repolarisation threshold
	float theta_ev; // Minimum depolarisation area
	float tau_ev; // Interval for depolarisation area
	float tau_event; // Characteristic event length
	float tau_coinc; // Window for coincident events
	int w_cs; // Center/surround weighting (5-channel interpolation) 
	float tau_pre; // Cut-out window before peak
	float tau_post; // Cut-out window after peak

	unsigned short max_out_threshold; // Outlier threshold (max)
	unsigned short min_out_threshold; // 

	int DEBUG_CH = 2310;
};
};