#include "dataStructures.h"
using namespace std;

namespace SpkDslowFilter {
class InterpDetection {	

public:
	InterpDetection(int cols, int rows);
	~InterpDetection();
	void detect(unsigned short* vm, long t0, long t1);

private:
	unsigned short interpolateFourChannels(unsigned short* vm, int ch);
	unsigned short* computeFourChInterp(unsigned short* vm, long tInc);
	unsigned short interpolateFiveChannels(unsigned short* vm, int ch);
	unsigned short* computeFiveChInterp(unsigned short* vm, long tInc);
	unsigned short* computeBaseline(unsigned short* vm, long tInc);
	inline bool isOutlier(unsigned short vm);
	void findSpike(unsigned short* fourChInterp, unsigned short* fiveChInterp, 
				   int i, int j, int t);

	CheckSpace registry;
	unsigned short* x_global;
	unsigned short* fourChInterp;
	unsigned short* fiveChInterp;

	int chCols;
	int chRows;
	int NChannels;

	// Algorithm parameters
	float f_s; // Sampling rate
	unsigned short d_pitch; // Electrode pitch [microV]
	float f_v; // Variability update rate
	float f_b; // Baseline update rate
	float theta; // Detection threshold
	float theta_b; // Repolarisation threshold
	float theta_ev; // Minimum depolarisation area
	float tau_ev; // Interval for depolarisation area
	float tau_event; // Characteristic event length
	float tau_coinc; // Window for coincident events
	unsigned short w_cs; // Center/surround weighting (5-channel interpolation) 
	float tau_pre; // Cut-out window before peak
	float tau_post; // Cut-out window after peak

	unsigned short max_out_threshold; // Outlier threshold (max)
};
};