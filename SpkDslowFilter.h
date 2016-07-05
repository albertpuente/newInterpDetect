#include <algorithm> // sort, max, min
#include <climits>

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

	int chCols;
	int chRows;
	int NChannels;
	unsigned short w_cs = 4; // Center/surround weighting (5-channel interpolation) 
};
};