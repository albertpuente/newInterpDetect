#include "SpkDslowFilter.h"

namespace SpkDslowFilter {

InterpDetection::InterpDetection(int rows, int cols) {
    chRows = rows;
    chCols = cols;
    NChannels = rows*cols;

	theta = 5; // Detection threshold
	theta_b = 0; // Repolarisation threshold
	theta_ev; // Minimum depolarisation area
	tau_ev; // Interval for depolarisation area
	tau_event = 7; // (1 ms) Characteristic event length
	tau_coinc = 2; // (0.27 ms) Window for coincident events
	w_cs = 4; // Center/surround weighting (5-channel interpolation) 
	tau_pre = 7; // (1 ms) Cut-out window before peak
	tau_post = 15; // (2.2 ms) Cut-out window after peak
}

unsigned short InterpDetection::interpolateFourChannels(unsigned short* vm, int ch) {
    // Average the three largest absolute signal amplitudes
    // A B -> E ·
    // C D    · ·
    // Where E = sum(A,B,C,D) - min(A,B,C,D)
    unsigned short values[] = {vm[ch],          vm[ch + 1],
                               vm[ch + chCols], vm[ch + chCols + 1]};        
    unsigned short interp = 0;
    unsigned short minValue = USHRT_MAX;
    for (auto v : values) { // Add all the values and find the minimum
        interp += v;
        minValue = min(minValue, v);        
    }
    return (interp - minValue)/3;
}

unsigned short* InterpDetection::computeFourChInterp(unsigned short* vm, long tInc) {
    auto fourChInterp = new unsigned short[NChannels*tInc]; 
    for (int t = 0; t < tInc; tInc++) {
        // Not all elements in fourChInterp will be used (last column and row)
        for (int i = 0; i < chRows - 1; i++) {
            for (int j = 0; j < chCols - 1; j++) {
                int ch = i*chRows + j;                
                fourChInterp[ch + t*NChannels] = interpolateFourChannels(vm, ch);
            }
        }        
    }
    return fourChInterp;  
}

unsigned short InterpDetection::interpolateFiveChannels(unsigned short* vm, int ch) {
    // Average using w_cs/(3 + w_cs)*the center and 1/(3 + w_cs)* the three largest
    // surrounding channels.
    // · A ·    · · ·
    // D E B -> · F ·
    // · C ·    · · ·
    // Where F = E*w_cs/(3 + w_cs) + (sum(A..D) - min(A..D))/(3 + w_cs)
    // (This is actually computed in another equivalent way).
    unsigned short values[] = {vm[ch - chCols], vm[ch + 1], 
                               vm[ch + chCols], vm[ch - 1]};        
    unsigned short interp = vm[ch]*w_cs/(3 + w_cs);
    unsigned short minValue = USHRT_MAX;
    for (auto v : values) { // Add all the values and find the minimum
        unsigned short weightedValue = v/(3 + w_cs);
        interp += weightedValue;
        minValue = min(minValue, weightedValue);        
    }
    return interp - minValue;
}

unsigned short* InterpDetection::computeFiveChInterp(unsigned short* vm, long tInc) {
    auto fiveChInterp = new unsigned short[NChannels*tInc]; // Has more elements than needed
    for (int t = 0; t < tInc; tInc++) {
        // Not all elements in fiveChInterp will be used (first/last columns and rows)
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                int ch = i*chRows + j;                
                fiveChInterp[ch + t*NChannels] = interpolateFiveChannels(vm, ch);
            }
        }        
    }
    return fiveChInterp;  
}

inline bool InterpDetection::isOutlier(unsigned short vm) {
    return vm >= max_out_threshold;
} 

unsigned short* InterpDetection::computeBaseline(unsigned short* vm, long tInc) {
    auto x_global = new unsigned short[tInc];
    for (int t = 0; t < tInc; tInc++) {
        // Take the mean value for each time-step
        int n = 0;
        long sum = 0;
        for (int ch = 0; ch < NChannels; ch++) {
            if (not isOutlier(vm[ch + t*NChannels])) {
                sum += vm[ch + t*NChannels];
                n++;
            } 
        }
        x_global[t] = sum/n;
    }
    return x_global;
}

void InterpDetection::findSpike(unsigned short* fourChInterp, 
                                unsigned short* fiveChInterp, 
                                int i, int j, int t) {
    int ch = i*chRows + j;
    unsigned short v = fiveChInterp[ch + t*NChannels];

    if (false) { // Spike condition
        if (not registry.collides(i, j)) {            
            registry.addSpike(t, i, j);
        }
    }
}

void InterpDetection::detect(unsigned short* vm, long t0, long t1) {    
    // vm is indexed as follows:
    //     vm[channel + tInc*nChannels]
    // or
    //     vm[i*chRows + j + tInc*nChannels]

    // Duration of the interval
    int tInc = t1 - t0;

    x_global = computeBaseline(vm, tInc);

    // The last row and column are empty.    
    fourChInterp = computeFourChInterp(vm, tInc);

    // The first/last rows and colums are empty.
    fiveChInterp = computeFiveChInterp(vm, tInc);

    // min distance (channels), min time (frames), rows, columns
    chunkSize = 8; // 64 chunks (8x8 chunks of 8x8 channels each)
    registry.initialise(5, 20, chRows, chCols, chunkSize);

    // Alternate sweep
    for (int t = 0; t < tInc; t++) {
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                findSpike(fourChInterp, fiveChInterp, i, j, t);
            }
        }
        registry.pruneOldSpikes(t);
    }  
}

}
