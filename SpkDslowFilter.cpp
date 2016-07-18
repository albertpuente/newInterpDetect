#include "SpkDslowFilter.h"

namespace SpkDslowFilter {

InterpDetection::InterpDetection(int rows, int cols, double samplingRate) {
    chRows = rows;
    chCols = cols;
    NChannels = rows*cols;

    scale = 64;

	theta = 5 * scale; // Detection threshold
    w_cs = 4; // Center/surround weighting (5-channel interpolation) 

    /**
	theta_b = 0; // Repolarisation threshold
	theta_ev; // Minimum depolarisation area
	tau_ev; // Interval for depolarisation area
	tau_event = 7; // (1 ms) Characteristic event length
	tau_coinc = 2; // (0.27 ms) Window for coincident events
	 
	tau_pre = 7; // (1 ms) Cut-out window before peak
	tau_post = 15; // (2.2 ms) Cut-out window after peak
    */

    f_s = samplingRate;


    f_v = 4; // ??? 
    f_b = 32; // ??? 

    movingWindowLength = 5; // In frames

    initialBaseline = 0;    
    minVariability = scale*0.5;
    initialVariability = minVariability; // ??? Provisional

    variability = new int[NChannels];
    baseline = new int[NChannels];
    std::fill_n(variability, NChannels, initialVariability);
    std::fill_n(baseline, NChannels, initialBaseline);

    max_out_threshold = 4000;
    min_out_threshold = 10;

    detectionInitialised = false;
}

inline int InterpDetection::interpolateFourChannels(unsigned short* vm, int ch) {
    // Average the three largest absolute signal amplitudes
    // A B -> E ·
    // C D    · ·
    // Where E = sum(A,B,C,D) - max(A,B,C,D)
    unsigned short values[] = {vm[ch],          vm[ch + 1],
                               vm[ch + chCols], vm[ch + chCols + 1]};        
    unsigned short interp = 0;
    unsigned short maxValue = 0;
    for (auto v : values) { // Add all the values and find the minimum
        interp += v;
        maxValue = max(maxValue, v);        
    }
    return ((int) interp - (int) maxValue)*scale/3;
}

int* InterpDetection::computeFourChInterp(unsigned short* vm, int* vGlobal, long tInc) {
    auto fourChInterp = new int[NChannels*tInc]; 
    for (int t = 0; t < tInc; t++) { // Parallelisable
        // Not all elements in fourChInterp will be used (last column and row)
        for (int i = 0; i < chRows - 1; i++) {
            for (int j = 0; j < chCols - 1; j++) { 
                int ch = i*chRows + j;                          
                fourChInterp[ch + t*NChannels] = interpolateFourChannels(vm, ch) - vGlobal[t];
            } 
        }
    } 
    return fourChInterp;  
}

inline int InterpDetection::interpolateFiveChannels(unsigned short* vm, int ch) {
    // Average using w_cs/(3 + w_cs)*the center and 1/(3 + w_cs)* the three largest
    // surrounding channels.
    // · A ·    · · ·
    // D E B -> · F ·
    // · C ·    · · ·
    // Where F = E*w_cs/(3 + w_cs) + (sum(A..D) - max(A..D))/(3 + w_cs)
    // (This is actually computed in another equivalent way).
    unsigned short values[] = {vm[ch - chCols], vm[ch + 1], 
                               vm[ch + chCols], vm[ch - 1]};        
    int interp = (scale*(int) vm[ch]*w_cs)/(3 + w_cs);
    int maxValue = INT_MIN;
    for (auto v : values) { // Add all the values and find the minimum
        int weightedValue = (scale*(int) v)/(3 + w_cs);
        interp += weightedValue;
        maxValue = max(maxValue, weightedValue);        
    }
    return interp - maxValue;
}

int* InterpDetection::computeFiveChInterp(unsigned short* vm, int* vGlobal, long tInc) {
    auto fiveChInterp = new int[NChannels*tInc]; // Has more elements than needed
    for (int t = 0; t < tInc; t++) { // Parallelisable
        // Not all elements in fiveChInterp will be used (first/last columns and rows)
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                int ch = i*chRows + j;
                
                fiveChInterp[ch + t*NChannels] = interpolateFiveChannels(vm, ch) - vGlobal[t];
            }
        }        
    }
    return fiveChInterp;  
}

inline bool InterpDetection::isOutlier(unsigned short v) {
    return v >= max_out_threshold || v <= min_out_threshold;
} 

int* InterpDetection::computeMedian(unsigned short* vm, long tInc) {
    auto vGlobal = new int[tInc];
    for (int t = 0; t < tInc; t++) {
        // Take the mean value for each time-step
        int n = 0;
        long sum = 0;
        for (int ch = 0; ch < NChannels; ch++) {
            if (not isOutlier(vm[ch + t*NChannels])) {
                sum += scale*(int) vm[ch + t*NChannels];
                n++;
            } 
        }
        vGlobal[t] = sum/n;
    }
    return vGlobal;
}

void InterpDetection::updateBaseline(int* fiveChInterp, int t) {
    for (int i = 1; i < chRows - 1; i++) {
        for (int j = 1; j < chCols - 1; j++) {
            int ch = i*chRows + j;
            int v = fiveChInterp[ch + t*NChannels];
            // Baseline
            if (v > baseline[ch] + variability[ch]) {
                baseline[ch] += f_b/2;
            }
            else if (v < baseline[ch] - variability[ch]) {
                baseline[ch] -= f_b;
            }
            
            // Variability
            // v in (baseline - variability, baseline] or (-inf, baseline - 6variability])
            if ((v > baseline[ch] - variability[ch] && v <= baseline[ch]) ||
                (v > INT_MIN && v < baseline[ch] - 6*variability[ch]) ) {

                variability[ch] -= f_v;
                variability[ch] = max(variability[ch], minVariability);
            } // v in (baseline - 5variability)
            else if (v > baseline[ch] - 5 *variability[ch] && 
                     v <=  baseline[ch] - variability[ch]) {
                variability[ch] += f_v;
            }

        }
    }
}

void InterpDetection::findSpike(int* fourChInterp, int* fiveChInterp, 
                                unsigned short* vm, int i, int j, int t) {
    int ch = i*chRows + j;
    int v = fiveChInterp[ch + t*NChannels];

    if (v < baseline[ch] - theta*variability[ch]) { // Spike condition
        // cout << "SPIKE Ch:" << ch << " t:" << t << " v:" << v << endl << flush;
        // bool relevant =  registry.collides(i, j);            
        // registry.addSpike(t, i, j, relevant);
    }
}

void InterpDetection::initialiseVMovingAvg(unsigned short* vm) {
    vMovingAvg = new int[NChannels];
    std::fill_n(vMovingAvg, NChannels, 0);
    for (int ch = 0; ch < NChannels; ++ch) {
        for (int t = 0; t < movingWindowLength; ++t) {
            vMovingAvg[ch] += (vm[ch + t*NChannels]*scale)/movingWindowLength;
        }
    }
}

void InterpDetection::initialiseMedianMovingAvg(int* vGlobal) {
    medianMovingAvg = 0;
    for (int t = 0; t < movingWindowLength; ++t) {
        medianMovingAvg += vGlobal[t]/movingWindowLength;
    }
}

/**
int* InterpDetection::computeMedianDerivative(int* vGlobal, int tInc) {
    auto vGlobalDerivative = new int[tInc];
    std::adjacent_difference(vGlobal, vGlobal + tInc, vGlobalDerivative);
    vGlobalDerivative[0] = 0;
    return vGlobalDerivative;
}
*/

void InterpDetection::detect(unsigned short* vm, long t0, long t1) {    
    // vm is indexed as follows:
    //     vm[channel + tInc*nChannels]
    // or
    //     vm[i*chRows + j + tInc*nChannels]

    // Duration of the interval
    int tInc = t1 - t0;

    cout << "Computing vGlobal ...\n" << flush;
    auto vGlobal = computeMedian(vm, tInc);
    cout << "Computing interpolations ...\n" << flush;
    // The last row and column are empty.    
    auto fourChInterp = computeFourChInterp(vm, vGlobal, tInc);

    // The first/last rows and colums are empty.
    auto fiveChInterp = computeFiveChInterp(vm, vGlobal, tInc);

    int startFrame = 0;

    if (not detectionInitialised) {
        initialiseVMovingAvg(vm);
        initialiseMedianMovingAvg(vGlobal);
        startFrame += movingWindowLength;
        detectionInitialised = true;
    }

    // min distance (channels), min time (frames), rows, columns
    chunkSize = 8; // 64 chunks (8x8 chunks of 8x8 channels each)
    registry.initialise(5, 20, chRows, chCols, chunkSize);
    // Alternate sweep
    for (int t = startFrame; t < tInc; t++) {
        updateBaseline(fiveChInterp, t);

        // COUT ====================
        cout << baseline[DEBUG_CH] << " " << variability[DEBUG_CH] << " " << fiveChInterp[DEBUG_CH + t*NChannels] << " " << ((int)vm[DEBUG_CH + t*NChannels])*scale - vGlobal[t]<< endl;
        //==========================
        for (int i = 1; i < chRows - 1; i++) { 
            for (int j = 1; j < chCols - 1; j++) {                
                findSpike(fourChInterp, fiveChInterp, vm, i, j, t);
            }
        }
        registry.pruneOldSpikes(t); 
    }  
}

}
