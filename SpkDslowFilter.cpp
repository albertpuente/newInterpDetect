#include "SpkDslowFilter.h"

namespace SpkDslowFilter {

InterpDetection::InterpDetection(int rows, int cols, double samplingRate) {
    chRows = rows;
    chCols = cols;
    NChannels = rows*cols;

    scale = 64;

	theta =  - 3.12 * scale; // Detection threshold
    theta_b = 0; // Repolarisation threshold

    w_cs = 4; // Center/surround weighting (5-channel interpolation) 

	theta_ev; // Minimum depolarisation area
	
	tau_pre = 7; // (1 ms) Cut-out window before peak
	tau_post = 15; // (2.2 ms) Cut-out window after peak

    tau_coinc = 2; // (0.27 ms) Window for coincident events
    tau_event = 7; // (1 ms) Characteristic event length	
    f_s = samplingRate;

    maxSpikeDelay = f_s/1000;

    f_v = 4; // ??? 
    f_b = 32; // ??? 

    movingWindowLength = 5; // In frames
    
    initialBaseline = 0;    
    minVariability = scale*0.5;
    initialVariability = minVariability; // ??? Provisional

    variability = new int[NChannels];
    baseline = new int[NChannels];
    fill_n(variability, NChannels, initialVariability);
    fill_n(baseline, NChannels, initialBaseline);

    max_out_threshold = 4000;
    min_out_threshold = 10;

    detectionInitialised = false;

    outlierWait = new int[NChannels];

    chunkSize = 8; // Default: 64 chunks (8x8 chunks of 8x8 channels each)
}

inline int InterpDetection::interpolateFourChannels(int* V, int ch) {
    // Average the three largest absolute signal amplitudes
    // A B -> E ·
    // C D    · ·
    // Where E = sum(A,B,C,D) - max(A,B,C,D)
    int values[] = {V[ch],          V[ch + 1],
                    V[ch + chCols], V[ch + chCols + 1]};        
    int interp = 0;
    int maxValue = 0;
    bool outlier = false;
    for (auto v : values) { // Add all the values and find the minimum
        if (v == outlierMark) {
            outlier = true;
            break;
        }
        interp += v;
        maxValue = max(maxValue, v);        
    }
    if (outlier) {
        return outlierMark;
    }
    else {
        return (interp - maxValue)*scale/3;
    }
}

int* InterpDetection::computeFourChInterp(int* V, int start, int tInc) {
    auto fourChInterp = new int[NChannels*tInc]; 
    for (int t = start; t < tInc; t++) { // Parallelisable
        // Not all elements in fourChInterp will be used (last column and row)
        for (int i = 0; i < chRows - 1; i++) {
            for (int j = 0; j < chCols - 1; j++) { 
                int ch = i*chRows + j;
                fourChInterp[ch + t*NChannels] = interpolateFourChannels(V + t*NChannels, ch);
            } 
        }
    } 
    return fourChInterp;  
}

inline int InterpDetection::interpolateFiveChannels(int* V, int ch) {
    // Average using w_cs/(3 + w_cs)*the center and 1/(3 + w_cs)* the three largest
    // surrounding channels.
    // · A ·    · · ·
    // D E B -> · F ·
    // · C ·    · · ·
    // Where F = E*w_cs/(3 + w_cs) + (sum(A..D) - max(A..D))/(3 + w_cs)
    
    int values[] = {V[ch - chCols], V[ch + 1], 
                    V[ch + chCols], V[ch - 1]};        
    int interp = (scale*V[ch]*w_cs)/(3 + w_cs);
    int maxValue = INT_MIN;
    bool outlier = false;
    for (auto v : values) { // Add all the values and find the minimum
        if (v == outlierMark) {
            outlier = true;
            break;
        }
        int weightedValue = (scale*v)/(3 + w_cs);
        interp += weightedValue;
        maxValue = max(maxValue, weightedValue);
    }
    if (outlier) {
        return outlierMark;
    }
    else {
        return interp - maxValue;
    }
}

int* InterpDetection::computeFiveChInterp(int* V, int start, int tInc) {
    auto fiveChInterp = new int[NChannels*tInc]; // Has more elements than needed
    for (int t = start; t < tInc; t++) { // Parallelisable
        // Not all elements in fiveChInterp will be used (first/last columns and rows)
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                int ch = i*chRows + j;                                
                fiveChInterp[ch + t*NChannels] = interpolateFiveChannels(V + t*NChannels, ch);
            }
        }        
    }
    return fiveChInterp;  
}

inline bool InterpDetection::isOutlier(unsigned short v) {
    return v >= max_out_threshold || v <= min_out_threshold;
} 

int* InterpDetection::computeMedian(unsigned short* vm, int tInc) {
    auto vGlobal = new int[tInc];
    for (int t = 0; t < tInc; t++) {
        // Take the mean value for each time-step
        int n = 0;
        int sum = 0;
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

inline void InterpDetection::updateBaseline(int v, int ch) {
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

void InterpDetection::initialiseVMovingAvg(unsigned short* vm, int* vGlobal) {
    vMovingAvg = new int[NChannels];
    fill_n(vMovingAvg, NChannels, 0);
    for (int ch = 0; ch < NChannels; ++ch) {
        for (int t = 0; t < movingWindowLength; ++t) {
            auto v = vm[ch + t*NChannels];
            if (isOutlier(v)) {
                // In case of outlier, use the global voltage
                vMovingAvg[ch] += vGlobal[t]/movingWindowLength;
            }
            else {
                vMovingAvg[ch] += (vm[ch + t*NChannels]*scale)/movingWindowLength;
            }
        }
    }
}

void InterpDetection::initialiseVGlobalMovingAvg(int* vGlobal) {
    vGlobalMovingAvgInterIt = 0;
    for (int t = 0; t < movingWindowLength; ++t) {
        vGlobalMovingAvgInterIt += vGlobal[t]/movingWindowLength;
        
    }
}

/**
int* InterpDetection::computeMedianDerivative(int* vGlobal, int tInc) {
    auto vGlobalDerivative = new int[tInc];
    adjacent_difference(vGlobal, vGlobal + tInc, vGlobalDerivative);
    vGlobalDerivative[0] = 0;
    return vGlobalDerivative;
}
*/

int* InterpDetection::preprocessData(unsigned short* vm, int* vGlobal, int start, int tInc) {
    auto Qdiff = new int[NChannels*tInc];
    auto Qmax = new int[NChannels*tInc];
    auto vGlobalMovingAvg = new int[tInc];
    vGlobalMovingAvg[start - 1] = vGlobalMovingAvgInterIt;
    // Compute global moving average
    for (int t = start; t < tInc; t++) {
        vGlobalMovingAvg[t] = vGlobalMovingAvg[t - 1] + 
            (vGlobal[t] - vGlobal[t - movingWindowLength])/movingWindowLength;
    }
 
    for (int ch = 0; ch < NChannels; ch++) { // Independent across channels

        for (int t = start; t < tInc; t++) {
            int QdiffValue;
              
            auto v = vm[ch + t*NChannels];
            bool outlier = isOutlier(v);
            v = v*scale;
            
            if (outlier) {
                // Update moving avg
                vMovingAvg[ch] += (vGlobal[t]-vGlobal[t - movingWindowLength])/movingWindowLength;

                // Mark outlier
                outlierWait[ch] = outlierWaitingTime;
                QdiffValue = outlierMark;
            }
            else {
                // Update moving avg
                auto v_old = vm[ch + (t-movingWindowLength)*NChannels];
                if (isOutlier(v_old)) {
                    vMovingAvg[ch] += ((int) vm[ch + t*NChannels]*scale - 
                                       vGlobal[t - movingWindowLength]) / movingWindowLength;
                }
                else {
                    vMovingAvg[ch] += ((int) vm[ch + t*NChannels]*scale - 
                                       (int) v_old*scale) / movingWindowLength;
                }                

                // Compute baseline
                updateBaseline(v, ch);
                if (outlierWait[ch] > 0) {
                    outlierWait[ch] -= 1;
                    QdiffValue = outlierMark;
                }
                else {
                    QdiffValue = vMovingAvg[ch] - vGlobalMovingAvg[t] - baseline[ch]; // +QmPreD[ch] -Vbias constant
                }
            }
            // Normalised by the variability estimate
            Qdiff[ch + t*NChannels] = QdiffValue / variability[ch]; 
            if (t == start) {
                Qmax[ch + t*NChannels] = Qdiff[ch + t*NChannels];
            }
            else {
                Qmax[ch + t*NChannels] = max(Qdiff[ch + (t-1)*NChannels], Qdiff[ch + t*NChannels]);
            }
        }
    }

    // Prepare putative next iteration
    vGlobalMovingAvgInterIt = vGlobalMovingAvg[tInc - 1];

    // Free memory
    delete[] Qdiff;
    delete[] vGlobalMovingAvg;

    return Qmax;
}

void InterpDetection::findSpikes(int* fourChInterp, int* fiveChInterp, int start, int tInc) {
    auto spikeDelay = new int[NChannels];
    fill_n(spikeDelay, NChannels, -1);

    // Spike delay indicates the number of frames since the spike:
    // -1 -> No recent spike
    //  0 -> Just updated
    //  X -> Number of frames since last spike (waiting for repolarisation)

    for (int t = 0; t < tInc; t++) {
        
        // Parallelisable in an alternate sweep
        
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                // For each channel
                int ch = i*chRows + j;
                auto v = fiveChInterp[ch + t*NChannels];
                // Detection threshold (spike)
                cout << "v:" << v << "\ttheta:"<< theta << "\ttheta_b:" << theta_b << endl;
                if (v < theta) {
                    spikeDelay[ch] = 1;
                }
                // Repolarisation threshold (with previous spike)
                else if (v > theta_b && spikeDelay[ch] != -1) {
                    bool collides = registry.collides(v, i, j);
                    int amplitude = fiveChInterp[ch + (t-spikeDelay[ch])*NChannels];
                    registry.addSpike(amplitude, t, i, j, collides);
                    spikeDelay[ch] = -1;
                }
                // Regular case
                else {
                    if (spikeDelay[ch] > 0) {
                        spikeDelay[ch] += 1;
                        if (spikeDelay[ch] > maxSpikeDelay) {
                            spikeDelay[ch] = -1;
                        } 
                    }
                }           
            }
        }

        // Print output
        registry.pruneOldSpikes(t);
    }

    delete[] spikeDelay;

}

void InterpDetection::detect(unsigned short* vm, int t0, int t1) {    
    // vm is indexed as follows:
    //     vm[channel + tInc*nChannels]
    // or
    //     vm[i*chRows + j + tInc*nChannels]

    // Duration of the interval
    int tInc = t1 - t0;
    int startFrame = 0;

    // Initialisation    
    auto vGlobal = computeMedian(vm, tInc);    
    if (not detectionInitialised) { // First detection
        cout << "Initialising variables..." << endl;
        initialiseVMovingAvg(vm, vGlobal);
        initialiseVGlobalMovingAvg(vGlobal);
        registry.initialise(tau_event, tau_coinc, chRows, chCols, chunkSize);
        startFrame += movingWindowLength;
        detectionInitialised = true;
    }
 
    // Compute baselines (with outlierMarks)
    cout << "Computing baselines..." << endl;
    auto Qmax = preprocessData(vm, vGlobal, startFrame, tInc);

    // Compute interpolations with the baselines (max over two consecutive frames)
    cout << "Computing interpolations..." << endl;
    auto fourChInterp = computeFourChInterp(Qmax, startFrame, tInc);
    auto fiveChInterp = computeFiveChInterp(Qmax, startFrame, tInc);

    cout << "Finding spikes..." << endl;
    findSpikes(fourChInterp, fiveChInterp, startFrame, tInc);    

    // Free memory
    delete[] vGlobal;
    delete[] Qmax;
    delete[] fourChInterp;
    delete[] fiveChInterp;    
}

} // End of namespace SpkDslowFilter