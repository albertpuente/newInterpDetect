#include "SpkDslowFilter.h"

namespace SpkDslowFilter {

InterpDetection::InterpDetection(int rows, int cols, double samplingRate) {

    // Detection size

    chRows = rows;
    chCols = cols;
    NChannels = rows*cols;

    // Algorithm parameters

    scale = 64;                 // Increase resolution in detection
    theta =  -14*scale;         // Detection threshold (PROVISIONAL)
    theta_b = 0;                // Repolarisation threshold
    w_cs = 4;                   // Center/surround weighting (5-channel interpolation)
	// theta_ev;                // Minimum depolarisation area	
    tau_pre = 7;                // (1 ms) Cut-out window before peak
    tau_post = 15;              // (2.2 ms) Cut-out window after peak
    tau_coinc = 3;              // (0.27 ms) Window for coincident events
    tau_event = 8;              // (1 ms) Characteristic event length	
    f_s = samplingRate;
    maxSpikeDelay = f_s/1000;
    f_v = 4;                    //  Variability increase factor <---------------------- PROVISIONAL
    f_b = 32;                   //  Baseline increase factor <------------------------- PROVISIONAL
    movingWindowLength = 5;     // In frames    
    minVariability = 2*scale; // <----------------------------------------------------- PROVISIONAL
    initialVariability = 3.12* scale; // <--------------------------------------------- PROVISIONAL
    max_out_threshold = 4000;
    min_out_threshold = 10;

    // Execution variables

    variability = new int[NChannels];
    fill_n(variability, NChannels, initialVariability);
    baseline = new int[NChannels];    
    detectionInitialised = false;
    outlierWait = new int[NChannels];
    chunkSize = 8; // Default: 64 chunks (8x8 chunks of 8x8 channels each)
    nSpikes = 0;
    nthreads = 8; // <------------------------------------------------------------------- PROVISIONAL
    threads = new std::thread[nthreads];
    startDetectionFrame = 50; // Provisional, equal to tCut
}

inline int InterpDetection::interpolateFourChannels(int* V, int t, int ch) {
    // Average the three largest absolute signal amplitudes
    // A B -> E ·
    // C D    · ·
    // Where E = sum(A,B,C,D) - max(A,B,C,D)
    int values[] = {V[ch + t*NChannels],          V[ch + 1 + t*NChannels],
                    V[ch + chCols + t*NChannels], V[ch + chCols + 1 + t*NChannels]};        
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
        return (interp - maxValue)/3;
    }
}

void InterpDetection::computeFourChInterpThread(int threadID, int* fourChInterp, int* V, 
                                                int start, int tInc) {

    int chunkSize = std::ceil( (float) (tInc - start)/ (float) nthreads);   

    // Loop accross all channels associated to this thread
    for (int t = threadID*chunkSize; t < tInc and t < (threadID+1)*chunkSize; t++) { 

        // Not all elements in fourChInterp will be used (last column and row)
        for (int i = 0; i < chRows - 1; i++) {
            for (int j = 0; j < chCols - 1; j++) { 
                int ch = i*chRows + j;
                fourChInterp[ch + t*NChannels] = interpolateFourChannels(V, t, ch);
            } 
        }
    }
}

int* InterpDetection::computeFourChInterp(int* V, int start, int tInc) {
    auto fourChInterp = new int[NChannels*tInc]; 
    /*
    // Call parallel preprocess and computation of Qmin and Qdiff
    for (int threadID = 0; threadID < nthreads; threadID++) {
        threads[threadID] = std::thread( [=] { 
            computeFourChInterpThread(threadID, fourChInterp, V, start, tInc); 
        });
    }

    // Wait for all threads
    for (int threadID = 0; threadID < nthreads; threadID++) { 
        threads[threadID].join();
    } 
    */
    return fourChInterp;  
}

int InterpDetection::interpolateFiveChannels(int* V, int t, int ch) {
    // Average using w_cs/(3 + w_cs)*the center and 1/(3 + w_cs)* the three largest
    // surrounding channels.
    // · A ·    · · ·
    // D E B -> · F ·
    // · C ·    · · ·
    // Where F = E*w_cs/(3 + w_cs) + (sum(A..D) - max(A..D))/(3 + w_cs)
    
    int values[] = {V[ch - chCols + t*NChannels], V[ch + 1 + t*NChannels], 
                    V[ch + chCols + t*NChannels], V[ch - 1 + t*NChannels]};        
    int interp = (V[ch + t*NChannels]*w_cs)/(3 + w_cs);
    int maxValue = INT_MIN;
    bool outlier = false;
    for (auto v : values) { // Add all the values and find the minimum
        if (v == outlierMark) {
            outlier = true;
            break;
        }
        int weightedValue = v/(3 + w_cs);
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

int* InterpDetection::computeFiveChInterpSYCL(int* V, int start, int tInc)  {
    auto fiveChInterp = new int[NChannels*tInc]; // Has more elements than needed    
   
    try {  // SYCL Block
        using namespace cl::sycl;

        // Choose the best available device 
        default_selector selector;
        
        // Queue to enqueue command groups (std::queue type conflict?)
        cl::sycl::queue Q(selector);
        
        // Data buffers

        // Scalars
        buffer<int, 1> NChannelsBuffer (&NChannels, range<1> (1));
        buffer<int, 1> chRowsBuffer (&chRows, range<1> (1));
        buffer<int, 1> chColsBuffer (&chCols, range<1> (1));
        buffer<int, 1> startBuffer (&start, range<1> (1));
        buffer<int, 1> tIncBuffer (&tInc, range<1> (1));
        buffer<int, 1> w_csBuffer (&start, range<1> (1));
        buffer<int, 1> outlierMarkBuffer (&start, range<1> (1));

        // Arrays
        buffer<int, 1> VBuffer(V, range<1> (NChannels*tInc));
        buffer<int, 1> fiveChInterpBuffer(fiveChInterp, range<1> (NChannels*tInc));

        // Command group 
        Q.submit([&] (cl::sycl::handler& cgh) {
            
            // Data access to the buffers from the device with different permissions

            // Scalars
            auto NChannelsPtr = NChannelsBuffer.get_access< access::mode::read >(cgh);  
            auto chRowsPtr = chRowsBuffer.get_access< access::mode::read >(cgh); 
            auto chColsPtr = chColsBuffer.get_access< access::mode::read >(cgh);
            auto startPtr = startBuffer.get_access< access::mode::read >(cgh);
            auto tIncPtr = tIncBuffer.get_access< access::mode::read >(cgh);
            auto w_csPtr = w_csBuffer.get_access< access::mode::read >(cgh);
            auto outlierMarkPtr = outlierMarkBuffer.get_access< access::mode::read >(cgh);
            
            // Arrays
            auto VPtr = VBuffer.get_access< access::mode::read >(cgh);   
            auto fiveChInterpPtr = fiveChInterpBuffer.get_access< access::mode::write >(cgh);

            // Work space definition: 1-dimensional: global range, work group range
            int workGroupRange = 256;
            int globalRange = tInc - start;
            globalRange += workGroupRange - globalRange%workGroupRange; // Make it multiple of the the workGroupRange
            
            auto workSpaceRange = nd_range<1>(range<1>(globalRange), range<1>(workGroupRange)); 
            
            auto interpKernel = ([=](nd_item<1> time) {
                // Kernel                
                                
                // Scalar vars   
                auto NChannels = NChannelsPtr[0];
                auto chRows = chRowsPtr[0];
                auto chCols = chColsPtr[0];
                auto start = startPtr[0];
                auto tInc = tIncPtr[0];
                auto w_cs = w_csPtr[0];
                auto outlierMark = outlierMarkPtr[0];

                auto t = time.get_global()[0] + start;
                
                if (t < tInc) { // Only for kernels in range                    
                
                    for (int i = 1; i < chRows - 1; i++) {
                        for (int j = 1; j < chCols - 1; j++) {
                            int ch = i*chRows + j;

                            int values[] = {VPtr[ch - chCols + t*NChannels], 
                                            VPtr[ch + 1 + t*NChannels], 
                                            VPtr[ch + chCols + t*NChannels], 
                                            VPtr[ch - 1 + t*NChannels]};

                            int interp = (VPtr[ch + t*NChannels]*w_cs)/(3 + w_cs);
                            int maxValue = INT_MIN;
                            bool outlier = false;
                            for (auto v : values) { // Add all the values and find the minimum
                                outlier = outlier || (v == outlierMark);     
                                int weightedValue = v/(3 + w_cs);
                                interp += weightedValue;
                                maxValue = cl::sycl::max(maxValue, weightedValue);
                            }
                            fiveChInterpPtr[ch + t*NChannels] = (!outlier)*(interp-maxValue) + outlier*outlierMark;
                        }
                    }
                }
                
                // End of Kernel
            });
            
            // Call
            cgh.parallel_for<class detect>(workSpaceRange, interpKernel);
        });

    } catch (cl::sycl::exception e) {
        std::cout << "SYCL exception caught: " << e.what();
    }      // End of SYCL Block    
    return fiveChInterp;
}

void InterpDetection::computeFiveChInterpThread(int threadID, int* fiveChInterp, int* V, 
                                                int start, int tInc) {

    int chunkSize = std::ceil( (float) (tInc - start)/ (float) nthreads);   

    // Loop accross all channels associated to this thread
    for (int t = threadID*chunkSize; t < tInc and t < (threadID+1)*chunkSize; t++) { 
        // Not all elements in fiveChInterp will be used (first/last columns and rows)
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                int ch = i*chRows + j;                                
                fiveChInterp[ch + t*NChannels] = interpolateFiveChannels(V, t, ch);
            }
        }
    }
}

int* InterpDetection::computeFiveChInterp(int* V, int start, int tInc) {
    auto fiveChInterp = new int[NChannels*tInc]; // Has more elements than needed
    
    // Call parallel preprocess and computation of Qmin and Qdiff
    for (int threadID = 0; threadID < nthreads; threadID++) {
        threads[threadID] = std::thread( [=] { 
            computeFiveChInterpThread(threadID, fiveChInterp, V, start, tInc); 
        });
    }

    // Wait for all threads
    for (int threadID = 0; threadID < nthreads; threadID++) { 
        threads[threadID].join();
    } 

    return fiveChInterp;  
}

inline bool InterpDetection::isOutlier(unsigned short v) {
    return v >= max_out_threshold or v <= min_out_threshold;
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
    if ((v > baseline[ch] - variability[ch] and v <= baseline[ch]) or
        (v > INT_MIN and v < baseline[ch] - 6*variability[ch]) ) {

        variability[ch] -= f_v;
        variability[ch] = max(variability[ch], minVariability);
    } // v in (baseline - 5variability)
    else if (v > baseline[ch] - 5 *variability[ch] and 
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

        // Initialise baseline with plausible values
        if (vMovingAvg[ch] != 0) {
            baseline[ch] = vMovingAvg[ch];
        }
        else {
            baseline[ch] = vGlobal[movingWindowLength];
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

void InterpDetection::preprocessDataThread(int threadID, unsigned short* vm, int* vGlobal, 
                        int start, int tInc, int* Qdiff, int* Qmin, int* vGlobalMovingAvg) {

    for (int t = start; t < tInc; t++) {
        // Number of channels associated to a thread
        int chunkSize = std::ceil( (float) NChannels/ (float) nthreads);        
        // Loop accross all channels associated to this thread
        for (int ch = threadID*chunkSize; ch < NChannels and ch < (threadID+1)*chunkSize; ch++) { 
             
            int v = (int) vm[ch + t*NChannels];
            bool outlier = isOutlier(v);
            v = v*scale;
            
            if (outlier) {
                // Update moving avg
                vMovingAvg[ch] += (vGlobal[t]-vGlobal[t - movingWindowLength])/movingWindowLength;

                // Mark outlier
                outlierWait[ch] = outlierWaitingTime;
                Qdiff[ch + t*NChannels] = outlierMark;
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
                    Qdiff[ch + t*NChannels] = outlierMark;
                }
                else {
                    // PROVISIONAL Without +QmPreD[ch] -Vbias constant - vGlobalMovingAvg[t] 
                    Qdiff[ch + t*NChannels] = (vMovingAvg[ch] - baseline[ch]); // Normalise? /variability[ch]; 
                    
                }
                
            }
            // Normalised by the variability estimate
            

            

            // Select maximum of consecutive frames
            if (t == start) {
                Qmin[ch + t*NChannels] = Qdiff[ch + t*NChannels];
            }
            else {
                Qmin[ch + t*NChannels] = min(Qdiff[ch + (t-1)*NChannels], Qdiff[ch + t*NChannels]);
            }
        }
    }
}

int* InterpDetection::preprocessData(unsigned short* vm, int* vGlobal, int start, int tInc) {
    auto Qdiff = new int[NChannels*tInc];
    auto Qmin = new int[NChannels*tInc];
    auto vGlobalMovingAvg = new int[tInc];
    vGlobalMovingAvg[start - 1] = vGlobalMovingAvgInterIt;
    // Compute global moving average
    for (int t = start; t < tInc; t++) {
        vGlobalMovingAvg[t] = vGlobalMovingAvg[t - 1] + 
            (vGlobal[t] - vGlobal[t - movingWindowLength])/movingWindowLength;
    }

    // Call parallel preprocess and computation of Qmin and Qdiff
    for (int threadID = 0; threadID < nthreads; threadID++) {
        threads[threadID] = std::thread( [=] { 
            preprocessDataThread(threadID, vm, vGlobal, start, 
                                 tInc, Qdiff, Qmin, vGlobalMovingAvg); 
        });
    }

    // Wait for all threads
    for (int threadID = 0; threadID < nthreads; threadID++) { 
        threads[threadID].join();
    }    

    // Prepare putative next iteration
    vGlobalMovingAvgInterIt = vGlobalMovingAvg[tInc - 1];

    // Free memory
    delete[] Qdiff;
    delete[] vGlobalMovingAvg;

    return Qmin;
}

void InterpDetection::writeOutput(const vector<Spike>& spikes, int* fiveChInterp, 
                                  unsigned short* vm, int t0) {

    for (auto s : spikes) {
        int start = s.t - tau_pre;
        int end = s.t + tau_post;
        int ch = s.chX*chRows + s.chY;

        output << s.t + t0 << " " << ch;        
        // Write 5-interp center
        for (int t = start; t <= end; ++t) {
            output << " " << fiveChInterp[ch + t*NChannels];
        } 
        output << endl;     

        // Write 9 surrounding real channels
        for (int i = s.chX - 1; i < s.chX + 2; i++) {
			for (int j = s.chY - 1; j < s.chY + 2; j++) {	
                ch = i*chRows + j;
                output << "VM";
                for (int t = start; t <= end; ++t) {
                    output << " " << vm[ch + t*NChannels];
                } 
                output << endl;
            }
        }	 
    }
}

void InterpDetection::findSpikes(unsigned short* vm, int* fourChInterp, int* fiveChInterp, 
                                 int start, int t0, int tInc) {

    auto spikeDelay = new int[NChannels];
    fill_n(spikeDelay, NChannels, -1);
    
    auto currentMin = new int[NChannels];

    // spikeDelay indicates the number of frames since the spike:
    // -1 -> No recent spike
    //  0 -> Just updated
    //  X -> Number of frames since last spike (waiting for repolarisation)

    for (int t = startDetectionFrame; t < tInc - tau_post; t++) {
        /* TO-DO
           Compatible with a parallelisation in an 4-step alternate sweep
           across chunks of channels. For example, if 4096 channels and 8x8
           chunks, it can be done in 4 steps with 16 simultaneous threads each.

           X · X ·      · X · X     · · · ·      · · · ·
           · · · ·      · · · ·     X · X ·      · X · X
           X · X ·      · X · X     · · · ·      · · · ·
           · · · ·      · · · ·     X · X ·      · X · X

           The size of the chunk determines a trade-off between number of chunks
           per step vs number of putative spike collisions to check. The size of 
           the chunks must be a least the tau_event characteristic event length so 
           that it is only necessary to check the surrounding chunks for collisions.

           TO-CHECK: Overhead?
        */
                
        for (int i = 1; i < chRows - 1; i++) {
            for (int j = 1; j < chCols - 1; j++) {
                // For each channel
                int ch = i*chRows + j;
                auto v = fiveChInterp[ch + t*NChannels];

                // Outlier filter
                if (v < - 10e7 or v > 10e7) {
                    spikeDelay[ch] = -1;
                }
                // Detection threshold (spike)
                else if (v < theta) { // Better spike in memory
                    if (spikeDelay[ch] > 0 and not (v < currentMin[ch])) {
                        spikeDelay[ch] += 1;
                    }                    
                    else { // New best spike
                        spikeDelay[ch] = 1;
                        currentMin[ch] = v;
                    }
                }
                // Repolarisation threshold (with previous spike)
                else if (v > theta_b and spikeDelay[ch] != -1) {
                    bool relevant = not registry.collides(v, i, j);
                    int peakTime = t-spikeDelay[ch];
                    int amplitude = fiveChInterp[ch + peakTime*NChannels];
                    
                    // Sanity check
                    if (v < - 10e7 or v > 10e7) {
                        cout << "ERROR!: Outlier is going to be used as a valid value!" << endl;
                        exit(1);
                    }
                    //

                    registry.addSpike(amplitude, peakTime, i, j, relevant);
                    spikeDelay[ch] = -1;
                }
                // Regular case
                else {
                    if (spikeDelay[ch] > 0) {
                        if (++spikeDelay[ch] > maxSpikeDelay) {
                            spikeDelay[ch] = -1;
                        } 
                    }
                }           
            }
        }

        // Print output (TO-DO: copy to thread to allow computations while printing)
        auto spikesToOutput = registry.pruneOldSpikes(t, fourChInterp, fiveChInterp);
        nSpikes += spikesToOutput.size();

        writeOutput(spikesToOutput, fiveChInterp, vm, t0);        
    }

    delete[] spikeDelay;
    delete[] currentMin;
}

void InterpDetection::detect(unsigned short* vm, int t0, int t1, int tCut) {    
    // vm is indexed as follows:
    //     vm[channel + tInc*nChannels]
    // or
    //     vm[i*chRows + j + tInc*nChannels]

    // Duration of the interval
    int tInc = t1 - t0;
    int startFrame;

    // Initialisation    
    cout << "Computing median..." << endl;
    auto vGlobal = computeMedian(vm, tInc);    
    
    if (not detectionInitialised) { // First detection
        cout << "Initialising variables..." << endl;
        output.open("spikes.txt");
        initialiseVMovingAvg(vm, vGlobal);
        initialiseVGlobalMovingAvg(vGlobal);
        registry.initialise(tau_event, tau_coinc, chRows, chCols, chunkSize);
        startFrame = movingWindowLength;
        detectionInitialised = true;
    }
    else {
        startFrame = tCut;
    }
 
    // Compute baselines (with outlierMarks)
    cout << "Computing baselines..." << endl;
    auto Qmin = preprocessData(vm, vGlobal, startFrame, tInc);

    // Compute interpolations with the baselines (max over two consecutive frames)
    cout << "Computing interpolations..." << endl;
    auto fourChInterp = computeFourChInterp(Qmin, startFrame, tInc);
    auto fiveChInterp = computeFiveChInterp(Qmin, startFrame, tInc);

    cout << "Finding spikes..." << endl;
    findSpikes(vm, fourChInterp, fiveChInterp, startFrame, t0, tInc);    

    // cout << nSpikes << " spikes." << endl;

    // Free memory between iterations
    delete[] vGlobal;
    delete[] Qmin;
    delete[] fourChInterp;
    delete[] fiveChInterp;

    registry.purge();
}

} // End of namespace SpkDslowFilter