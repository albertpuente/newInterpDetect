#include "SpkDslowFilter.h"

namespace SpkDslowFilter {

InterpDetection::InterpDetection(int rows, int cols) {
    chRows = rows;
    chCols = cols;
    NChannels = rows*cols;
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
        minValue = std::min(minValue, v);        
    }
    return (interp - minValue)/3;
}

unsigned short* InterpDetection::computeFourChInterp(unsigned short* vm, long tInc) {
    auto fourChInterp = new unsigned short[chRows*chCols*tInc]; 
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
        minValue = std::min(minValue, weightedValue);        
    }
    return interp - minValue;
}

unsigned short* InterpDetection::computeFiveChInterp(unsigned short* vm, long tInc) {
    auto fiveChInterp = new unsigned short[chRows*chCols*tInc]; // Has more elements than needed
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

void InterpDetection::detect(unsigned short* vm, long t0, long t1) {    
    // vm is indexed as follows:
    //     vm[channel + tInc*nChannels]
    // or
    //     vm[i*chRows + j + tInc*nChannels]

    // Duration of the interval
    int tInc = t1 - t0;

    // The first row and column are empty.    
    auto fourChInterp = computeFourChInterp(vm, tInc);

    // The first/last rows and colums are empty.
    auto fiveChInterp = computeFiveChInterp(vm, tInc);

}

}
