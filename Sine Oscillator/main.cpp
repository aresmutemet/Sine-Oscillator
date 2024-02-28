//
//  main.cpp
//  Sine Oscillator
//
//  Created by Areş Mutemet on 28.02.2023.
//

#include <iostream>
#include <vector>
#include <array>
#include <cstdint>
#include <math.h>
#include "stdio.h"

using namespace std;

//                          Sine Wave Oscillator Class

//                          x[n] = Acos(2πfnTs + ϕ) -----> Discrete Time Signal Representation

class SineOsc {
    public :
    SineOsc ( double sampleRate , double amplitude , double phase) {
        
//        Set the time step, amplitude and phase of oscillator
        mTs = 1 / sampleRate;
        mAmp = amplitude;
        mPhase = phase;
    }
    double process (double freq) {
        double newPhase = mPhase + freq * 2 * M_PI * mTs;
        mPhase = fmod(newPhase , 2 * M_PI);
        return mAmp * cos(mPhase);
    }
private:
    double mTs;
    double mAmp;
    double mPhase;
};

/*
 Butterworth filter Transfer Function
 :
                B0 + B1*z^(-1) + B2*z^(-2)
 H(z) = ---------------------------------------
           1 + A1*z^(-1) + A2*z^(-2)

 */

class ButterworthFilter {
public:
    
//    Initialize the Butterworth filter coefficients based on the desired cutoff sampling frequencies
    ButterworthFilter(double cutoff_freq, double sampling_freq) {
        double cutoff_norm = cutoff_freq / sampling_freq;
        double c = tan(M_PI * cutoff_norm);
        
//        Calculate the feedforward coefficients (B0, B1, B2) and feedback coefficients (A0, A1)
        mB0 = 1.0 / (1.0 + sqrt(2.0) * c + pow(c, 2.0));
        mB1 = -2.0 * mB0;
        mB2 = mB0;
        mA0 = 2.0 * mB0 * (pow(c, 2.0) - 1.0);
        mA1 = mB0 * (1.0 - sqrt(2.0) * c + pow(c, 2.0));
        mXnz1 = mXnz2 = mYnz1 = mYnz2 = 0.0;
    }

    void process_samples(double *wave, int length) {
        for (int i = 0; i < length; i++) {
            double x = wave[i];
            
//            Calculate the output sample using the filter difference equation and update values for the next iteration
            double y = mB0 * x + mB1 * mXnz1 + mB2 * mXnz2 - mA0 * mYnz1 - mA1 * mYnz2;
            mXnz2 = mXnz1;
            mXnz1 = x;
            mYnz2 = mYnz1;
            mYnz1 = y;
            wave[i] = y;
        }
    }
    

private:
    double mB0, mB1, mB2, mA0, mA1;
    double mXnz1, mXnz2, mYnz1, mYnz2;
};

//Wav File Data Definition
struct WAV_HEADER {
    uint8_t RIFF[4];            // RIFF Header Magic header
    uint32_t ChunkSize;         // RIFF Chunk Size
    uint8_t WAVE[4];            // WAVE Header
    uint8_t fmt[4];             // FMT header
    uint32_t Subchunk1Size;     // Size of fmt chunk
    uint16_t AudioFormat;       // Audio format 1=PCM,6=mulaw,7=alaw
    uint16_t NumOfChan;         // Number of Channels 1=mono, 2 =stereo
    uint32_t SamplesPerSec;     // Sampling Frequency in Hz
    uint32_t bytesPerSec;       // bytes per second
    uint16_t blockAlign;        // 2=16-bit mono, 4=16-bit stereo
    uint16_t bitsPerSample;     // Number of bits per sample
    uint8_t Subchunk2ID[4];     // "data" string
    uint32_t Subchunk2Size;     // Sampled data length
} wav_hdr;

int main() {
    
    WAV_HEADER wavHeader;
    
//    Define Sample Rate
    int sampleRate = 44100;
    
//    WAV Data Parameters Updated(Mono, 44.1 kHz)
    std::strncpy((char*)wavHeader.RIFF,"RIFF",4);
    wavHeader.ChunkSize=(sampleRate*2)+36;
    std::strncpy((char*)wavHeader.WAVE,"WAVE",4);
    std::strncpy((char*)wavHeader.fmt,"fmt ",4);
    wavHeader.Subchunk1Size=16;
    wavHeader.AudioFormat=1;
    wavHeader.NumOfChan=1;
    wavHeader.SamplesPerSec=sampleRate;
    wavHeader.bytesPerSec=2*sampleRate;
    wavHeader.blockAlign=2;
    wavHeader.bitsPerSample=16;
    wavHeader.Subchunk2Size=sampleRate*2;
    std::strncpy((char*)wavHeader.Subchunk2ID,"data",4);
    
//    Declare Signal Frequency and amount of harmonics of sine and Nyquist
    double frequency = 300;
    int nHarmonics = 10;
    double Nyquist = sampleRate / 2;
    
//    Declare Amplitude and Phase of the Signal
    double amplitude = 0.4;
    double phase = 0.3;
    std::vector<SineOsc> oscillators;
    for(int i = 0; i < nHarmonics; i++){
//        Call the SineOsc Class
        oscillators.push_back(SineOsc(sampleRate, amplitude, phase));
    }
    
    std::array<double, 44100> wave = { 0.0 };
    int N = wave.size();
    for(int i = 0; i < N; i++){
        wave[i] = 0.0;
        for(int j = 1; j < nHarmonics; j ++){
            int h = 2 * j - 1;
            if((h * frequency) < Nyquist){
//                Apply process with Signal frequency 100
                wave[i] = wave[i] + (1.0 / h)*oscillators[j].process(h * frequency);
            }
        }
    }
    
    //    Unfiltered WAV File Creation
        FILE* fptr;
        fptr = fopen("Orginal_Sine_Wave.wav", "wb");
        fwrite(&wavHeader, sizeof(wavHeader), 1, fptr);
        for (int i = 0; i < 44100; i++) {
            int value = (int)(32767 * wave[i]);
            int size = 2;
                for (; size; --size, value >>= 8) {
                    fwrite(&value, 1, 1, fptr);
                }
        }
        fclose(fptr);
    
//    Apply Filter with 3kHz Cut-off Frequency
    ButterworthFilter filter(3000.0, 44100.0);
    int length_bytes = sizeof(wave);
    int length_samples = length_bytes / sizeof(double);
    filter.process_samples(wave.data(), length_samples);
    
//    Filtered WAV File Creation
    FILE* fptr_filtered;
    fptr_filtered = fopen("Filtered_Sine_Wave.wav", "wb");
    fwrite(&wavHeader, sizeof(wavHeader), 1, fptr_filtered);
    for (int i = 0; i < 44100; i++) {
        int value = (int)(32767 * wave[i]);
        int size = 2;
            for (; size; --size, value >>= 8) {
                fwrite(&value, 1, 1, fptr_filtered);
            }
    }
    fclose(fptr_filtered);
    
    
    return 0;
}
