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
    
//    Define Sample Rate & Sine Wave Frequencies
    int sampleRate = 44100;
    double frequency = 100;
    
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
    std::array<double, 44100> wave;
    
//    Declare Amount of Harmonics of Sine
    int nHarmonics = 20;
    double Nyquist = 44100 / 2;
    std::vector<SineOsc> oscillators;
    for(int i = 0; i < nHarmonics; i++){
//        Call the SineOsc Class(Fs=44.1kHz, Amplitude=0.5, Phase = -π/2)
        oscillators.push_back(SineOsc(44100, 0.5, -M_PI / 2));
    }
    for(int i = 0; i < 44100; i++){
        wave[i] = 0;
        for(int j = 1; j < nHarmonics; j ++){
            int h = 2 * j - 1;
            if((h * frequency) < Nyquist){
//                Apply process with Signal frequency 100
                wave[i] = wave[i] + (1.0 / h)*oscillators[j].process(h * frequency);
            }
        }
    }
    
//    WAV File Creation
    FILE* fptr;
    fptr = fopen("sine_100_0.5_0.3_out.wav", "wb");
    fwrite(&wavHeader, sizeof(wavHeader), 1, fptr);
    for (int i = 0; i < 44100; i++) {
        int value = (int)(32767 * wave[i]);
        int size = 2;
            for (; size; --size, value >>= 8) {
                fwrite(&value, 1, 1, fptr);
            }
    }
    fclose(fptr);
    
//    cout final data parameters
    cout<<"RIFF[4]:        "<<wavHeader.RIFF<<endl;
    cout<<"Chunksize:      "<<wavHeader.ChunkSize<<endl;
    cout<<"WAVE[4]:        "<<wavHeader.WAVE<<endl;
    cout<<"fmt[4]:         "<<wavHeader.fmt<<endl;
    cout<<"Subchunk1Size:  "<<wavHeader.Subchunk1Size<<endl;
    cout<<"AudioFormat:    "<<wavHeader.AudioFormat<<endl;
    cout<<"NumOfChan:      "<<wavHeader.NumOfChan<<endl;
    cout<<"SamplesPerSec:  "<<wavHeader.SamplesPerSec<<endl;
    cout<<"bytesPerSec:    "<<wavHeader.bytesPerSec<<endl;
    cout<<"blockAlign:     "<<wavHeader.blockAlign<<endl;
    cout<<"bitsPerSample:  "<<wavHeader.bitsPerSample<<endl;
    cout<<"Subchunk2ID[4]: "<<wavHeader.Subchunk2ID<<endl;
    cout<<"Subchunk2Size:  "<<wavHeader.Subchunk2Size<<endl;
    
    return 0;
}
