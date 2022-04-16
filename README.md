# mhd-simulation
Simulation code for electron paramagnetic resonance (EPR) spectroscopy using multiple harmonic detection (MHD)

We developed the Fortran code to simulate a multiple harmonic detection (MHD) receiver system used for continuous-wave (CW) electron paramagnetic resonance (EPR) spectroscopy.

This code was used for the article published in the Journal of X (to be announced).
DOI: to be announced.

We developed and tested our codes in the following environment:

Apple Mac mini (Intel Core i5, 3GHz, memory 8 GB)
Fortran 95 Compiler: Absoft Pro Fortran 2020
Numerical Library: IMSL 2018 Fortran Numerical Libraries

Three programs and sample data are deposited here. The purposes and functions of each code are given below.

(1) noise-snr1_3.f95
This code computes the signal-to-noise ratio (SNR) of the first-derivative EPR absorption spectrum under a given modulation ratio (input). In this computation, the additional noise amplitude (gamma in the article) is set to the peak height of the given EPR absorption profile.

By considering the SNR obtained by this code, the additional noise level can be roughly determined. The SNR for the first-derivative spectrum (SNR1) is inversely proportional to the additional noise level (gamma) approximately.

The output data file ‘gamma-snr2.csv’ is written in the directory ‘./gradput’.
The data file gamma-snr2.csv contains two parameters:
Given modulation ratio, SNR1

(2) noisecheck3.f95
This program simulates multiple harmonic signals with and without additional noise. The SNR of the first-derivative spectrum should be given when the program runs.
This program needs the data file ‘./gradput/gamma-snr2.csv’.

Output:
./gradput/harmonics.csv (harmonic signals with additional noise)
./gradput/nharmonics.csv (noise-free harmonic signals)

(3) snrlw3_rline.95
This code computes the SNR and the peak-to-peak linewidth of the reconstructed spectrum.
When this program runs, the SNR of the first-derivative spectrum given to the program noisecheck3 should be used again.

Output:
./gradput/snrlw2.csv (SNR and the peak-to-peak linewidth as a function of filter passband)
./gradput/reconstructed.csv (the reconstructed spectrum using MHD)

Additional files

The measured spectral data files are given for the program test. These files are loaded in the above programs and should be placed in the same directory as the executable files.

a1nakaoka.csv
a2nakaoka.csv
a3nakaoka.csv
3a1nakaoka.csv
3a3nakaoka.csv


These codes were written by Yamato Mori, Magnetic Resonance Engineering Laboratory, Hokkaido University, Sapporo, Japan.
