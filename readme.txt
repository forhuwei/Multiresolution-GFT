% *************************************************************
% author:   Wei Hu                                                   **
% date:     12/20/2012                                             **
% modified: 01/01/2015                                         **
% purpose:  codes for multi-resolution GFT          **
% *************************************************************


Introduction:

This directory contains Matlab files that implement the multi-resolution graph Fourier transform coding as described in the following paper:

Wei Hu, Gene Cheung, Antonio Ortega, and Oscar C. Au, “Multi-resolution Graph Fourier Transform for Compression of Piecewise Smooth Images,” IEEE Transactions on Image Processing (TIP), vol. 24, no. 1, pp. 419-433, Jan. 2015.

You are free to use / modify / distribute the code for non-profit academic research purposes. We only ask that you cite the above paper in your publications. 

If you have comments / suggestions about the code, please send emails to me (forhuwei@pku.edu.cn). 


Code Overview:

main files:

1. Training
"collectTableGBT": train the lookup table for GTs

2. Testing
"main_test": two transform modes (8x8 DCT, 4x4 LR-GFT)
