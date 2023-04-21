# Speckle-REMI
Matlab code and sample dataset for fast estimation of multi-exposure speckle imaging data
Test dataset uses the 50us-80ms exposure times from literature. The fullMESIestimates.m and fullMESIestimates_spatial.m scripts allow for application of REMI and sREMI, respectively.
Speckle matrices are to be passed in as an MxNxE matrix with E being the number of exposures with values as K^2 from basic speckle contrast.
T is to be passed in as a 1x1xE matrix.
To call the functions, use the format:
  estimates = fullMESIestimates(K2,T) or estimates = fullMESIestimates_spatial(K2,T)
tc = estimates(:,:,1), beta = estimates(:,:,2), rho = estimates(:,:,3), Dmu = estimates(:,:,4)

dataset can be found here:
https://drive.google.com/file/d/1DhWFWuH9RY2P7_ktaK7oEu2G1wA4Ozvk/view?usp=share_link
