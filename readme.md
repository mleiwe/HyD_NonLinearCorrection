# Non-Linear  Correction Readme
## Download Example Data (if needed)
If you want to try with the demo data provided in the sub-folder here . Alternatively you can create your own image with a single channel and increase the laser power incrementally.
Requirements
## Works with MATLAB v9.11 and above, you may also need…
* Image Processing Toolbox (v11.4)
* Statistics and Machine Learning Toolbox (v12.2) 
* [Bioformats MATLAB toolbox](https://www.openmicroscopy.org/bio-formats/downloads/)
  
## Calculate the equation weights
In the MATLAB Command Window type
  
  `[modelfunc,beta,OverallRMSE]=mnl_CalculateNonLinearEvaluation;`

You will then be asked to select your reference tiff stack.

Several figures will be produced, 
* Figure 1: The intensity changes of the same pixel as the laser power increases (threshold, baseline intensity >500 AU)
* Figure 2: Figure to show the correction of known values and differences in normalisation parameters
* Figure 3: Final comparison of original intensity values to the predicted values if true
* Figure 4: Fit to an exponential, with the residuals provided to view the accuracy of the fit
* Figure 5: Further exploration of the residuals binned according to their original intensity.

There will also be three variables produced…
* modelfunc : This is the equation used to fit the correction
* beta : The covariates for the modelfunc equation
* OversallRMSE – The root mean square error of the fit

*Make sure to save the workspace.* In the next code you will be required to upload this workspace (NB make sure not to change the variable names). If this will be used as part of QDyeFinder saving this workspace is particularly important.
Our example workspace is called ‘BetaAndFunction’

## Correct the Image
Finally, you can run
  
  `mnl_CorrectImageForNonLinearity(Data,dimOrder,ThreshVal,func,beta,fn);`

This will automatically save the images as tiffs. (NB 4D, xycz, images will be saved as separate xyz images).
Variables are described within the code documentation.
