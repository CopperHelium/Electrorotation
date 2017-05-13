# Electrorotation

**
Chi-squared methods 

codes are inclueded in the "Chi_test.m" and "stepfit_test.m"(function)

***
**
Wavelet Transform


When choosing the files, see the upper comments of the jumped window for selecting to choose the right file formats.

**
First, run the "convertAVI.m".

Second, run the "TdmsWavelet.m" and "VideoWavelet.m" to get the wavelet coefficients for tdms and avi files (I think it is better store these parameters for further use.)
 
Third, run the "CombWavelet.m" to combine these two data formats. 

Fourth, run the 'PeakSteps.m' to get the peaks and steps for combined data.

Finally, run the "Com_Plot.m" to plot the combined data. 


***

I only test very few data sets up to now. Errors may happen. Ask me whenever bugs arise.

***

For Chi-squared step finding methods, just run the "Chi_test.m" 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Here are two main .m files, one is for tdms files(Freq_tdms.m) and another one is for the videos(Freq_avi.m).
Including the main codes, I call several functions as followings:

Center.m:
A function for determing the center position for the cells in the video. 

peakfinder.m: 
A function downloaded from mathwork package to find the peak for the magnitude scalogram of the wavelet transform.

ps_plot.m:
A function to plot the peaks, white dashed lines and the steps.

smooth1.m: 
A step modified function. I didn't use it in the updated analysis.

stepmodify.m:
A function for step modifying. The main tasks it does is to find the spikes in the results of peakfinder. 

wavelet_spectum.m (Well, I spelled it wrong...spectrum):  
A function to plot the magnitude scalogram of the wavelet transform.

wltpeakfinder.m:
A modified peakfiner function, which call peakfinder.m in the codes. 

For run:
Run the "installPottslab.m" in the Pottslab0.5 folder first, then just run the "Freq_tdms.m" and "Freq_avi.m" by choosing the data folder simply. The dissociation and recovery plots will be saved automatedly.  
