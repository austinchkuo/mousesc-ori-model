# Mouse SC orientation selectivity model simulations

Code associated with the study: Shared computational principles for mouse superior colliculus and primate population orientation selectivity

### Usage:

Run each of the functions in the directory "figsNoGUI" to generate the population plots in the desired figure. \
(e.g., running ```fig3_noGUI.m``` generates the orientation preference plots in Figure 3.)

To run the bar stimulus simulations, run the script "simulateBarStimuli" in the fig6_barSims directory.

You will need to download the pre-generated simulated neural responses from [OSF](https://osf.io/bu9cm/) and add the directory to your MATLAB path.

Select different stimulus and receptive field parameters from the pre-generated set of neural responses by changing the loadfile specified in the comments of each figure script.

All code was validated to work with MATLAB 2023b and 2024a.

Please contact Austin Kuo (achkuo@stanford.edu) with any issues in running simulations.
