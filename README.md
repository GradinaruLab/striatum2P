# striatum2P

Code and sample data used to perform analysis and generate figures in Badimon et al., 2020

The repository contains three Matlab scripts used to calculate the three main parameters reported in our paper: event rate, spatial correlation and synchronicity of calcium transients in the striatum.

The folder Required_Scripts contains dependencies required for calculation of these parameters while the folder Example_Cell_Traces contains sample data on which the scripts may be run.

Basic_Extraction_Analysis.m calculates event rates per minute and average magnitude of calcium traces obtained after source extraction from Suite2P. Use any of the sample data files in Required_Scripts to run this script.


Batch_Process_Spatial_Correlation_Analysis.m and Batch_Process_Synchrony_Analysis.m are used to obtain spatial correlation and synchrony measurements reported in our paper for any number of files obtained from Suite2P. Point them to the directory Required_Scripts to run this script.

If there are any questions/problems, please contact Adi: adi.nair@caltech.edu and Xinhong Chen xchen3@caltech.edu


