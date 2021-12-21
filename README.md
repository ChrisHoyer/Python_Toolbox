 Python Toolbox

This is a library with different sub-libraries to provide basic functions, which can be used in data and signal processing and to compare measurement results against simulation.

# Content

## Basic_Toolbox
Major focus on plotting and import

 - CSV2Dict: script for importing CSV Files
 - CSV2Area: script for calculating polygon area based on coordinates
 - FitFct_Exp: script that fits data against an Exp Function
 - FindPoint_FitFct: fits higher order Polynom against Dataset to find point
 - Linearization_Point: linearization around one point
 - Linear_Plot: linear plot function with automated labeling
 - Polar_Plot: Polar plot function with automated labeling
 - Histogram_Plot: Plot Histogram
 - Box_Plot: Boxplots with automated labeling
 - SemiLogX_Plot: semilog x plot function with automated labeling
 - Vline_Plot: generates vertical line in plot with label
 - Hline_Plot: generates horizontal line in plot with label
 - Rectangle_Plot: generates rectangle inside plot
 - Align_YAxis: Align two YGrids in one plot
 - FindPoint_NextValue: find nearest point
 - Digitalize_Data: Generates binary stream from data and clock from
 - Extract_DigitalClock: Extracts clock from data and clock from waveform
 - CMPLX2Format: converts complex numbers in different formats
 - Average: Average a Dataset
 - FrequencyFiltering: Frequency Filtering of Data
 - String2List: Generates a List out of a csv with one line
 - XYZ_Plot: generates 3D Plot (X,Y,Z) for i.e. waterfall diagrams
 - Spectrum_Minimizer: Generates Mean and Max Spectrum from Dataset
 - MovingFilter: Moving Average, Max, Min Filter
 - printProgressBar: Print progress Bar in Console
 - EngNot: Engineering Notation (Si-Prefix)

## ADS_Toolbox
Functionality to import Keysight ADS Simulation Data (ASCII-File) and calculate different figure of merrits:

 - ImportData: Imports an ASCII file from ads
 - Calc_HarmonicBalance_Single: Calculates IP3, Psat, 1dB Compression Point based on HarmonicBalance Simulation
 - Calc_HarmonicBalance: Uses Calc_HarmonicBalance_Single for multiple frequencies
 - Calc_StabGain: Calculates MSG, MAG, K-Fact, Masons Gain U and Gmax based on 2-Port s2p
 - Calculate_StabCircle: Generates stability circles based on 2-Port s2p
 - MixedModeSparam: Calculates standart and mixed moded matricies based on 4-Port s4p
 - ImportS2P: Data conversion (change in future versions)

## ControlTheory_Toolbox
Major fucus on control theorie and plotting

 - Extract_Sympy_1Var: Substitutes Sympy and generates numeric solution
 - BodePlot_FBCTRL: Generate BodePlot out of symbolic feedback transfer function
 - BodePlot: Generate BodePlot out of symbolic transfer function
 - ZeroPole_Plot: Generate ZeroPole Plot
 - StepResponse: Generate Step Response with Heaviside Fct from symbolic transfer function
 - Substitute_Datatype: Substitute constant values with symboles
 - ReSubstitute_Datatype: Resubstitute constant values with symboles