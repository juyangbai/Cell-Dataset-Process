%% Objective of the code: Convert PWS image to D-map, Lmax-map and  Nf_map. Refer to:
% https://opg.optica.org/ol/fulltext.cfm?uri=ol-45-17-4810&id=437557
% Written by: Vasundhara (07-08-22)

%% Inputs:
% 1. cubeRMS:Image cube from PWS; 
% 2. BW: Nuclear ROI drawn using PWS analysis GUI
% 3. noise: Estimates using background signal using PWS analysis GUI

%% How to call the function:
% [cubeRms, BW] = load_pws(folder, analysis, roi);
% D_map = RMSMapToDmap(cubeRms,BW,noise);

%% Output:
% Map of D-values Lmax and Nf per pixel.

function D_map = RMSMapToDmap(cubeRms, BW, noise)
%UNTITLED3 Summary of this function goes here
%noise =0.04; % this needs to be estimated for every batch of cell from the background

%% Noise subtraction:
RMSImage=sqrt(abs(cubeRms.^2-noise.^2));

%% RMS to D conversion:
phi=0.35;Nf=1e6;thickIn = 2;Sigma_min=0.02;Sigma_max=0.6; 
%Phi is estimated based on A549 cells (between 0.3-0.4).Thickness was
%measured for HCT116 to be 2 micron by 3D optical profilometer.
ncCenterLambda = 2 * pi / ((2 * pi / 450 + 2 * pi / 700) / 2);
sigma = linspace(0.02, Sigma_max);
sMin = min(sigma);
sMax = max(sigma);
liveCellRI = S2D.RIDefinition.createFromGladstoneDale(1.337, phi);
nuSys = S2D.SystemConfiguration(liveCellRI, 0.52, 1.49, 585, true, true);
system_config=nuSys;
[dOut, dCorrected, Nf_expected, lmax_corrected] = SigmaToD_AllInputs(RMSImage, system_config, Nf, thickIn);
D_map=dCorrected;

end