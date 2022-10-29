%% Set File Path:
PWS_Original_FilesPath='\\backmanlabnas.myqnapcloud.com\Public\Ali_Vasundhara_Confocal_Data\NEW_TEST_CODE_IMAGES\ConvertToD\PWS_Fluor_Data_All\Cell1';

%% Load PWS Image and Nuclear handrawn-ROI:
[cubeRms_OriginalPWS,BW_OriginalPWS] = load_pws(PWS_Original_FilesPath, 'p0', 'nuc');
figure;imagesc(cubeRms_OriginalPWS);

%% RMS to D conversion:
noise_OriginalPWS=0.02; %estimated from the backgroud of the image, i.e. draw ROI outside cell and find RMS
D_map = RMSMapToDmap(cubeRms_OriginalPWS,BW_OriginalPWS,noise_OriginalPWS);

%% Multiply with nucleus mask
NucMask=sum(BW_OriginalPWS,3);
D_map(NucMask==0)=0;

%% Scaled Color Image of D-map:
figure;imagesc(D_map);