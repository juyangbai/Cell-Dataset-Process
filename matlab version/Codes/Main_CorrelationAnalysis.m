%% Set Filepaths:
clc;clear all;close all;
PWS_Original_FilesPath='C:\Users\Juyang\Project_IVPL\Aggelos_Codes_PWSD_ConfocalDV\Run\PWS_Data';
DV_ConfocalZStack_FilesPath='C:\Users\Juyang\Project_IVPL\Aggelos_Codes_PWSD_ConfocalDV\Run\DV_Data';
Analysis_Cell_FilesPath='C:\Users\Juyang\Project_IVPL\Aggelos_Codes_PWSD_ConfocalDV\Run\Results';

% Set Parameters:
% estimated from the backgroud of the PWS image, i.e. draw ROI(Return on Investment) outside cell and find RMS(Root Mean Square)
noise_PWS=0.02;                         

% Loop through all Cells:
cd(PWS_Original_FilesPath)
Num_Nucleus_Analyzed = dir('Cell*');

for fileindexing =1:length(Num_Nucleus_Analyzed)

    %% PWS Image
    % 1. Load PWS:
    cd(PWS_Original_FilesPath)
    PWS_New_FilesPath = append(PWS_Original_FilesPath, ['\Cell' num2str(fileindexing)]);
    [cubeRms_PWS, BW_Nuc_PWS] = load_pws(PWS_New_FilesPath, 'p0', 'nuc');
    
    % 2. Create RMS-map:
    % locate the single cell
    RMS_map_PWS = cubeRms_PWS;                                      
    figure;imagesc(RMS_map_PWS);
    RMS_map_PWS = RMS_map_PWS.*single(BW_Nuc_PWS);                  % Eliminate non nuclear region
    figure;imagesc(RMS_map_PWS);

    filename = strcat('PWS','.png');
    t = Tiff(filename,'w');
    tagstruct.ImageLength = size(RMS_map_PWS,1);
    tagstruct.ImageWidth = size(RMS_map_PWS,2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct)
    t.write(RMS_map_PWS);
    t.close;

    
    
    %% Mapping PWS and Deltavision stacks
    cd(DV_ConfocalZStack_FilesPath) % Load DV Image Slices and Nuc Mask Slices
    DV_ConfocalZStack_New_FilesPath = append(DV_ConfocalZStack_FilesPath,['\Cell' num2str(fileindexing)]);
    cd(DV_ConfocalZStack_New_FilesPath)
    Image_name = [dir('*_D3D-processed.tif')];
    Nuc_name = [dir('*_NUClabels.mat')];
    Stack_DV = (tiffreadVolume(Image_name.name));
    load(Nuc_name.name);
    BW_Nuc_DV = labels;
    
    % 3. Transform RMS-map after doing correlation analysis for the mask of DV and PWS:
    [Output_cubeRMS_PWS_Transformed, OutputNucMask_PWS_Transformed, Output_DVStackImage_Transformed, OutputNucMask_DV_Transformed, tform, index_slice_maxcorrelation_begin, index_slice_maxcorrelation_end] = MatchPWS_Confocal(cubeRms_PWS, BW_Nuc_PWS, Stack_DV, BW_Nuc_DV, Analysis_Cell_FilesPath, fileindexing);
%     figure; imshow(Output_cubeRMS_PWS_Transformed)
%     x_width=9 ;y_width=7; set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%     Fig6=strcat('PWS','.png');
%     print(gcf,Fig6,'-dpng','-r300');
%     figure; imshow(Output_DVStackImage_Transformed)
%     x_width=9 ;y_width=7; set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%     Fig6=strcat('Confocal','.png');
%     print(gcf,Fig6,'-dpng','-r300');
    filename = strcat('PWS','.png');
    t = Tiff(filename,'w');
    tagstruct.ImageLength = size(Output_cubeRMS_PWS_Transformed,1);
    tagstruct.ImageWidth = size(Output_cubeRMS_PWS_Transformed,2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.BitsPerSample = 32;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    t.setTag(tagstruct)
    t.write(Output_cubeRMS_PWS_Transformed);
    t.close;
    
    Confocal = single(sum(Output_DVStackImage_Transformed, 3));
    filename = strcat('Confocal','.png');
    t = Tiff(filename,'w');
    tagstruct.ImageLength = size(Output_DVStackImage_Transformed,1);
    tagstruct.ImageWidth = size(Output_DVStackImage_Transformed,2);
    tagstruct.Compression = Tiff.Compression.None;
    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.BitsPerSample = 32;
    t.setTag(tagstruct)
    t.write(Confocal);
    t.close;
    
    % Print Most Correlated index
    sprintf("Index: %d, %d", index_slice_maxcorrelation_begin, index_slice_maxcorrelation_end)
    
    % 4. Create D-map:
    RMS_map_PWS=Output_cubeRMS_PWS_Transformed;
    D_map_PWS=RMSMapToDmap(Output_cubeRMS_PWS_Transformed,OutputNucMask_PWS_Transformed,noise_PWS);     
    figure;imagesc(D_map_PWS);
    D_map_PWS=D_map_PWS.*single(OutputNucMask_PWS_Transformed);       
%     
%     % 5. plot image scale
%     dataRangeDensity = [1 max(D_map_PWS(:))];
%     figure;imagesc(D_map_PWS, dataRangeDensity);
%     colormap(jet);colorbar;axis off; set(gcf, 'PaperUnits', 'inches');
%     x_width=9 ;y_width=7; set(gcf, 'PaperPosition', [0 0 x_width y_width]);
%     Fig5_D_map_PWS=strcat('Fig5_D_map_PWS','.tif');
%     print(gcf,Fig5_D_map_PWS,'-dpng','-r300'); clear dataRangeDensity
end