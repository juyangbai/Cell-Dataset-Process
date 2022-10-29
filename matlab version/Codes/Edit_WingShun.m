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
    
    Confocal = single(sum(Output_DVStackImage_Transformed,3));
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