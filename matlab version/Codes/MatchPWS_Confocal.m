%% Code Description: Match PWS with Confocal Slices using nuc mask correlation (DNA Density)

% Input: (1) PWS image and (2) PWS nuclear mask (hand-drawn); (3) Delta vision (DV)  DNA density Z-stack, and (4) DV nuclear mask (hand-drawn);
% Output: Transformed Images of all (4); (5) tranformation matrices; (6)
% slices with max correlation

function [Output_cubeRMS_PWS_Transformed,OutputNucMask_PWS_Transformed,Output_DVStackImage_Transformed,OutputNucMask_DV_Transformed,tform,index_slice_maxcorrelation_begin,index_slice_maxcorrelation_end]=MatchPWS_Confocal(cubeRms_PWS,BW_Nuc_PWS,Stack_DV,BW_Nuc_DV,Analysis_Cell_FilesPath,fileindexing)

%% A. Load Orginal PWS Image: 
distorted=logical(flipdim(rot90(BW_Nuc_PWS,1),1));    figure;imshow(distorted);       % Image A

%% B. Load Original DV confocal image stack:
BW_Nuc_DV= logical(double(bwareaopen(BW_Nuc_DV,7)));
se2 = strel('disk',1);
BW_Nuc_DV = imerode(BW_Nuc_DV,se2);

Labels_name=[dir('*_fixedpoints2.mat')];
load(Labels_name.name);
original=BW_Nuc_DV;                                  figure;imshow(original(:,:,4));  % Image B

%% C. Transform A to match B:
% documentation: https://www.mathworks.com/help/images/find-image-rotation-and-scale.html
correlation_matcher=zeros(1,size(Stack_DV,3));
cd(Analysis_Cell_FilesPath)
mkdir(['Cell' int2str(fileindexing)])
cd(['Cell' int2str(fileindexing)])
mkdir Mapping_DV_PWS
cd ("Mapping_DV_PWS\")

for i=1:size(original,3)
    % I selected points in fixedpoints.mat using:
    % cpselect(distorted,BW_Nuc_DV(:,:,4),movingPoints,fixedPoints); to get movingPoints1,fixedPoints1
    tform(i) = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
    tformInv = invert(tform(i));
    Tinv = tformInv.T;
    ss = Tinv(2,1);
    sc = Tinv(1,1);
    scale_recovered = sqrt(ss*ss + sc*sc);
    theta_recovered = atan2(ss,sc)*180/pi;
    Roriginal(:,:,i)= imref2d(size(original(:,:,i)));
    recovered(:,:,i) = imwarp(distorted,tform(i),'OutputView',Roriginal(:,:,i));

    montage({original(:,:,i),logical(recovered(:,:,i))})

    correlation_matcher(i)=corr2(original(:,:,i),recovered(:,:,i));
end
correlation_matcher(isnan(correlation_matcher))=0;

%% D. Plot correlation between PWS and DV per slice:
figure;plot(correlation_matcher(1:size(original,3)),'*r','MarkerSize',12,'linewidth',2);xlabel('Confocal slice #');ylabel('Correlation');title('Binary Mask Correlation between Confocal Slice and PWS Image');ylim([0.8 1]);set(gca,'FontSize',18) % Creates an axes and sets its FontSize to 18
set(gcf, 'PaperUnits', 'inches');
x_width=9 ;y_width=7;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
Fig1_Name_DV=strcat('Fig1_BinaryMaskCorrelation_PWS_DV','.tif');
print(gcf,Fig1_Name_DV,'-dpng','-r300');
close all;

%% E. Find slices with max correlation:
slice_maxcorrelation=find(correlation_matcher==max(correlation_matcher));
[~,Top4slice_maxcorrelation]= maxk(correlation_matcher,3);   % TOP 3 slices with max correlation
index_slice_maxcorrelation_begin=min(Top4slice_maxcorrelation);
index_slice_maxcorrelation_end=max(Top4slice_maxcorrelation);
disp('Max Correlation from DV confocal slices:');disp(mean(correlation_matcher(index_slice_maxcorrelation_begin:index_slice_maxcorrelation_end)));

%% F. Save output in directory of cell #:

Output_PWSImage_pretransformation=flip(rot90(cubeRms_PWS,1),1);

% Output PWS image and nuclear mask:
OutputNucMask_PWS_Transformed=recovered(:,:,index_slice_maxcorrelation_end-1); % OutputNucMask_PWS_Transformed=recovered(:,:,3);
Output_cubeRMS_PWS_Transformed=imwarp(Output_PWSImage_pretransformation,tform(index_slice_maxcorrelation_end-1),'OutputView',Roriginal(:,:,index_slice_maxcorrelation_end-1));
Output_cubeRMS_PWS_Transformed(OutputNucMask_PWS_Transformed==0)=0;

figure;imagesc(Output_cubeRMS_PWS_Transformed);colormap(jet);colorbar
set(gcf, 'PaperUnits', 'inches');
x_width=9 ;y_width=7;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); axis off;
Fig2_PWSTransformed_DV=strcat('Fig2_PWSTransformed','.tif');
print(gcf,Fig2_PWSTransformed_DV,'-dpng','-r300');
close all;

% Output DV image and nuclear mask:
OutputNucMask_DV_Transformed=(original(:,:,index_slice_maxcorrelation_begin:index_slice_maxcorrelation_end)); % figure;imshow(sum(OutputNucMask_DV_Transformed,3));
Output_DVStackImage_Transformed=Stack_DV(:,:,index_slice_maxcorrelation_begin:index_slice_maxcorrelation_end);
Output_DVStackImage_Transformed(OutputNucMask_DV_Transformed==0)=0;                                           % figure;imagesc(Output_DVStackImage_Transformed);

figure;imagesc(sum(Output_DVStackImage_Transformed,3));colormap(jet);
colorbar
set(gcf, 'PaperUnits', 'inches');axis off;
x_width=9 ;y_width=7;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
Fig3_DVTransformed=strcat('Fig3_DVTransformed','.tif');
print(gcf,Fig3_DVTransformed,'-dpng','-r300');
close all;

figure;imshowpair(imbinarize(sum(OutputNucMask_DV_Transformed,3)),OutputNucMask_PWS_Transformed);
set(gcf, 'PaperUnits', 'inches');
x_width=9 ;y_width=7;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
Fig4_Name_DV=strcat('Fig4_ImshowPair_PWS_DV','.tif');
print(gcf,Fig4_Name_DV,'-dpng','-r300');
end
