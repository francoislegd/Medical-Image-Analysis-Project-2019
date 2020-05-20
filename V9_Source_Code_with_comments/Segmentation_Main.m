%% Project in Medical Image Analysis conducted by 
%  Olivier  AGBOHOUI
%  Hardik   DARBAR
%  François LEGRAND
%% Based on the original work of 
%  Marie-Pierre JOLLY
%% Subject of the project being
% Left ventricle epicardium and endocardium 
% sefmentation on 3D + T image set
% using training dataset from ACDC challenge for MICCAI
%% The following program will be our base code for the project
% it will be also reshaped as executable script with embedded GUI if
% possible
%% progam
clearvars
close all;
clc;

% unzip the target nifty file 
% filetoread = gunzip('patient001_4d.nii.gz');
% read and store the resulting nifty file
I = niftiread('patient002_4d.nii');

% One third of data have to be excluded
% exclude 4,6,11 maybe 8,10

% show dimensions of our image, I is 10 slices of 30 images in time
%whos I
[A, B, C, D] = size(I)
% A = x axis
% B = y Axis
% C = z axis
% D = time domain
%% Bloodpooldetection

% compute Fourier transform on temporal domain. 
% Main outputs are First rank harmonic on 10 slices (Harmonic_1)
% and DC_AVG which are the average images on the temporal cycle.
[FTimg, DC_AVG, Harmonic_1] = HeartROI_2(I, A, B, C, D);


%DC_AVG = int16(DC_AVG);
H = int16(FTimg);
% figure
% size(FTimg)
% displaying ft results timeframe by timeframe
% for slicenumber = 1:300
% %     for timeframe = 1:Time
%         %subplot(2,5,slicenumber);
%         imshow(H(:,:,slicenumber),[])
% %     end
% end
%% Test
% displaying slice by slice the temporal results
% There is noise artifacts in the middle of each sequences, maybe because
% of fourier harmonics
% We have to extract the first harmonics, it seems the second set of 10
% images
% coresponds to the targeted result

% % figure
% for slicenumber = 1:C
%     shifted = [];
%     for timeframe = 1:D
%         %subplot(2,5,slicenumber);
%         imshow (H(:,:,slicenumber,timeframe), [])
%     end
% end
%% ROI targeting

% In this function, we compute weighted centroids on Harmonic_1 images,
% apply distance mask around, take transformed images and recompute 5 times
% to obtain the best possible results, ideally centered on the most moving
% regions
[ROI_slices, H1_ROI, H1_post, Mask_ROI] = Heart_detector_3(Harmonic_1, DC_AVG);
%[ROI_slices_2, H1_ROI_2, H1_post_2] = Heart_detector_part2(H1_ROI, ROI_slices);


% The extracted ROI is refined in this function using active contour
% in-built function. the Extractd_ROI mask is reused for further
% processing, without superfluous areas which might disturb the program
[refined_mask,Extracted_ROI,H1_Extracted] = ROI_extractor(ROI_slices,H1_ROI, Mask_ROI);

%%


figure(16);
for i = 1:C
    subplot(2,6,i)
    imshow(DC_AVG(:,:,i),[]);
    hold on;
    visboundaries(refined_mask(:,:,i), 'LineWidth', 1, 'Color', 'y');
end

% The following function uses a dependancy found online at 
% https://fr.mathworks.com/matlabcentral/fileexchange/17933-polar-to-from-rectangular-transform-of-images
% We compute again weighted centroids on all slices using our refined ROI
% masked image, and use them as coordinate centers to perform polar
% coodinates
% transform. Althought originally intended for segmentation according to
% Jolly .et all, we intead use it to discard very dark values which can
% interfere with further thresholding to extract borders.
Pol_ROI = Centroid_pol(Extracted_ROI, DC_AVG);
%Pol_ROI = zeros(size(DC_AVG));
% Do polar transform wrt weighted centroid from H1_Extracted

figure(17);
for i = 1:C
    %Pol_ROI(:,:,i) = ImToPolar (DC_AVG(:,:,i), 0, 1, A, B);
    subplot(2,6,i)
    imshow(Pol_ROI(:,:,i),[]);
end

% Here, we use two functions to perform padding. The first one is used as
% basis for endocardial segmentation, the second one for epicardial
% segmentation.
ROI_padd = Niveleur(Extracted_ROI, Pol_ROI);
ROI_padd2 = Nivelage(ROI_slices, Pol_ROI, H1_ROI);
%ROI_padd2 = histeq(ROI_padd2);

% Binary images are created to compute Otsu's method on histogams of each
% level padded images
BW = zeros(size(ROI_padd));
BW2 = BW;
figure (18)
for i = 1:C
    subplot(2,6,i)
    %imshow(ROI_padd(:,:,i),[]);
    [counts,x] = imhist(uint8(ROI_padd(:,:,i)), 16);
    stem(x,counts)
    T1 = otsuthresh(counts);
    BW(:,:,i) = imbinarize(uint8(ROI_padd(:,:,i)),T1);
    imshow(BW(:,:,i));
end

figure (19)
for i = 1:C
    subplot(2,6,i)
    %imshow(ROI_padd2(:,:,i),[]);
    [counts,x] = imhist(histeq(uint8(ROI_padd2(:,:,i)), 16));
    stem(x,counts)
    T2 = otsuthresh(counts);
    BW2(:,:,i) = imbinarize(uint8(ROI_padd2(:,:,i)), T2);
    imshow(BW2(:,:,i),[]);
end

% There each functions to perform endocardium and epicardium segmentation
% are applied on the previously created binary images.
Endo = Endocardium(BW);
Epi = Epicardium(BW2);

% Here we put in evidence boundaries linked with various segmentations
figure(21);
for i = 1:C
    %Pol_ROI(:,:,i) = ImToPolar (DC_AVG(:,:,i), 0, 1, A, B);
    subplot(2,6,i)
    imshow(DC_AVG(:,:,i),[]);
%     hold on;
%     visboundaries(Mask_ROI(:,:,i), 'LineWidth', 0.2, 'Color', 'y');
    hold on;
    visboundaries(Endo(:,:,i), 'LineWidth', 0.2, 'Color', 'g');
    hold on;
    visboundaries(Epi(:,:,i), 'LineWidth', 0.2, 'Color', 'r');
end

figure(22);
for i = 1:C
    %Pol_ROI(:,:,i) = ImToPolar (DC_AVG(:,:,i), 0, 1, A, B);
    subplot(2,6,i)
    imshow(Epi(:,:,i),[]);
end
