function [epicard] = Epicardium(BW)

epicard = zeros(size(BW));
mask = zeros(size(BW));
[A, B, C] = size(BW);

% Due to scarcity of kept areas, we performed a closing to fuse them and
% eas process. They correspond to movement and left ventricle wall. 
se = strel('disk',5);
I = imclose(logical(BW), se);
% Then we keep only principal erea which corresponds roughly to left
% ventricle wall.
for i = 1:C
    epicard(:,:,i) = bwpropfilt(I(:,:,i),'FilledArea',1); % 'ConvexArea' 'FilledArea'    
end

% We compute centroid to further remove walls that do not belong to left
% ventricle with a smaller radius mask around a centroid midway between weighted one and
% real one (due to heart orientation).


Centroid_vector = zeros(2,C);
% for iter = 1:5
for i = 1:C
binaryImage = true(size(epicard(:,:,i)));
labeledImage = logical(binaryImage);
measurements1 = regionprops(labeledImage, epicard(:,:,i), 'WeightedCentroid');
measurements2 = regionprops(labeledImage, epicard(:,:,i), 'centroid');
centerOfMass = measurements1.WeightedCentroid;
center = measurements2.Centroid;
Centroid_vector(:,i) = (centerOfMass + center)/2;
% subplot(2,6,i);
% imshow(BW(:,:,i),[]);
hold on;
plot(Centroid_vector(1,i), Centroid_vector(2,i), 'r+', 'LineWidth', 1, 'MarkerSize', 16);
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(BW, 1);
    l = 1:size(BW, 2);   
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    subplot(2,6, k);
    imshow(Dist_hist,[]);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 35)
               mask(i,j,k) = 0;
            else 
               mask(i,j,k) = 1;
            end
            epicard(i,j,k) = epicard(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end

% Further processes are done to refine  area to extract. 
% Most importantly, we close and fill it using bwconvhull function, so we
% have a single boundary corresponding to left ventricle epicardium.
% Done
se2 = strel('disk',5);
epicard = imclose(epicard, se2);

for i = 1:C
    % https://fr.mathworks.com/matlabcentral/answers/243961-how-to-fill-the-black-regions-inside-the-white-region-imfill-is-not-working
    epicard(:,:,i) = bwconvhull(logical(epicard(:,:,i)));
    if i >= 8
        se3 = strel('square',3);
        epicard(:,:,i) = imerode(epicard(:,:,i), se3);
    end
end