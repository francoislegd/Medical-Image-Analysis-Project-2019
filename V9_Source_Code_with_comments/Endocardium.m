function [endocardium] = Endocardium(BW, DC_AVG)

endocardium = zeros(size(BW));
epicard = zeros(size(BW));
mask = zeros(size(BW));
[A, B, C] = size(BW);

A2 = BW;
A3 = BW;

% Centroids are computed and refined once on binary  slices, so a new
% distance mask is computed with lower radius, to remove all unecessary
% areas around the centroids, which in the ideal case must located in the
% left vntricle bloodpool.

figure(20);
Centroid_vector = zeros(2,C);
% for iter = 1:5
for i = 1:C
binaryImage = true(size(BW(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, BW(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
subplot(2,6,i);
imshow(BW(:,:,i),[]);
hold on;
plot(centerOfMass(1), centerOfMass(2), 'r+', 'LineWidth', 1, 'MarkerSize', 16);
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
            if (mask(i,j,k) > 10)
               mask(i,j,k) = 0;
            else 
               mask(i,j,k) = 1;
            end
            A2(i,j,k) = BW(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end
Centroid_vector = (Centroid_vector)';
% end
for i = 1:C
binaryImage = true(size(A2(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, A2(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
subplot(2,6,i);
imshow(A2(:,:,i),[]);
hold on;
plot(centerOfMass(1), centerOfMass(2), 'r+', 'LineWidth', 1, 'MarkerSize', 16);
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';

for k = 1:C    
    m = 1:size(A2, 1);
    l = 1:size(A2, 2);   
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    subplot(2,6, k);
    imshow(Dist_hist,[]);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 10)
               mask(i,j,k) = 0;
            else 
               mask(i,j,k) = 1;
            end
            A3(i,j,k) = A2(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end
Centroid_vector = (Centroid_vector)';
for i = 1:C
binaryImage = true(size(A3(:,:,i)));
labeledImage = logical(binaryImage);
measurements = regionprops(labeledImage, A3(:,:,i), 'WeightedCentroid');
centerOfMass = measurements.WeightedCentroid;
Centroid_vector(:,i) = centerOfMass;
subplot(2,6,i);
imshow(A3(:,:,i),[]);
hold on;
plot(centerOfMass(1), centerOfMass(2), 'r+', 'LineWidth', 1, 'MarkerSize', 16);
%     https://fr.mathworks.com/matlabcentral/answers/368901-calculate-center-of-image
%     C = round([A B]/2) ;
%     plot(C(1),C(2),'*r');
%     plot(stat(i).Centroid(1),stat(i).Centroid(2),'ro');
end
Centroid_vector = (Centroid_vector)';
for k = 1:C    
    m = 1:size(A3, 1);
    l = 1:size(A3, 2);
    % maybe an if statement
    Dist_hist = sqrt((m.' - Centroid_vector(k,2)) .^ 2 + (l - Centroid_vector(k,1)) .^ 2);%./Storage1(y,x,k);
    subplot(2,6, k);
    imshow(Dist_hist,[]);
    mask(:,:,k) = Dist_hist;
    
end

for i = 1:A
    for j = 1:B
        for k = 1:C
            if (mask(i,j,k) > 45)
                mask(i,j,k) = 0;
            else 
                mask(i,j,k) = 1;
            end
            BW(i,j,k) = BW(i,j,k) * mask(i,j,k);
            %Storage_num5(i,j,k) = mask(i,j,k) * Storage_num4(i,j,k);
        end
    end
end 

% There, we perform an opening to remove papillar muscles, and use the
% function bwpropfilt with filled area as argument as the left ventricle is
% the most compact one in theory. This keep the bigesst compact area in the image, while removing the others. Sometimes it instead targets right
% ventricle, in particular when it is bigger tha the left one.
se = strel('disk',5);
I = imopen(logical(BW), se);
for i = 1:C
    endocardium(:,:,i) = bwpropfilt(I(:,:,i),'FilledArea',1); % 'ConvexArea' 'FilledArea'    
%     if i  <=8
%         endocardium(:,:,i) = bwpropfilt(I(:,:,i),'EquivDiameter',[25 35]); % 'ConvexArea' 'FilledArea'
%     else
%         endocardium(:,:,i) = bwpropfilt(I(:,:,i),'EquivDiameter',[20 30]);
%     %endocardium(:,:,i) = imclose(BW(:,:,i), se);
%     end   
end




