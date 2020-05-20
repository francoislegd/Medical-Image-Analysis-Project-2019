function [nivelee] = Nivelage(Extracted_ROI, DC_AVG, H1_ROI)

% A similar operation is done here, but instead we use the minimum value as
% limit. We make an interval of tolerance to remove all high values and
% lower ones to keep only the ones similar to left ventricle myocardium
[A, B, C] = size(DC_AVG)

%mini = max(DC_AVG,[],C);

for i = 1:C
    minimum = min(DC_AVG);
    for j = 1:A
        for k = 1:B
            if Extracted_ROI(j,k,i) > 2.6 * minimum(i)
                Extracted_ROI(j,k,i) = 1 - minimum(i);
            end
            if Extracted_ROI(j,k,i) == 0 
                Extracted_ROI(j,k,i) = 1 - 1.3*minimum(i);
            end
        end
    end
end

% We multiply the extracted region by 2 times the H1_ROI image, to see with
% more evidences the boundaries.
nivelee = Extracted_ROI .* (2*H1_ROI);
