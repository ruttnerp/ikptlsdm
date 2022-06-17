function [img_conc] = imCat(A,B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% concatenate two images of different size next to each other

% allocate image blank
if size(A,1) > size(B,1)
    img_conc = zeros(size(A,1),size(A,2)+size(B,2));
else
    img_conc = zeros(size(B,1),size(A,2)+size(B,2));
end

% allocating
 img_conc(1:size(A,1),1:size(A,2)) = A; 
 img_conc(1:size(B,1),size(A,2)+1:end) = B;
 
 img_conc = single(img_conc);
end

