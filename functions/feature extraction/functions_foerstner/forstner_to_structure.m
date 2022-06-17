function [ pts ] = forstner_to_structure(corner,circle,ftype)

%% Extracting data from Förstner operator structures to arrays
% Input:
%       - corner [1x3] Structure containting corner feature points
%                  corner.r - row positions of the feature points
%                  corner.c - column positions of the feature points
%                  corner.cov   - covariance matrices of feature points
%       - circle [1x3] Structure containting circle feature points
%                  circle.r - row positions of the feature points
%                  circle.c - column positions of the feature points
%                  circle.cov   - covariance matrices of feature points
%       - ftype - which feature types you want to use
%               - 'all' - corners and circles (default)
%               - 'corners' - corners only
%               - 'circles' - circles only
% Output:
%       - pts [2x1] Structure containting feature points information
%       - pts.rc [nx2] - row and column of the image features
%       -pts.cov [nx3] - sigma_r^2 | sigma_c^2 | sigma_r*sigma_c
% author:   Tomislav Medic
% date:     10.12.2018

    if strcmp(ftype,'all')
        
        % Extracting data for Corners
        hv = size(corner,2); % help variable
        cor = zeros(hv,2);
        cor_cov = zeros(hv,3);
        for i = 1:hv
            cor(i,1) = corner(i).r;
            cor(i,2) = corner(i).c;
            cov = corner(i).cov;
            cor_cov(i,1) = cov(1,1);
            cor_cov(i,2) = cov(2,2);
            cor_cov(i,3) = cov(2,1);
        end
        
        % Extracting data for Circles
        hv = size(circle,2); % help variable
        cir = zeros(hv,2);
        cir_cov = zeros(hv,3);
        for i = 1:hv
            cir(i,1) = circle(i).r;
            cir(i,2) = circle(i).c;
            cov = circle(i).cov;
            cir_cov(i,1) = cov(1,1);
            cir_cov(i,2) = cov(2,2);
            cir_cov(i,3) = cov(2,1);
        end
        
        % Combining and storing data
        pts.rc = [cor;cir];
        pts.cov = [cor_cov;cir_cov];
        
    elseif strcmp(ftype,'corners')
        
        % Extracting data for Corners
        hv = size(corner,2); % help variable
        cor = zeros(hv,2);
        cor_cov = zeros(hv,3);
        for i = 1:hv
            cor(i,1) = corner(i).r;
            cor(i,2) = corner(i).c;
            cov = corner(i).cov;
            cor_cov(i,1) = cov(1,1);
            cor_cov(i,2) = cov(2,2);
            cor_cov(i,3) = cov(2,1);
        end
        % Storing data
        pts.rc = cor;
        pts.cov = cor_cov;
        
    elseif strcmp(ftype,'circles')
        
         % Extracting data for Circles
        hv = size(circle,2); % help variable
        cir = zeros(hv,2);
        cir_cov = zeros(hv,3);
        for i = 1:hv
            cir(i,1) = circle(i).r;
            cir(i,2) = circle(i).c;
            cov = circle(i).cov;
            cir_cov(i,1) = cov(1,1);
            cir_cov(i,2) = cov(2,2);
            cir_cov(i,3) = cov(2,1);
        end
        % Storing data
        pts.rc = cir;
        pts.cov = cir_cov;
    end
    
end

