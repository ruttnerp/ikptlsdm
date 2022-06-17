function [feat] = featExtraction(I,parameter,type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% FUNCTION BEGINS
%% FEATURE DETECTION

switch type
    case 'foerstner'    % F??rstner Operator
        
        % Splitting intensity image according to vertical angle:
        swin_sz = deg2rad(15); % 15 [deg] Scan window size
        nr_swin = ceil((I.V(end,1)-I.V(1,1))/swin_sz); % Nr. of scan windows
        bottom = zeros(nr_swin,1); bottom(1) = 1;
        top = zeros(nr_swin,1); top(end) = numel(I.V(:,1));
        for i = 1:nr_swin-1
            bottom(i+1) = find(I.V(:,1) > I.V(1,1) + swin_sz*i ,1);
            top(i) = bottom(i+1) -1;
        end
        limits = [bottom,top];
        clear bottom top;
        
        % Defining tunable parameters
        par2 = linspace(parameter{2},parameter{2}+2,nr_swin);
        
        % Storing features of all images together
        pts_all.rc = []; pts_all.cov = [];
        for i = 1:nr_swin
            
            [~, corner, circle, ~]=ip_fop(I.int(limits(i,1):limits(i,2),:), ...
                'DETECTION_METHOD'         ,parameter{1},...
                'SIGMA_N'                  ,par2(i),...
                'DERIVATIVE_FILTER'        ,parameter{3},...
                'INTEGRATION_FILTER'       ,parameter{4},...
                'SIGMA_DERIVATIVE_FILTER'  ,parameter{5},...
                'SIGMA_INTEGRATION_FILTER' ,parameter{6},...
                'PRECISION_THRESHOLD'      ,parameter{7},...
                'ROUNDNESS_THRESHOLD'      ,parameter{8},...
                'SIGNIFICANCE_LEVEL'       ,parameter{9},...
                'VISUALIZATION'            ,parameter{10});
            
            % Extracting data from forstner_structures to simple structures
            ftype_forstner = 'all';
            % 'all' - corners and circles (default)
            % 'corners' - corners only
            % 'circles' - circles only
            
            [ pts ] = forstner_to_structure(corner,circle,ftype_forstner);
            pts.rc(:,1) = pts.rc(:,1) + limits(i,1)-1;
            
            pts_all.rc = [pts_all.rc;pts.rc];
            pts_all.cov = [pts_all.cov;pts.cov];
        end
        pts = pts_all;
        
    case 'SURF'     % SURF:Speeded Up Robust Features
        pts = detectSURFFeatures(I.int,...
            'MetricThreshold'   ,parameter{1},...
            'NumOctaves'        ,parameter{2},...
            'NumScaleLevels'    ,parameter{3});
        
    case 'KAZE'     % KAZE Detector
        pts = detectKAZEFeatures(I.int,...
            'Diffusion'     ,parameter{1},...
            'Threshold'     ,parameter{2},...
            'NumOctaves'    ,parameter{3},...
            'NumScaleLevels',parameter{4});
        
    case 'ORB'      % Oriented FAST and Rotated BRIEF
        pts = detectORBFeatures(I.int,...
            'ScaleFactor', 1.2,...
            'NumLevels', 8);
        
end

%% FEATURE EXTRACTION
if strcmp(type,'foerstner')
    [features, validPts] = extractFeatures(I.int,[pts.rc(:,2) pts.rc(:,1)],'Method','BRISK','Upright',true);
    features = features.Features;
else
    [features, validPts] = extractFeatures(I.int,pts,'Upright',true);
end

% extract row col location of keypoint
if strcmp(type,'foerstner')
    rc = [validPts(:,2),validPts(:,1)];
else
    rc = validPts.Location;
    %     rc = [validPts.Location(:,2), validPts.Location(:,1)];
end

% 3D coordinates of each keypoint
% POLAR
if strcmp(type,'foerstner')
    rhv = [interp2(I.range, rc(:,2), rc(:,1)),...
        interp2(I.H, rc(:,2), rc(:,1)), ...
        interp2(I.V, rc(:,2), rc(:,1))];
else
    rhv = [interp2(I.range, rc(:,1), rc(:,2)),...
        interp2(I.H, rc(:,1), rc(:,2)), ...
        interp2(I.V, rc(:,1), rc(:,2))];
end

% CARTESIAN
xyz = rhv2xyz(rhv);

%% OUTPUT
feat.features = features;
feat.pts = rc;
feat.rhv = rhv;
feat.xyz = xyz;

% remove features with NaN ranges
% they still have good angular measurements and could be saved, by
% changing interp2(I.range) with something (line 93)
test = isnan(feat.rhv(:,1));
feat.rhv(test,:) = []; feat.xyz(test,:) = [];
feat.features(test,:) = []; feat.pts(test,:) = [];

end

