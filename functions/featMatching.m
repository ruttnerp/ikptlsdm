function [A, B, metrics, indexPairs, cos_c ] = featMatching(A_feat, B_feat,T_type, T, type)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% FUNCTION BEGINS
% feature type based tuning
switch type
    case 'ORB'
    MatchingPara = [10, 0.6];
    case {'foerstner','SURF','KAZE'}
    %MatchingPara = [1, 0.6];  
    MatchingPara = [100, 1];
end

% DISTANCES BETWEEN POINTS
switch T_type
    case 'euclidean'    % euclidean 3D-distance between feature points
        D = pdist2(A_feat.xyz,B_feat.xyz,'euclidean');
        
    case 'angular' % Angle ("angular distance") between feature points
        cos_c = zeros(size(A_feat.xyz,1),size(B_feat.xyz,1));
        
        % Law of Spherical Cosine
        for i = 1:size(A_feat.xyz,1)
            cos_c(i,:) = cos(A_feat.rhv(i,3)).*cos(B_feat.rhv(:,3)) + sin(A_feat.rhv(i,3)).*sin(B_feat.rhv(:,3)).*cos(A_feat.rhv(i,2) - B_feat.rhv(:,2));
        end
        D = acos(cos_c);    % in radians
end

indexPairs = zeros(size(A_feat.xyz,1),2);
metrics = 1000*ones(size(A_feat.xyz,1),1);  % cant be zeros because zero equals perfect match
for i=1:size(A_feat.xyz,1)
    
    % find feature indices within distance threshold
    idx1 = find(D(i,:) < T); idx1 = idx1';
    
    switch type % slightly different data structure for different types
        case 'foerstner'
            %feat_check = ~isempty(B_feat.features.Features(idx1,:));
            feat_check = ~isempty(B_feat.features);
        case {'SURF','KAZE'}
            feat_check = ~isempty(B_feat.features);
    end
    
    if feat_check % Formerly B_feat.features -> probably not compatible with F?rstner or BRISK
        % Matching Corners
        switch type
            case 'foerstner'
                %feat_A = A_feat.features.Features(i,:);
                %feat_B = B_feat.features.Features(idx1,:);
                feat_A = A_feat.features(i,:);
                feat_B = B_feat.features(idx1,:);
            case {'SURF','KAZE'}
            	feat_A = A_feat.features(i,:);
                feat_B = B_feat.features(idx1,:);
        end
        [indexPair, metric] = matchFeatures(feat_A, feat_B,...
            'Method','Exhaustive',... Exhaustive Approximate
            'Unique',true,...
            'MatchThreshold',MatchingPara(1),...
            'MaxRatio',MatchingPara(2));
        
        if ~isempty(indexPair)
            indexPairs(i,:) = [i, idx1(indexPair(1,2))];
            metrics(i,:) = metric;
        end
    end
end

% clear Pairs of non Matches
indexPairs(indexPairs(:,1) == 0,:) = [];
metrics(metrics == 1000) = [];
%% Output
% Matched Keypoints
A.kp = A_feat.pts(indexPairs(:,1),:);
B.kp = B_feat.pts(indexPairs(:,2),:);
A.xyz = A_feat.xyz(indexPairs(:,1),:);
B.xyz = B_feat.xyz(indexPairs(:,2),:);
A.rhv = A_feat.rhv(indexPairs(:,1),:);
B.rhv = B_feat.rhv(indexPairs(:,2),:);
end

