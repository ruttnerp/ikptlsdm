clc; clear all; addpath(genpath(cd)); addpath(genpath('../data'));
close all;

%% 2D-Feature-based Deformation Analysis
% Authors: Manuel Weber, Pia Ruttner, Tomislav Medic
% date: 16.06.2022

%% Variables and definitions
refscan = 0; % reference scan
n_scans = 25; % total number of scans
pillar_zmax = 148.9; % from bridge plans

% load correction values
corr_vals = load('correction_values.mat');

% initialize deformation vectors
deformation_vectors = cell(n_scans,2);

use_corr = false; % if correction values should be applied to scan values
get_pillar = true; % if the bounding box of the pillar has to be estimated
get_roi = false; % if the pillar area has to be defined in the images
filter_roi = true; % if the defined pillar area is used to filter matches

scans_to_cut = [0 5 13 23 24 25]; % some scan images have to be cut

% outlier detection options
image_based = true;
polar_based = true;
vector_based = true;
magnitude_based = true;
orientation_based = true;

Threshold_Factor = 2; % outlier threshold factor

% Feature Type
type = 'SURF';     % 'foerstner' 'SURF' 'KAZE' 'ORB'

% Parameters
switch type
    case 'foerstner'
        % Förstner
        parameter{1} = 'foerstner'; % method for optimal search window: 'foerstner' (default) or 'koethe'
        parameter{2} = 2; % standard deviation of (constant) image noise (default: 2.0)
        parameter{3} = 'gaussian2d'; %filter for gradient: 'gaussian2d'(default) oder 'gaussian1d'
        parameter{4} = 'gaussian'; % integration kernel: 'box_filt' (default) oder 'gaussian'
        parameter{5} = 1.5; % size of derivative filter (sigma) (default: 1.0)
        parameter{6} = 7; % size of integration filter (default: 1.41 * SIGMA_DIFF)
        parameter{7} = 1; % threshold for precision of points (default: 0.5 Pixel)
        parameter{8} = 0.3; % threshold for roundness (default: 0.3)
        parameter{9} = 0.999; % significance level for point classification (default: 0.999)
        parameter{10} = 'off'; % visualization on or off (default : 'off')

    case 'KAZE'
        parameter{1} = 'edge';      % 'Diffusion': Method to compute conductivity 'region' (default) 'sharpedge' 'edge'
        parameter{2} = 0.001;       % 'Threshold': Local extrema (default: 0.0001)
        parameter{3} = 3;           % 'NumOctaves': Multiscale detection factor (default: 3)
        parameter{4} = 4;           % 'NumScaleLevels': Scale levels (default: 4)

    case 'SURF'
        parameter{1} = 2000;        % 'MetricThreshold': Strongest feature threshold (default: 1000)
        parameter{2} = 3;           % 'NumOctaves': Number of octaves (default: 3)
        parameter{3} = 4;           % 'NumScaleLevels': Number of scale levels per octave (default: 4)
        parameter{6} = 7;           % size of integration filter (default: 1.41 * SIGMA_DIFF)

    case 'ORB'
        parameter{1} = 1.2;         % 'ScaleFactor': Scale factor for image decomposition (default: 1.2)
        parameter{2} = 8;           % 'NumLevels': Number of decomposition levels (default: 8)

end

%% Data Input
% Loading data:
% Image A - EPOCH 1
tic; % start measureing running time
load(sprintf('p50 data not registered/S%d_full_images.mat', refscan));
A_full = I;
A_full.int = im2uint8(mat2gray(A_full.int_raw));
clear I; A_full = rmfield(A_full,'int_raw');

% Handling scans to cut
A_full.int = A_full.int(140:end,:);
A_full.range = A_full.range(140:end,:);
A_full.H = A_full.H(140:end,:);
A_full.V = A_full.V(140:end,:);

if get_pillar
    % Visualizing point cloud
    % subsample ptCloud for segmentation and convert to ptCloud
    ptCloud = pointCloud(datasample(rhv2xyz([A_full.range(:),A_full.H(:),A_full.V(:)]),1000000));
    % segment pointcloud, parameter = min Dist to next cluster
    [labels,numClusters] = pcsegdist(ptCloud,0.2);
    % get correct label
    [GC, GR] = groupcounts(labels);
    GC = GC(2:end); GR = GR(2:end);
    [~, hv] = max(GC);
    label_pillar = GR(hv);
    % extract pillar point cloud
    pillar = ptCloud.Location(labels==label_pillar,:);
    % fit plane and estimate pillar normal direction
    pillar_plane = pcfitplane(pointCloud(pillar), 1);
    pillar_normal = pillar_plane.Normal;

    % bounding box of pillar (add some margin to make sure to include all points)
    bbp_margin = 0.1; % [m]
    bbp_min = min(pillar) + bbp_margin;
    bbp_max = max(pillar) + bbp_margin;
else % in case the pillar bounding box coordinates and pillar normal is known
    bbp_min = [78.2660 77.5115 6.0316];
    bbp_max = [84.1899 84.0868 141.1045];
    pillar_normal = [0.7910 0.6116 -0.0162];
end

% in case matching features should be filtered additionally by a region of
% interest (roi)
if filter_roi
    % determine roi
    if get_roi
        figure;imshow(A_full.int)
        roi = drawpolygon;
        save(sprintf("roi%d.mat", refscan),roi)
    % in case roi is already known and saved as .mat variable
    else
        roi = load(sprintf("roi%d.mat", refscan));
        roi = roi.roi;
    end
end

results.z_shift = pillar_zmax - max(pillar(:,3)); % get shift to refer all results to top of pillar

% Extract features of reference epoch/image (A)
% Image A
[A_feat_full_all] = featExtraction(A_full,parameter,type);

% filter for features on the pillar with pillar bounding box
test = ~((A_feat_full_all.xyz < bbp_max) & (A_feat_full_all.xyz > bbp_min));
test = test(:,1) | test(:,2) | test(:,3);
A_feat_full_all.features(test,:) = [];A_feat_full_all.pts(test,:) = [];
A_feat_full_all.rhv(test,:) = [];A_feat_full_all.xyz(test,:) = [];

if filter_roi
    % filter for features on the pillar with graphical roi
    test = ~(inpolygon(A_feat_full_all.pts(:,1),A_feat_full_all.pts(:,2), roi.Position(:,1), roi.Position(:,2)));
    A_feat_full_all.features(test,:) = [];A_feat_full_all.pts(test,:) = [];
    A_feat_full_all.rhv(test,:) = [];A_feat_full_all.xyz(test,:) = [];
end

% STORING RESULTS:
% Results stored as per "height segment" vector magnitudes
% Defining n segments
res = (8e-5)*180/pi;    % scanner resolution in degree
segment = 5;            % segments of 5 degree vertial
n_segments = floor(size(A_full.V,1)/(segment/res));
% Empty results:
% mean of each segment [height of segment, dist3d, dist1d, v(XYZ)]
results.mxyz = nan(n_segments,6,n_scans);
results.mrhv = nan(n_segments,6,n_scans);
results.fnr = nan(n_segments,n_scans);

results.vxyz = nan(n_segments,3,n_scans);
results.vrhv = nan(n_segments,3,n_scans);

results.stdvxyz = nan(n_segments,3,n_scans);
results.stdvrhv = nan(n_segments,3,n_scans);

results.v1dh = nan(n_segments,n_scans);
results.v2d = nan(n_segments,n_scans);
results.v3d = nan(n_segments,n_scans);
results.z = nan(n_segments,n_scans);
results.vangle = nan(n_segments,n_scans);

% Deformation magnitudes projected onto a plane
results.v3dp = nan(n_segments,n_scans);
results.v2dp = nan(n_segments,n_scans);

% All matched features (before outlier removal)
matches_noor = cell(n_segments,n_scans);
% All matched features (after outlier removal)
matches_or = cell(n_segments,n_scans);
% Some info out of it...

%% MAIN LOOP
% ----------

start_scan = refscan+1; % start reference scan + 1
for i=start_scan:n_scans
    fprintf('Scan Pair %i of %i\n',i,n_scans)

    % Image B
    nameB = strcat('S',num2str(i));
    load(strcat('p50 data not registered\',nameB,'_full_images.mat'));
    B_full = I;
    B_full.int = im2uint8(mat2gray(B_full.int_raw));
    clear I; B_full = rmfield(B_full,'int_raw');

    if ismember(i,scans_to_cut)
        B_full.int = B_full.int(140:end,:);
        B_full.range = B_full.range(140:end,:);
        B_full.H = B_full.H(140:end,:);
        B_full.V = B_full.V(140:end,:);
    end

    %% FEATURE EXTRACTION

    % Image A, start again with all features for feature matching and
    % outlier removal
    [A_feat_full] = A_feat_full_all;

    % Image B
    [B_feat_full] = featExtraction(B_full,parameter,type);

    % add corrections values
    if use_corr
        B_feat_full.rhv(:,2:3) = B_feat_full.rhv(:,2:3) + [corr_vals.h_corr(i+1), corr_vals.v_corr(i+1)];
        B_feat_full.xyz = rhv2xyz(B_feat_full.rhv);
    end


    %% FEATURE MATCHING
    %T_type = 'euclidean';  % 'euclidean' 'angular'
    %T = 0.1;   % euclidean threshold
    T_type = 'angular';
    T = deg2rad(0.05); % angular threshold
    [A_feat_matched, B_feat_matched, all_matches, metrics, ~] = featMatching(A_feat_full, B_feat_full, T_type, T, type);

    % Refining distance measurements
    if strcmp(type,'foerstner')
        [A_feat_matched,B_feat_matched,std_r,nr_d] = distance_msm_refinement(A_feat_matched,B_feat_matched,A_full,B_full,parameter{6});
    end

    % do outlier removal section based - IMPOSING LOCAL RIGIDITY
    start = 1;

    for k=1:n_segments

        if k < n_segments
            stop = k*ceil((segment/res));
            % Segment A
            A.int = A_full.int(start:stop,:);
            A.range = A_full.range(start:stop,:);
            A.H = A_full.H(start:stop,:);
            A.V = A_full.V(start:stop,:);
            A.xyz = rhv2xyz([A.range(:),A.H(:),A.V(:)]);
            %             A.xyz = A_full.xyz(start:stop,:);

            % Segment B
            B.int = B_full.int(start:stop,:);
            B.range = B_full.range(start:stop,:);
            B.H = B_full.H(start:stop,:);
            B.V = B_full.V(start:stop,:);
        else
            % Segment A
            A.int = A_full.int(start:end,:);
            A.range = A_full.range(start:end,:);
            A.H = A_full.H(start:end,:);
            A.V = A_full.V(start:end,:);
            A.xyz = rhv2xyz([A.range(:),A.H(:),A.V(:)]);
            % Segment B
            B.int = B_full.int(start:end,:);
            B.range = B_full.range(start:end,:);
            B.H = B_full.H(start:end,:);
            B.V = B_full.V(start:end,:);
        end

        % check again that only points are taken from the pillar
        test = ~((A.xyz < bbp_max) & (A.xyz > bbp_min));
        test = test(:,1) | test(:,2) | test(:,3);
        A.xyz(test,:) = [];

        if isempty(A.xyz)
            disp('No pillar points in segment')
            continue
        end

        % start with all features
        A_feat = A_feat_matched;
        B_feat = B_feat_matched;

        % take only features that are within section
        feat_mask_sec = ~((A_feat_matched.xyz(:,3) < max(A.xyz(:,3))) & (A_feat_matched.xyz(:,3) > min(A.xyz(:,3))));

        A_feat.kp(feat_mask_sec,:) = []; A_feat.xyz(feat_mask_sec,:) = []; A_feat.rhv(feat_mask_sec,:) = [];
        B_feat.kp(feat_mask_sec,:) = []; B_feat.xyz(feat_mask_sec,:) = []; B_feat.rhv(feat_mask_sec,:) = [];

        % Original number of matches:
        org_match_nr = size(A_feat.xyz,1);

        % Storing matches with no outlier removal
        results.matches_noor.xyz{k,i} = [A_feat_matched.xyz,B_feat_matched.xyz];
        results.matches_noor.kp{k,i} = [A_feat_matched.kp,B_feat_matched.kp];

        if org_match_nr < 3
            disp('not enough features in segment')
            A_feat.xyz = nan(2,3); B_feat.xyz = nan(2,3);
            start = start+ceil((segment/res));
            continue
        end

        %% OUTLIER REMOVAL

        % IMAGE BASED (MAD, row/column)
        if image_based
            dist = A_feat.kp - B_feat.kp;
            test = isoutlier(dist,'ThresholdFactor',Threshold_Factor);
            test = test(:,1) | test(:,2);
            A_feat.kp(test,:) = []; A_feat.xyz(test,:) = []; A_feat.rhv(test,:) = [];
            B_feat.kp(test,:) = []; B_feat.xyz(test,:) = []; B_feat.rhv(test,:) = [];
        end

        % 3D Polar-based (MAD, r,h,v)
        if polar_based
            dist = A_feat.rhv - B_feat.rhv;
            test = isoutlier(dist,'ThresholdFactor',Threshold_Factor);
            test = test(:,1) | test(:,2) | test(:,3);
            A_feat.kp(test,:) = []; A_feat.xyz(test,:) = []; A_feat.rhv(test,:) = [];
            B_feat.kp(test,:) = []; B_feat.xyz(test,:) = []; B_feat.rhv(test,:) = [];
        end

        % 3D VECTOR BASED (MAD, x,y,z)
        if vector_based
            dist = B_feat.xyz-A_feat.xyz;
            test = isoutlier(dist,'ThresholdFactor',Threshold_Factor);
            test = test(:,1) | test(:,2) | test(:,3);
            A_feat.kp(test,:) = []; A_feat.xyz(test,:) = []; A_feat.rhv(test,:) = [];
            B_feat.kp(test,:) = []; B_feat.xyz(test,:) = []; B_feat.rhv(test,:) = [];
        end

        % 3D Magnitude
        if magnitude_based
            dist = B_feat.xyz-A_feat.xyz;
            dist = vecnorm(dist,2,2);
            test = isoutlier(dist,'ThresholdFactor',Threshold_Factor);
            A_feat.kp(test,:) = []; A_feat.xyz(test,:) = []; A_feat.rhv(test,:) = [];
            B_feat.kp(test,:) = []; B_feat.xyz(test,:) = []; B_feat.rhv(test,:) = [];
        end

        % 3D Orientation
        if orientation_based
            v = B_feat.xyz-A_feat.xyz;
            v_median = repmat(median(v,1),size(v,1),1);
            angle = acos( dot(v_median,v,2) ./ ( vecnorm(v_median,2,2) .* vecnorm(v,2,2) ) );
            c = -1/(sqrt(2)*erfcinv(3/2)); % correction factor between MAD and std
            angle_std = c* median(abs(angle)); % std around zero
            test = angle > Threshold_Factor * angle_std;
            A_feat.kp(test,:) = []; A_feat.xyz(test,:) = []; A_feat.rhv(test,:) = [];
            B_feat.kp(test,:) = []; B_feat.xyz(test,:) = []; B_feat.rhv(test,:) = [];
        end

        v  = B_feat.xyz-A_feat.xyz;

        % make sure more than 10 features are left
        if size(A_feat.xyz,1) < 10 % if not enough matches -> want 10 or more after outlier removal
            disp('after outlier removal not enough features')
            A_feat.xyz = nan(2,3); B_feat.xyz = nan(2,3);
            start = start+ceil((segment/res));
            continue
        end

        % Store filtered feature matches:
        results.matches_or.xyz{k,i} = [A_feat.xyz,B_feat.xyz];
        results.matches_or.kp{k,i} = [A_feat.kp,B_feat.kp];

        % Store results:
        r = median(A_feat.rhv(:,1));
        h = median(B_feat.rhv(:,2) - A_feat.rhv(:,2));
        v1dh = (r * sin(h))*1000;

        % filter/check if median vector is significantly off zero
        v = B_feat.xyz - A_feat.xyz;
        v_mean = mean(v);

        % store result values
        results.mxyz(k,:,i) = [mean(A_feat.xyz,1),mean(B_feat.xyz,1)];
        results.mrhv(k,:,i) = [mean(A_feat.rhv,1),mean(B_feat.rhv,1)];
        results.fnr(k,i) = size(A_feat.xyz,1);

        results.vxyz(k,:,i) = v_mean;
        results.vrhv(k,:,i) = mean(B_feat.rhv - A_feat.rhv,1);

        results.stdvxyz(k,:,i) = std(v);
        results.stdvrhv(k,:,i) = std(B_feat.rhv - A_feat.rhv,1,1);

        results.v1dh(k,i) = v1dh;
        results.v3d(k,i) = norm(v_mean);
        results.v2d(k,i) = norm(v_mean(1:2));
        results.z(k,i) = 0.5*(min(A.xyz(:,3)) + max(A.xyz(:,3)));
        results.vangle(k,i) = median(A_feat.rhv(:,3));

        % Projecting 2d on a plane
        n_westwall = [0.617713	-0.78633	0.012437];
        v3dp = dot(v_mean,n_westwall)*n_westwall;

        % get angle between v3dp and normal to check if they point in same
        % direction
        angle = atan2(norm(cross(v3dp,n_westwall)), dot(v3dp,n_westwall));
        angle_th = 1e-3;
        if abs(angle) < angle_th
            direction = 1;
        elseif abs(angle-pi) < angle_th
            direction = -1;
        else
            direction = nan;
        end
        v2dp = norm(v3dp(1:2)) * direction;
        v3dp = norm(v3dp) * direction;

        results.v3dp(k,i) = v3dp;
        results.v2dp(k,i) = v2dp;

        % Help variable for looping through segments
        start = start+ceil((segment/res));

    end

    t = toc;
    fprintf('(in %f sec) \n',t)
end

%% data preparation for further analysis
% load correction values from corresponding M3C2 analysis
results.v3dp_corr = nan(n_segments,n_scans);

if refscan == 5
    resultsP50M5 = load('resultsP50M5.mat');
    resultsP50M5 = resultsP50M5.resultsP50M5;
    corr_vals = resultsP50M5.m3c2(end,:); % take last row of m3c2 as correction
elseif refscan == 0
    resultsP50M0 = load('resultsP50M0.mat');
    resultsP50M0 = resultsP50M0.resultsP50M0;
    corr_vals = resultsP50M0.m3c2(end,:);
end

cols = 1+refscan:refscan+numel(corr_vals);
results.v3dp_corr(:,cols) = -(results.v3dp(:,cols) - corr_vals);

% do not use last scan
if refscan == 0
    results.v3dp_corr(:,end) = nan;
end

% shift z and use same window as for v2sp
results.z_corr = results.z + results.z_shift;

%% reshape for curve fitting
z_corr_all = reshape(results.z_corr, [],1);
v3dp_corr_all = reshape(results.v3dp_corr, [],1);

%time
load scantimesP50;
scantimes = scantimesP50{2};
time_s = seconds(scantimes(2:end)-scantimes(refscan+1));
time_all = reshape(repmat(time_s', size(results.z_corr,1), 1), [],1);

[results.curvefit, ~] = createFit(z_corr_all, time_all, v3dp_corr_all);

% manually do curvefit with curve fitting toolbox and save all outputs as
% results.curvefit.fittedmodel
% results.curvefit.goodness
% results.curvefit.output

%% PLOT
% -----------

plot_nr = 1:25;
cmap = turbo(numel(plot_nr));
figure;
grid;
hold on;
for p=plot_nr
    plot(results.v3dp_corr(:,p), results.z_corr(:,p), '*-', 'Color', cmap(p,:))
end


load scantimesP50
scantimes = scantimesP50{2}(arrayfun(@(x) find(scantimesP50{1}==x,1),plot_nr));
legend_hm = [];
for i=1:numel(scantimes)
    [h,m] = hms(scantimes(i));
    legend_hm = [legend_hm, sprintf("%02d:%02d", [h,m])];
end
legend(legend_hm);


% plot_corr5_fit
plot_nr = 5:25;
% cmap = turbo(numel(plot_nr));
cmap = turbo(plot_nr(end));
figure;
hold on;
for p=plot_nr
    plot(results.fit_2d_corr5_2d(:,p), results.z(:,p), '*-', 'Color', cmap(p,:))
end
load scantimesP50
scantimes = scantimesP50{2}(arrayfun(@(x) find(scantimesP50{1}==x,1),plot_nr));
legend_hm = [];
for i=1:numel(scantimes)
    [h,m] = hms(scantimes(i));
    legend_hm = [legend_hm, sprintf("%02d:%02d", [h,m])];
end
legend(legend_hm);
grid;

clearvars -except results matches_noor matches_or;
% save foerstner_results.mat;

%% more plots
% plot timeseries

load scantimesP50
scantimes = scantimesP50{2}(arrayfun(@(x) find(scantimesP50{1}==x,1),plot_nr));
plot_seg = 1:size(results.z,1);
cmap = turbo(numel(plot_seg));
figure;
hold on;
for s=plot_seg
    plot(scantimes, results.v2dp(s,:), 'Color', cmap(s,:))
end
mean_h = mean(results.z,2);
legend_h = [];
for i=1:numel(mean_h)
    [h,m] = hms(scantimes(i));
    legend_h = [legend_h, sprintf("%.0fm",mean_h(i))];
end
legend(legend_h);
grid;

%% plot matches on pillar (before/after outlier removal)

% get all matches
epoch_B = 18;
matches_or_kp_A = [];
matches_or_kp_B = [];
matches_noor_kp_A = [];
matches_noor_kp_B = [];
for i=1:size(matches_or,1)
    kp_segment_or = results.matches_or.kp{i,epoch_B};
    kp_segment_noor = results.matches_noor.kp{i,epoch_B};
    if ~isempty(kp_segment_or)
        matches_or_kp_A = [matches_or_kp_A; kp_segment_or(:,1:2)];
        matches_or_kp_B = [matches_or_kp_B; kp_segment_or(:,3:4)];
    end
    if ~isempty(kp_segment_noor)
        matches_noor_kp_A = [matches_noor_kp_A; kp_segment_noor(:,1:2)];
        matches_noor_kp_B = [matches_noor_kp_B; kp_segment_noor(:,3:4)];
    end
end

figure; ax=axes; showMatchedFeatures(A_full.int, B_full.int, matches_noor_kp_A, matches_noor_kp_B,'montage','Parent',ax);
figure; ax=axes; showMatchedFeatures(A_full.int, B_full.int, matches_or_kp_A, matches_or_kp_B,'montage','Parent',ax);

