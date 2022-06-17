function [A_feat,B_feat,std_r,nr_d] = distance_msm_refinement(A_feat,B_feat,img1,img2,sigma_integration_filter)
    % Cut-out windows with distances
    % ------------------------------
    % Windows radius:
    wd_rad = ceil( sqrt( 12 * sigma_integration_filter ^ 2 + 1 ) / 2 ) + 1;
    % Estimated number of pixels:
    wd_n_px = (2 * wd_rad + 1) ^ 2;
    % Central pixels of each feature window:
    c_px_img1 = [ceil(A_feat.kp(:,1)) ceil(A_feat.kp(:,2))];
    c_px_img2 = [ceil(B_feat.kp(:,1)) ceil(B_feat.kp(:,2))];
    %Extract distances windows:
    wd_d_img1 = zeros(2 * wd_rad + 1,2 * wd_rad + 1,size(c_px_img1,1));
    wd_d_img2 = zeros(2 * wd_rad + 1,2 * wd_rad + 1,size(c_px_img2,1));

    for i = 1:size(c_px_img1,1)
        wd_d_img1(:,:,i) = img1.range(c_px_img1(i,1)-wd_rad : c_px_img1(i,1)+wd_rad,c_px_img1(i,2)-wd_rad:c_px_img1(i,2) + wd_rad);
    end

    for j = 1:size(c_px_img2,1)
        wd_d_img2(:,:,j) = img2.range(c_px_img2(j,1)-wd_rad:c_px_img2(j,1) + wd_rad,c_px_img2(j,2)-wd_rad:c_px_img2(j,2) + wd_rad);
    end
    % Number of pairs
    n_pairs = size(A_feat.kp,1);
    % Calculate distance differences:
    wd_d_diff = wd_d_img1 -wd_d_img2;
    % Transform them from 2D image to 1D vector
    d_diff = zeros(wd_n_px,size(wd_d_diff,3));
    for i = 1:n_pairs
        d_diff(:,i) = reshape(wd_d_diff(:,:,i),wd_n_px,1);
    end
    % Find outliers in distance differences (wrong correspondences)
    for i = 1:n_pairs
        test = ones(wd_n_px,1);
        if nnz(test) > 0
            test = isoutlier(d_diff(:,i));
            d_diff(test == 1,i) = NaN;
        end
        wd_d_diff(:,:,i) = reshape(d_diff(:,i),wd_rad*2+1,wd_rad*2+1);
    end

    % STD of distance differences (correct because all outliers removed)
    std_r = zeros(n_pairs,1);
    nr_d = zeros(n_pairs,1); % number of distances to compute nanstd
    for i = 1:n_pairs
        std_r(i) = nanstd(reshape(wd_d_diff(:,:,i),wd_n_px,1));
        nr_d(i) = sum(~isnan(reshape(wd_d_diff(:,:,i),wd_n_px,1)));
    end
    % Mark outliers in distance images as NaN
    for i = 1:n_pairs
        test = isnan(wd_d_diff(:,:,i));
        wd_d_current = wd_d_img1(:,:,i);
        wd_d_current(test) = NaN;
        wd_d_img1(:,:,i) = wd_d_current;

        wd_d_current = wd_d_img2(:,:,i);
        wd_d_current(test) = NaN;
        wd_d_img2(:,:,i) = wd_d_current;
        % can be optimized
    end
    clear wd_d_current test i;

    % Robust mean of distances for each feature (r_img1 & r_img2)
    r_img1 = zeros(n_pairs,1);
    r_img2 = zeros(n_pairs,1);

    for i = 1:n_pairs
        r_img1(i) = nanmean(reshape(wd_d_img1(:,:,i),wd_n_px,1));
        r_img2(i) = nanmean(reshape(wd_d_img2(:,:,i),wd_n_px,1));
    end
    % Storing new distances
    A_feat.rhv(:,1) = r_img1;
    B_feat.rhv(:,1) = r_img2;
end

