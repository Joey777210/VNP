function [transmission] = cal_transmission_of_cluster(img, roi, roi_cluster, flag, air_light, gamma)


%% Validate input
[h,w,n_colors] = size(img);
if (n_colors ~= 3) % input verification
    error(['Non-Local Dehazing reuires an RGB image, while input ',...
        'has only ',num2str(n_colors),' dimensions']);
end

if ~exist('air_light','var') || isempty(air_light) || (numel(air_light)~=3)
    error('Dehazing on sphere requires an RGB airlight');
end

if ~exist('gamma','var') || isempty(gamma), gamma = 1; end

img = im2double(img);
img_hazy_corrected = img.^gamma; % radiometric correction


%% Find Haze-lines
% Translate the coordinate system to be air_light-centric (Eq. (3))
dist_from_airlight = double(zeros(h,w,n_colors));
for color_idx=1:n_colors
    dist_from_airlight(:,:,color_idx) = img_hazy_corrected(:,:,color_idx) - air_light(:,:,color_idx);
end

% Calculate radius (Eq. (5))
radius = sqrt( dist_from_airlight(:,:,1).^2 + dist_from_airlight(:,:,2).^2 +dist_from_airlight(:,:,3).^2 );

% Cluster the pixels to haze-lines
% Use a KD-tree impementation for fast clustering according to their angles
dist_unit_radius = reshape(dist_from_airlight,[h*w,n_colors]);
dist_norm = sqrt(sum(dist_unit_radius.^2,2));
dist_unit_radius = bsxfun(@rdivide, dist_unit_radius, dist_norm);
n_points = 1000;
% load pre-calculated uniform tesselation of the unit-sphere
fid = fopen(['TR',num2str(n_points),'.txt']);
points = cell2mat(textscan(fid,'%f %f %f'));
fclose(fid);
mdl = KDTreeSearcher(points);
ind = knnsearch(mdl, dist_unit_radius);
cluster = reshape( ind, h, w);

% Estimating Initial Transmission
 
% Estimate radius as the maximal radius in each haze-line (Eq. (11))
K = accumarray(ind, radius(:), [n_points,1], @max);
radius_new = reshape( K(ind), h, w);

% 计算目标roi的透射率
radius_roi = [];
radius_cluster = [];

for i = 1:h
    for j = 1 : w
        if roi(i, j) == 1
            radius_roi = [radius_roi radius(i,j)];
        end
        if cluster(i, j) == roi_cluster
            radius_cluster = [radius_cluster radius(i, j)];
        end
    end
end

radius_max_roi = max(radius_cluster);

trans = radius_roi / radius_max_roi;
for i = 1 : length(trans)
    if trans(i) == 1
        trans(i) = 0.99999;
    end
end
transmission = trans;

% transmission(transmission==1)=0.9;
%% 关键的透射率选择问题
% x = find(trans == 1);
% trans(x) = [];
% transmission = mean(trans);
% 
% if transmission == 1 
%     transmission = 0.99999;
% end
