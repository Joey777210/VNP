function [roi_2, roi_cluster] = find_roi(img_hazy, air_light, gamma)

[h,w,n_colors] = size(img_hazy);
if (n_colors ~= 3) % input verification
    error(['Non-Local Dehazing reuires an RGB image, while input ',...
        'has only ',num2str(n_colors),' dimensions']);
end

if ~exist('air_light','var') || isempty(air_light) || (numel(air_light)~=3)
    error('Dehazing on sphere requires an RGB airlight');
end

if ~exist('gamma','var') || isempty(gamma), gamma = 1; end

img_hazy = im2double(img_hazy);
img_hazy_corrected = img_hazy.^gamma; % radiometric correction

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

% 统计
cluster_stat = zeros(1, 1000);
[row, col] = size(cluster);
for i = 1:row
    for j = 1:col
        cluster_stat(cluster(i, j)) = cluster_stat(cluster(i, j)) + 1;
    end
end

% 求前三
[len, index] = sort(cluster_stat, 'descend');
disp(['长度： ', num2str(len(1:3)), ' index: ', num2str(index(1:3))]);

% 计算标准差
% Estimate radius as the maximal radius in each haze-line (Eq. (11))
K = accumarray(ind, radius(:), [n_points,1], @max);

radius_new = reshape( K(ind), h, w);

radius_1 = [];
radius_2 = [];
radius_3 = [];

dist_1 = [];
dist_2 = [];
dist_3 = [];

center_x = row;
center_y = 0;

for i = 1 : row
    for j = 1 : col
        if cluster(i, j) == index(1)
            radius_1 = [radius_1 radius(i, j)];
            dist_1 = [dist_1 sqrt((i-center_x)^2 + (j-center_y)^2)];
        elseif cluster(i, j) == index(2)
            radius_2 = [radius_2 radius(i, j)];
            dist_2 = [dist_2 sqrt((i-center_x)^2 + (j-center_y)^2)];
        elseif cluster(i, j) == index(3)
            radius_3 = [radius_3 radius(i, j)];
            dist_3 = [dist_3 sqrt((i-center_x)^2 + (j-center_y)^2)];
        end    
    end
end

standard_division(1,1) = kurtosis(dist_1);
standard_division(1,2) = kurtosis(dist_2);
standard_division(1,3) = kurtosis(dist_3);

disp(['标准差： ', num2str(standard_division)]);

[max_s_d, index2] = max(standard_division);

roi_cluster = index(index2);
% roi_cluster = 633;
roi_1 = zeros(row, col);

for i = 1:row
    for j = 1:col
        if cluster(i, j) == roi_cluster
            roi_1(i, j) = 1;
        end           
    end
end

%% 寻找roi第二阶段
% 1. 找半径最大的点的半径邻域[rmax - p, rmax]中的点
% 2. 这些点的最大连通域作为最终的roi

% 设置p为rmax的95%


radius_max_roi = 0;
for i = 1:h
    for j = 1 : w
        if cluster(i, j) == roi_cluster && radius_max_roi == 0
            radius_max_roi = radius_new(i, j);
        end
    end
end

bound = radius_max_roi - 0.07;
roi_2 = zeros(row, col);
count = 0;
for i = 1 : row
    for j = 1 : col
        if radius(i, j) >= bound && cluster(i, j) == roi_cluster
            roi_2(i, j) = 1;
            count = count + 1;
        end
    end
end

imLabel = bwlabel(roi_2);                %对各连通域进行标记
stats = regionprops(imLabel,'Area');    %求各连通域的大小
area = cat(1,stats.Area);
index = find(area == max(area));        %求最大连通域的索引
roi_2 = ismember(imLabel,index);          %获取最大连通域图像

%% 打印roi图片
cluster_img=zeros(row, col, 3);
for i = 1:row
    for j = 1 : col
        % 602的雾线非常牛逼
        % 124
        % 188 189
        if cluster(i,j) == roi_cluster
            cluster_img(i, j, 1)=0;
            cluster_img(i, j, 2)=50;
            cluster_img(i, j, 3)=0;
        else
            cluster(i,j)=20;
            cluster_img(i, j, 1)=255;
            cluster_img(i, j, 2)=255;
            cluster_img(i, j, 3)=255;
        end
        if roi_2(i,j) == 1
            cluster_img(i, j, 1)=255;
            cluster_img(i, j, 2)=0;
            cluster_img(i, j, 3)=187;
        end
    end
end

imwrite(cluster_img, 'roi_2.bmp')
