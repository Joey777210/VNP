
haze_img = 'img/WILD0.5_input.jpg';
clear_img = 'img/WILD0.5.jpg';
k0 = 2.2518*10^(-5);

% 读图像
img_hazy = imread(haze_img);
imshow(img_hazy);
img_clear = imread(clear_img);
fid = fopen(['img/WILD1_params.txt'],'r');
[C] = textscan(fid,'%s %f');
fclose(fid);  
gamma = C{2}(1);

air_light = reshape(estimate_airlight(im2double(img_hazy).^(gamma)),1,1,3);

% 找ROI
[roi, roi_cluster] = find_roi(img_hazy, air_light, gamma);
disp(['roi_cluster: ', num2str(roi_cluster)])
imwrite(roi, 'roi.jpg')

% 计算透射率
transmission_hazy = cal_transmission_of_cluster(img_hazy, roi, roi_cluster, 1, air_light, gamma);
disp(['transmission_hazy: ', num2str(transmission_hazy)]);
% 去雾
% [img_dehazed, t, c] = non_local_dehazing(img_hazy, air_light, gamma);

% 计算透射率
transmission_clear = cal_transmission_of_cluster(img_clear, roi, roi_cluster, 2, air_light, gamma);
% disp(['transmission_clear: ', num2str(transmission_clear)]);

% 计算消光系数
k = k0 * (log(transmission_hazy(1,:))/log(transmission_clear(1,:)));
k = mean(k);
disp(['k: ', num2str(k)])

% 计算能见度
vis = 3 ./ k;
disp(['vis:', num2str(vis)]);

