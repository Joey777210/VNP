
win_size = 3;

image = double(imread('img/WILD0.5_input.jpg'))/255;

% image = imresize(image, 0.1);

[m, n, ~] = size(image);

pad_size = floor(win_size/2);

padded_image = padarray(image, [pad_size pad_size], Inf);

dark_channel = zeros(m, n); 

for j = 1 : m
    for i = 1 : n
        patch = padded_image(j : j + (win_size-1), i : i + (win_size-1), :);
        dark_channel(j,i) = min(patch(:));
     end
end
figure,imshow(dark_channel)

%% 找目标域roi
dark_channel = dark_channel * 255;
roi = zeros(m, n);
min_dark_channel = min(dark_channel(:));
s = 15;

for i = 1 : m
    for j = 1 : n
        if dark_channel(i, j) >= min_dark_channel && dark_channel(i, j) <= (min_dark_channel + s)
            roi(i, j) = 1;
        end
    end
end

imLabel = bwlabel(roi);                %对各连通域进行标记
stats = regionprops(imLabel,'Area');    %求各连通域的大小
area = cat(1,stats.Area);
index = find(area == max(area));        %求最大连通域的索引
roi = ismember(imLabel,index);          %获取最大连通域图像

%% 打印roi图片
dark_channel = dark_channel / 255;
dark_channel_new = zeros(m, n, 3);
dark_channel_new(:, :, 1) = dark_channel;
dark_channel_new(:, :, 2) = dark_channel;
dark_channel_new(:, :, 3) = dark_channel;
dark_channel = dark_channel_new;

for i = 1:m
    for j = 1 : n
        if roi(i,j) == 1
            dark_channel(i, j, 1)=255;
            dark_channel(i, j, 2)=0;
            dark_channel(i, j, 3)=187;
        end
    end
end
imwrite(dark_channel, 'roi_DCP.jpg')
figure,imshow(dark_channel)
        
