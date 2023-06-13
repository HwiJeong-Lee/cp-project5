%% Homework Assignment5
% Student ID: 20233613
% 
% Student Name: 이휘정

%% 1. Initials (5 points)
% Load the light field image
image = imread('data/chessboard_lightfield.png');

lensletSize = 16; % each lenslet covers 16 x 16 pixels
u = lensletSize;
v = lensletSize;
s = size(image, 1) / lensletSize;
t = size(image, 2) / lensletSize;
c = 3;

% Initialize the 5-dimensional array
L = zeros(u, v, s, t, c);

for s_ = 1:s
    for t_ = 1:t
        for u_ = 1:u
            for v_ = 1:v
                for c_ = 1:c
                    L(u_, v_, s_, t_, c_) = image((s_-1)*u+u_, (t_-1)*v+v_, c_);
                end
            end
        end
    end
end

%% 2. Sub-aperture views (20 points)
% Initialize the mosaic image
mosaic = zeros(u*s, v*t, c);
mosaic = uint8(mosaic);

% Generate sub-aperture views and create the mosaic
for u_ = 1:u
    for v_ = 1:v
        mosaic(s*(u_-1)+1:s*(u_-1)+s, t*(v_-1)+1:t*(v_-1)+t, :) = L(u_, v_, :, :, :);
    end
end

% resize the image due to size limit of github
resized_mosaic = imresize(mosaic, 0.5, 'bilinear');
imwrite(resized_mosaic, 'results/mosaic.png');

%% 3. Refocusing and focal-stack generation (40 points)
% Generate focal stack for a range of depth values
maxUV = (lensletSize - 1) / 2;
u = (1:lensletSize) - 1 - maxUV;
v = (1:lensletSize) - 1 - maxUV;

depths = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2];
d = length(depths);

focal_stack = zeros(s, t, c, d);

for d_index = 1:d
    depth = depths(d_index);
    refocused_image = zeros(s, t, c);
    for i = 1:lensletSize
        for j = 1:lensletSize
            du = round(u(i) * depth) * -1;
            dv = round(v(j) * depth) * 1;
            img = squeeze(L(i, j, :, :, :));

            % Perform shifting
            img_shifted = zeros(size(img));
            for c_ = 1:c
                img_shifted(:, :, c_) = circshift(img(:, :, c_), [du, dv]);
            end
            refocused_image = refocused_image + img_shifted / lensletSize^2;    
        end
    end
    focal_stack(:, :, :, d_index) = refocused_image;
    imwrite(uint8(refocused_image), strcat('results/refocused_', num2str(depth), '.png'));
end

%% 4. All-focus image and depth from defocus (35 points)
all_focus_image = zeros(s, t, c);
depth_map = zeros(s, t);
w_sum = zeros(s, t);

% best parameters
std_1 = 3;
std_2 = 5;

for d_index = 1:d
    img = squeeze(focal_stack(:, :, :, d_index));
    img_xyz = rgb2xyz(img, 'ColorSpace', 'srgb');
    img_L = img(:,:,2);

    img_low = imgaussfilt(img_L, std_1);
    img_high = img_L - img_low;
    w_sharp = imgaussfilt(img_high.^2,std_2); 
    depth = depths(d_index);

    for c_ = 1:c
        all_focus_image(:, :, c_) = all_focus_image(:, :, c_) + img(:, :, c_) .* w_sharp;
    end

    depth_map = depth_map + w_sharp * depth;
    w_sum = w_sum + w_sharp;
end

all_focus_image = all_focus_image ./ w_sum;

max_depth = max(depths);
depth_map = max_depth - (depth_map ./ w_sum);

imwrite(uint8(all_focus_image), strcat('results/all_focus_', num2str(std_1), '_', num2str(std_2), '.png'));
imwrite(depth_map, strcat('results/depth_', num2str(std_1), '_', num2str(std_2),  '.png'));
