%%
clc;
%%
tic;
GRAPH = false;
SAVE = false;

%% load images
img_name = '1.jpg';
imgs = load_images('s1.png', 'csource.jpg', 'source.jpg');
target = {};
gsource = {};
csource = {};
target.image = imgs.target_image;
final_image = target.image;
gsource.image = imgs.gsource_image;
csource.image = imgs.csource_image;
tic;
figure,
imshow(target.image);
title('input');
%% to LAB color space

target_lab = rgb2lab(target.image);
target.image = rgb2gray(target.image);
skymask = target_lab(:,:,3);
new_target= ones(size(skymask));
[rowT, colT] = size(skymask);

for row=1:rowT
    for col=1:colT
        if (skymask(row,col)<0)
            pixel = target.image(row,col);
            new_target(row,col)=pixel;
        end
    end
end

target.image = new_target;
csource.lab = rgb2lab(csource.image);
if ndims(gsource.image) == 3
    gsource.lab = rgb2lab(gsource.image);
else
    gsource.lab = gsource.image;
end
if ndims(target.image) == 3
    target.lab = rgb2gray(target.image); 
else
    target.lab = target.image;
end

%% map luminance to target luminance
csource.luminance = luminance_remap(csource.lab, target.lab);
gsource.luminance = luminance_remap(gsource.lab, target.lab);
target.luminance = target.lab;
% pixel values are luminance

%% 
transferred = image_colorization_jitter_sampling(target, csource, GRAPH);
new_name = strcat('jitter_', img_name);
%%
new_image = lab2rgb(transferred);
for row=1:rowT
    for col=1:colT
        if (skymask(row,col)<0)
            pixel = new_image(row,col,:);
            final_image(row,col,:)=pixel;
        end
    end
end
toc;
%%
figure,
imshow(final_image);
%% save images
success = save_image(new_image, new_name, SAVE);
toc;