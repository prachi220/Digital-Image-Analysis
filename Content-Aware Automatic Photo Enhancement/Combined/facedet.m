clear all
clc
%Detect objects using Viola-Jones Algorithm

%To detect Face
FDetect = vision.CascadeObjectDetector;

%Read the input image
input_image = imread('im6.png');
figure,
imshow(input_image);
%Returns Bounding Box values based on number of objects
BB = step(FDetect,input_image);

%coordinates for the face rectangle
x_face = BB(1);
y_face = BB(2);
width_face = BB(3);
height_face = BB(4);

cropped_face= imcrop(input_image,BB(1,:));
% figure,
% imshow(cropped_face);
% title('cropped face');

%Read the image, and capture the dimensions
VidImage = cropped_face;
height = size(VidImage,1);
width = size(VidImage,2);
%Initialize the output images
out = VidImage;
bin = zeros(height,width);
%Convert the image from RGB to YCbCr
img_ycbcr = rgb2ycbcr(VidImage);
Cb = img_ycbcr(:,:,2);
Cr = img_ycbcr(:,:,3);
%Detect Skin
[r,c,v] = find(Cb>=77 & Cb<=127 & Cr>=133 & Cr<=173);
numind = size(r,1);
%Mark Skin Pixels
for i=1:numind
    out(r(i),c(i),:) = [0 0 255];
    bin(r(i),c(i)) = 1;
end
binaryImage=im2bw(bin,graythresh(bin));
B = bwboundaries(binaryImage);
binaryImage = imfill(binaryImage, 'holes');
% Remove tiny regions.
binaryImage = bwareaopen(binaryImage, 5000);
%---------------------------------------------------------------------------
% Extract the largest area using ImageAnalyst's custom function ExtractNLargestBlobs().
biggestBlob = zeros(size(VidImage));
biggestBlob = ExtractNLargestBlobs(binaryImage, 1);
% figure,
% imshow(biggestBlob);
%--------------------------------------------------------------------------
[labeledImage, numberOfBlobs] = bwlabel(biggestBlob, 8);
% Get all the blob properties.
blobMeasurements = regionprops(labeledImage, 'BoundingBox','Area');
allBlobAreas = [blobMeasurements.Area];
% Loop through all blobs, putting up Bounding Box.
hold on; 
input_image_hsv = rgb2hsv(input_image);
input_image = hsv2rgb(input_image_hsv);

I = rgb2ycbcr(input_image);
P = I(:,:,1);
input_image_gray = im2double(I(:,:,1));
max_norm = max(input_image_gray(:));
input_image_gray = input_image_gray./max_norm;
base = wlsFilter(input_image_gray);
disp(size(base));
Iout = base;

% figure,
% imshow(base);
Iout = im2uint8(base);
Iout = im2double(Iout);
detail = imsubtract(input_image_gray,base);

%finding intersection -  not doubtful
S = zeros(size(cropped_face));
[rows,cols] = size(cropped_face);
[rowsB,colsB] = size(biggestBlob);
[rowsC,colsC] = size(Iout);
biggestBlob = im2uint8(biggestBlob);
S = im2uint8(S);
for col = 1:colsB
    for row = 1:rowsB
        if biggestBlob(row, col)> 0
            S(row, col) =  cropped_face(row,col);
        end
    end
end

%converted to gray image
S = rgb2gray(S);
% figure,
% % imshow(S);
% title('skin crop face');

%converted to uint8
Iout = im2uint8(Iout);
%applying sidelit correction
%get histogram of Iout/Base of the part of S
% smoothen only the skin region of the image
sub_imgae_Iout = zeros(size(width_face,height_face));
array = zeros(1,256,'double');
array1 = zeros(1,256,'double');
total_pixels=0;
for col = 1:colsB
    for row = 1:rowsB
        if biggestBlob(row, col)> 0 && y_face+row<=rowsC && x_face+col<=colsC  %this is the skin region
            val = Iout(y_face+row,x_face+col);
            val1 = cropped_face(row, col);
             array(val+1) = array(val+1) + 1 ;
             array1(val1+1) = array1(val1+1) + 1 ;
             total_pixels=total_pixels+1;
              
        end
    end
end
%smoothning
smooth_h2 = smooth(array);
% figure,
% hist(smooth_h2);
%now to check if the image is bimodal
[num_pixels,indices] = sort(smooth_h2,'descend');
max_num_pixels= num_pixels(1);
max_num_pixels_index= indices(1);
minimum_ind = 0;
minimum_val=0;
bright_mode_ind=0;
dark_mode_ind=0;
for i = 2:size(indices)
    next_max_index=indices(i);
    next_max=num_pixels(i);
    if(max_num_pixels_index>next_max_index)
        bright_mode_ind = max_num_pixels_index;
        dark_mode_ind = next_max_index;
    else
        dark_mode_ind = max_num_pixels_index;
        bright_mode_ind = next_max_index;
    end
    if(bright_mode_ind-dark_mode_ind >1)
        [min,min_indices] = sort(smooth_h2(dark_mode_ind + 1:bright_mode_ind - 1),'ascend'); %minimum between them
        minimum_val = min(1);
        minimum_ind = min_indices(1)+dark_mode_ind;
        if(minimum_val < 0.2*next_max)
            break;
        end
    end
end
% find the local minimum between them
if(~(minimum_ind==0 && minimum_ind==0))

        f = (bright_mode_ind - dark_mode_ind)/(1.8*(minimum_ind - dark_mode_ind));
        disp(f);
        % %create a scaling of Iout
        mask = zeros(size(Iout));
        black = Iout;
        A = Iout;

        for col = 1:colsB
            for row = 1:rowsB
                if biggestBlob(row, col)> 0 && y_face+row<=rowsC && x_face+col<=colsC && Iout(row,col)<=minimum_ind   %this is the skin region
                    A(row + y_face,col + x_face) =Iout(row + y_face,col + x_face) * f;
                    mask(row + y_face,col + x_face) = 255;

                end
            end
        end
        
        white = A;
        black = cat(3,black,black,black);
        white = cat(3,white,white,white);
        mask = cat(3,mask,mask,mask);
        Bout = LaplacianBlend(black, white, mask);
        Bout = LaplacianBlend(black, Bout, mask);
        Bout = LaplacianBlend(black, Bout, mask);
        Bout = LaplacianBlend(black, Bout, mask);

        A= rgb2gray(Bout);
        Iout= A;
end


%%%%%%%%%exposure correction%%%%%%%%%%%

percentile75_ind1 = prctile(array1,75);
disp(percentile75_ind1);
if (percentile75_ind1 < 120)
    f1 = (percentile75_ind1+120)/(2*percentile75_ind1);
    disp(f1);
    if(f1 >= 1 && f1 <=2)
        black = Iout;
        % %create a scaling of Iout
        for col = 1:colsB
            for row = 1:rowsB
                if biggestBlob(row, col)> 0 && y_face+row<=rowsC && x_face+col<=colsC && Iout(row,col)<=minimum_ind   %this is the skin region
                    A(row + y_face,col + x_face) =Iout(row + y_face,col + x_face) * f1;
                end
            end
        end
        white = A;
        black = cat(3,black,black,black);
        white = cat(3,white,white,white);
        mask = cat(3,mask,mask,mask);
        Bout = LaplacianBlend(black, white, mask);
        Bout = LaplacianBlend(black, Bout, mask);
        Bout = LaplacianBlend(black, Bout, mask);
        Bout = LaplacianBlend(black, Bout, mask);
%         figure,
%         imshow(Bout);
%         title('Bout exposed');
        A= rgb2gray(Bout);
        Iout= A;
    end
end

detail = im2uint8(detail);
% Iout = imadd(Iout,detail);
new_img = (Iout+detail).*max_norm;
I(:,:,1)= im2double(new_img);
I = ycbcr2rgb(I);
figure, imshow(I);
title('newimage');



