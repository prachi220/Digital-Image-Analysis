input_image = imread('samplepics/kid.png');
figure,
imshow(input_image);
title('input image');

 map = gbvs(input_image); % map.master_map contains the actual saliency map 
 map_itti = ittikochmap(input_image); % map_itti.master_map contains the actual saliency map 
 hm = heatmap_overlay(input_image, map_itti.master_map_resized);
 
input_img_lab = rgb2lab(input_image);
lab_l = input_img_lab(:,:,1);
lab_a = input_img_lab(:,:,2);
lab_b = input_img_lab(:,:,3);

mask = zeros(size(input_image(:,:,1)));

%split the image into dark and bright

[rowsL,colsL] = size(lab_l);
dark_matrix = ones(size(lab_l)).*51;
bright_matrix = zeros(size(lab_l));
total_dark = 0;
total_bright = 0;
for col = 1:colsL
    for row = 1:rowsL
%         
        if lab_l(row,col) < 51
                dark_matrix(row,col) = lab_l(row,col);
                total_dark = total_dark + 1;   
        else
            bright_matrix(row,col) = lab_l(row,col);
            total_bright = total_bright + 1;
        end
    end
end

dark_array = zeros(1,total_dark, 'double');
bright_array = zeros(1,total_bright, 'double');
dark_counter=1;
bright_counter=1;

for col = 1:colsL
    for row = 1:rowsL
        if lab_l(row,col) < 50
                dark_array(dark_counter) = lab_l(row,col);
                dark_counter=dark_counter+1;  
        else
             bright_array(bright_counter) = lab_l(row,col);
             bright_counter=bright_counter+1;
        end
    end
end

%calculate the factor
percentile35 = prctile(bright_array,35);
percentile95 = prctile(dark_array,95);
disp(percentile35);
disp(percentile95);
f_sal = min(2,(percentile35/percentile95));
disp(f_sal);

%apply wlsfilter to dark

lab_l_double = im2double(lab_l);
base = wlsFilter(lab_l_double);
detail = lab_l_double - base;

for col = 1:colsL
    for row = 1:rowsL
        Mi = hm(row,col);
        Bi = base(row,col);
        if(dark_matrix(row,col) <51)
            base(row,col) = (f_sal*Mi*Bi +((1-Mi)*Bi));
        end
    end
end
out_image = base+detail;
input_img_lab(:,:,1)= out_image;
final_image=lab2rgb(input_img_lab);
figure,
imshow(final_image);
title('outimage');