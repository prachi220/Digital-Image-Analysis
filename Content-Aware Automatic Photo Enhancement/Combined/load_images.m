function imgs = load_images(target_image_name, gsource_image_name, source_image_name) 
    target_image = im2double(imread(target_image_name));
    gsource_image = im2double(imread(gsource_image_name));
    csource_image = im2double(imread(source_image_name));
    imgs.target_image = target_image;
    imgs.gsource_image = gsource_image;
    imgs.csource_image = csource_image;
end