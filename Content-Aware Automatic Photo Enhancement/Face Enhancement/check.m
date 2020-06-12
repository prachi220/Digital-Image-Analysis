
white = imread('white.jpg');
figure,
imshow(white);
black = imread('black.jpg');
figure,
imshow(black);
mask = imread('mask.png');
figure,
imshow(mask);
out = LaplacianBlend(black, white, mask);
figure,
imshow(out);
