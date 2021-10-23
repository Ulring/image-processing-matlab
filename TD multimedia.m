R    = [0 128 255;255 0 0;255 255 0]
G    = [0 128 255;0 255 0;255 0 255]
B    = [0 128 255;0 0 255;0 255 255]
IMG1 = cat(3, R, G, B) / 255;
figure
imshow(IMG1, 'InitialMagnification', 10000);
figure
imshow(rgb2gray(IMG1), 'InitialMagnification', 10000);
matrice_de_passage = [1/3 1/3 1/3;1/2 0 1/2;-1/4 1/2 -1/4]
IMG2 = cat(3, R * matrice_de_passage, G * matrice_de_passage, B * matrice_de_passage) / 255;
figure
imshow(IMG2, 'InitialMagnification', 10000);