% 1) Créer un fichier .m qui lit et affiche l’image ’gravier.png’.
[img, map] = imread('gravier.png');
img_rgb = ind2rgb(img, map);
% figure();
% imshow(img_rgb)
[img, map] = imread('gravier2.png');
img_rgb = ind2rgb(img, map);
% figure();
% imshow(img_rgb)
info = imfinfo('gravier2.png');
% display(info);

% 2) Convertir l’image ’gravier.png’ en une image d’intensité et donner les caractéristiques de cette image.
[img, map] = imread('gravier2.png');
img_gray = ind2gray(img, map);
imwrite(img_gray,'gravier2_gris.png');
info = imfinfo('gravier2_gris.png');
% figure();
% imshow(img_rgb);
% display(info);

% 3) Tester les outils imtool et impixelinfo sur l’image couleur et l’image d’intensité obtenue précédemment.
% [X,map] = imread('gravier2.png');
% imtool(X,map)
% h = imshow(img_rgb);
% hp = impixelinfo;
% [X,map] = imread('gravier2_gris.png');
% imtool(X,map)
% h = imshow(img_rgb);
% hp = impixelinfo;

% 4) Extraire et afficher les images correspondant aux composantes R, G et B.
% [img, map] = imread('gravier2.png');
% img = ind2rgb(img, map);
% red = img(:,:,1); % Red channel
% green = img(:,:,2); % Green channel
% blue = img(:,:,3); % Blue channel
% % Que peut-on observer ? : R,G,B on des valeurs entre 0 et 1.

% 5) En utilisant les fonctions subplot et title, organiser le visuel.
% a = zeros(size(img, 1), size(img, 2));
% just_red = cat(3, red, a, a);
% just_green = cat(3, a, green, a);
% just_blue = cat(3, a, a, blue);
% back_to_original_img = cat(3, red, green, blue);
% subplot(3, 3, 5);
% imshow(img);
% fontSize = 10;
% title('Image couleur originale', 'FontSize', fontSize)
% subplot(3, 3, 7);
% imshow(just_red);
% title('R', 'FontSize', fontSize)
% subplot(3, 3, 8);
% imshow(just_green)
% title('G', 'FontSize', fontSize)
% subplot(3, 3, 9);
% imshow(just_blue);
% title('B', 'FontSize', fontSize)
% subplot(3, 3, 2);
% imshow(img_gray);
% title('Image d’intensité', 'FontSize', fontSize)
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1])
% set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')

% 6) Transformer l’image couleur dans l’espace HSV puis Lab et observer le résultat obtenu.
% hsv = rgb2hsv(img_rgb);
% figure();
% imshow(hsv);
% lab = rgb2lab(img_rgb);
% figure();
% imshow(lab);

% 7) La fonction imhist permet de calculer l’histogramme d’une image en niveaux de gris.
img_lena = imread('lena_std.tif');
% figure();
% imhist(img_lena);

% 8) Effectuer une égalisation d’histogramme sur l’image ’lena.png’ en utilisant la fonction histeq.
% img_lena_eq = histeq(img_lena);
% figure();
% imhist(img_lena_eq);
% figure();
% imshow(img_lena);
% figure();
% imshow(img_lena_eq);

% 9) Afficher l’image ’lena_noisy.png’ ainsi que son histogramme.
img_lena_noisy = imread('lena_noisy.png');
% figure();
% imshow(img_lena_noisy);
% figure();
% imhist(img_lena_noisy);

% 10) Calculer l’erreur quadratique moyenne entre l’image originale et l’image bruitée.
% img_lena_gray = rgb2gray(img_lena);
% MSE = reshape(mean(mean((double(img_lena_gray) - double(img_lena_noisy)).^2,2),1),[1,3]);
% display(MSE);
img_lena_noisy = imnoise(img_lena,'salt & pepper', 0.01);
Mean_Square_Error =immse(img_lena_noisy,img_lena);
% display(Mean_Square_Error);

% 11) Appliquer un filtre gaussien sur l’image bruitée.

% Malheureusement ça ne marche pas...
% TRIVIA(Anass): IEEE 745 + erreur de logic =
% digits(32);
% NBR1 = 0.01;
% NBR2 = 1;
% Mean_Square_Error1 =immse(img_lena_noisy,img_lena);
% while MSE1 > 1
%     img_lena_gauss = imgaussfilt(img_lena_noisy,NBR1);
%     Mean_Square_Error1 =immse(img_lena_noisy,img_lena);
%     img_lena_gauss = imgaussfilt(img_lena_noisy,NBR2);
%     Mean_Square_Error2 =immse(img_lena_noisy,img_lena);
%     if Mean_Square_Error1 <= Mean_Square_Error2
%         NBR2 = (NBR1 + NBR2) / 2;
%         display("A est proche: ");
%         display(num2str(NBR1,'%.16f'));
%         display("Mean_Square_Error1 is:");
%         display(num2str(Mean_Square_Error1,'%.10f'));
%     else
%         NBR1 = (NBR1 + NBR2) / 2;
%         display("B est proche:");
%         display(num2str(NBR2,'%.16f'));
%         display("Mean_Square_Error2 is:");
%         display(num2str(Mean_Square_Error2,'%.10f'));
%     end 
% end
% img_lena_gauss = imgaussfilt(img_lena_noisy,1);
% Mean_Square_Error1 =immse(img_lena_noisy,img_lena);
% % Note(Anass): MSE different chaque compilation????
% display(Mean_Square_Error1);
% figure();
% imshow(img_lena_gauss);

% 12) Même question avec le filtre moyenneur.
h = fspecial('average', [3 3]);
img_lena_moy = imfilter(img_lena_noisy, h);
% figure();
% imshow(img_lena_moy);
% Conclure.
% Le filtre moyenneur supprime le bruit tout en gardant les bords nets.

% 13) Appliquer le filtre médian sur l’image ’lena_abimee.png’ et observer le résultat avant et après restauration.
img_lena_noisy = imnoise(rgb2gray(img_lena),'salt & pepper',0.03);
img_lena_med = medfilt2(img_lena_noisy);
% figure();
% imshowpair(img_lena_noisy,img_lena_med,'montage')

% 14) Enregistrer l’image ’lena.png’ sous les formats suivants : JPEG, BMP, GIF et le format TIFF
% imwrite(img_lena,'img_lena.jpeg');
% 
% imwrite(img_lena,'img_lena.bmp');
% 
% filename = 'img_lena.gif';
% for n = 1:0.5:5
%     [imind,cm] = rgb2ind(img_lena,256);
%     if n == 1
%         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%         imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
% end
% 
% imwrite(img_lena,'img_lena.tiff');
% 
% % 15) Ouvrir et afficher chacune de ces images.
% figure();
% imshow('img_lena.jpeg');
% 
% figure();
% imshow('img_lena.bmp');
% 
% fullFileName = 'img_lena.gif';
% [gifImage, cmap] = imread(fullFileName, 'Frames', 'all');
% [rows, columns, numColorChannels, numImages] = size(gifImage);
% rgbImage = zeros(rows, columns, 3, numImages, 'uint8');
% hFig = figure;
% for k = 1 : numImages
%   Frame = gifImage(:,:,:, k);
%   RGB = uint8(255 * ind2rgb(Frame, cmap));
%   imshow(RGB);
%   rgbImage(:,:,:,k) = RGB;
%   drawnow;
% end
% 
% figure();
% imshow('img_lena.tiff');

% Observer leurs différences en taille et en qualité et comparer les avec l’image d’origine.
% img_lena_jpeg = imread('img_lena.jpeg');
% figure();
% imshowpair(img_lena,img_lena_jpeg,'montage')
% 
% img_lena_bmp = imread('img_lena.bmp');
% figure();
% imshowpair(img_lena,img_lena_bmp,'montage')
% 
% for k = 1 : numImages
%   Frame = gifImage(:,:,:, k);
%   RGB = uint8(255 * ind2rgb(Frame, cmap));
%   imshowpair(RGB,img_lena,'montage')
%   rgbImage(:,:,:,k) = RGB;
%   drawnow;
% end
% 
% img_lena_tiff = imread('img_lena.tiff');
% figure();
% imshowpair(img_lena,img_lena_tiff,'montage')

% Mesurer ces différences avec l’image d’origine en calculant l’erreur quadratique moyenne.
% Mean_Square_Error1 =immse(img_lena_jpeg,img_lena);
% Mean_Square_Error2 =immse(img_lena_bmp,img_lena);
% Mean_Square_Error3 =immse(RGB,img_lena);
% Mean_Square_Error4 =immse(img_lena_tiff,img_lena);
% display(Mean_Square_Error1);
% display(Mean_Square_Error2);
% display(Mean_Square_Error3);
% display(Mean_Square_Error4);

% 16) Enregistrer l’image ’lena.png’ au format JPEG avec différents niveaux de compression.
% imwrite(img_lena,'img_lena0.jpeg','Quality',0);
% imwrite(img_lena,'img_lena25.jpeg','Quality',25);
% imwrite(img_lena,'img_lena50.jpeg','Quality',50);
% imwrite(img_lena,'img_lena75.jpeg','Quality',75);
% imwrite(img_lena,'img_lena100.jpeg','Quality',100);

% Ouvrir et afficher ensuite chacune de ces images, observer leurs différences en taille et en qualité.
% img_lena0 = imread('img_lena0.jpeg');
% figure();
% imshowpair(img_lena,img_lena0,'montage')
% 
% img_lena25 = imread('img_lena25.jpeg');
% figure();
% imshowpair(img_lena,img_lena25,'montage')
% 
% img_lena50 = imread('img_lena50.jpeg');
% figure();
% imshowpair(img_lena,img_lena50,'montage')
% 
% img_lena75 = imread('img_lena75.jpeg');
% figure();
% imshowpair(img_lena,img_lena75,'montage')
% 
% img_lena100 = imread('img_lena100.jpeg');
% figure();
% imshowpair(img_lena,img_lena100,'montage')

% Comparer les avec l’image d’origine en calculant également l’erreur quadratique moyenne.
% Mean_Square_Error1 =immse(img_lena0,img_lena);
% Mean_Square_Error2 =immse(img_lena25,img_lena);
% Mean_Square_Error3 =immse(img_lena50,img_lena);
% Mean_Square_Error4 =immse(img_lena75,img_lena);
% Mean_Square_Error5 =immse(img_lena100,img_lena);
% display(Mean_Square_Error1);
% display(Mean_Square_Error2);
% display(Mean_Square_Error3);
% display(Mean_Square_Error4);
% display(Mean_Square_Error5);

% 17) Lire et afficher l’image ’jetons.bmp’, puis tester l’algorithme "flood fill" sur cette image.
% imwrite(imread('coins.png'),'coins.bmp');
% coins = imfill(imbinarize(imread('coins.bmp')),'holes');
% imshow(coins);

% 18) Transformer l’image ’jetons.bmp’ en image en niveaux de gris et en extraire son histogramme.
% coins = imread('coins.png');
% % Coins.png est déjà converti en niveaux de gris
% if size(coins, 3) > 1
%   I = rgb2gray(coins);
% end
coins = imread('coins.bmp');
% imhist(coins);

% 19) Binariser l’image en niveaux de gris de telle sorte à obtenir des objets en blanc et un fond en noir.
coins_bin = imbinarize(coins);
% figure();
% imshow(coins_bin);

% 20) Appliquer une érosion sur l’image binaire en utilisant comme élément structurant :
% SE = strel('square',3);
% coins_eros_square_3 = imerode(coins_bin,SE);
% SE = strel('square',3);
% coins_eros_square_5 = imerode(coins_bin,SE);
% SE = strel('sphere',3);
% coins_eros_sphere_3 = imerode(coins_bin,SE);
% SE = strel('sphere',5);
% coins_eros_sphere_5 = imerode(coins_bin,SE);
% multi = cat(3,coins_eros_square_3,coins_eros_square_5,coins_eros_sphere_3,coins_eros_sphere_5);
% montage(multi);

% 21) Appliquer une dilatation sur l’image binaire en utilisant les éléments structurants utilisés dans la question précédente.
% SE = strel('square',3);
% coins_eros_square_3 = imdilate(coins_bin,SE);
% SE = strel('square',3);
% coins_eros_square_5 = imdilate(coins_bin,SE);
% SE = strel('sphere',3);
% coins_eros_sphere_3 = imdilate(coins_bin,SE);
% SE = strel('sphere',5);
% coins_eros_sphere_5 = imdilate(coins_bin,SE);
% multi = cat(3,coins_eros_square_3,coins_eros_square_5,coins_eros_sphere_3,coins_eros_sphere_5);
% montage(multi);

% 22) Traiter l’image binarisée afin d’obtenir une image dans laquelle les formes correspondent au mieux aux objets de la scène réelle.
% coins_clear = bwareaopen(imclearborder(coins),1900);
% SE = strel('sphere',5);
% coins_eros_sphere_5 = imerode(coins_clear,SE);
% imshowpair(coins,coins_eros_sphere_5,'montage');

% 23) Détecter les contours des jetons présents dans l’image en niveaux de gris par une approche de type Sobel.
% coins_sobel = edge(coins,'Sobel');
% imshow(coins_sobel);

% 24) Affiner la segmentation en utilisant les opérateurs morphologiques.
% coins_sobel_fill = imfill(coins_sobel,'holes');
% imshow(coins_sobel_fill);

% 25) Utiliser les fonctions bwlabel, bwperim, bwarea.
coins_bwlabel = bwlabel(coins_bin);
coins_bwperim = bwperim(coins_bin);
coins_bwarea = bwarea(coins_bin);
% vislabels(coins_bwlabel);
% vislabels(coins_bwperim);
% % Fix(Anass) : vislabels(coins_bwarea)
% vislabels(coins_bwarea);
% str = {'coins_bwlabel','coins_bwperim','coins_bwarea'};
% text(int2str(coins_bwlabel),int2str(coins_bwperim),int2str(coins_bwarea),str);

% 26) Exécuter ce code et expliquer ce que chaque ligne effectue en vous aidant de l’aide Matlab.
% chemin = '.\';
% nb_images = 9;
% nb_classes = 3;
% num_classe = zeros(nb_classes,1);
% Att1 = zeros(nb_images,1);
% Att2 = zeros(nb_images,1);
% for i=1:nb_images
% 	fichier_image = [chemin int2str(i) '.bmp'];
% 	num_classe(i) = floor((i-1)/nb_classes) + 1;
% 	I = imread(fichier_image);
% 	I_ndg = rgb2gray(I);
% 	Att1(i) = mean(mean(I_ndg));
% 	I_ndg = double(I_ndg);
% 	Att2(i) = std(std(I_ndg));
% end
% plot(Att1(num_classe==1),Att2(num_classe==1),'c*',Att1(num_classe==2),Att2(num_classe==2),'b*',Att1(num_classe==3),Att2(num_classe==3),'r*');

% 27) Reprendre le code précédent en extrayant maintenant la corrélation et l’énergie de la matrice de co-occurrences.

chemin = '.\';
nb_images = 9;
nb_classes = 3;
num_classe = zeros(nb_classes,1);
Att1 = zeros(nb_images,1);
Att2 = zeros(nb_images,1);
for i=1:nb_images
	fichier_image = [chemin int2str(i) '.bmp'];
	num_classe(i) = floor((i-1)/nb_classes) + 1;
	I = imread(fichier_image);
    offsets0 = [zeros(40,1) (1:40)'];
    MCC = graycomatrix(rgb2gray(I), 'offset', [0 1], 'Symmetric', true);
    Haralick = graycoprops(MCC,'Energy');
    Att1(i) = [Haralick.Energy];
    Haralick = graycoprops(MCC,'Correlation');
    Att2(i) = [Haralick.Correlation];
end
plot(Att1(num_classe==1),Att2(num_classe==1),'c*',Att1(num_classe==2),Att2(num_classe==2),'b*',Att1(num_classe==3),Att2(num_classe==3),'r*');

%---------------------------------------------------------------------------------------------------------
function vislabels(L)
background_shade = 200;
foreground_shade = 240;
I = zeros(size(L), 'uint8');
I(L == 0) = background_shade;
I(L ~= 0) = foreground_shade;
imageHandle = imshow(I, 'InitialMagnification', 'fit');
axesHandle = ancestor(imageHandle, 'axes');
s = regionprops(L, 'Extrema');
hold(axesHandle, 'on');
for k = 1:numel(s)
   e = s(k).Extrema;
   text(e(1,1), e(1,2), sprintf('%d', k), ...
      'Parent', axesHandle, ...
      'Clipping', 'on', ...
      'Color', 'b', ...
      'FontWeight', 'bold');
end
hold(axesHandle, 'off');
end