% Manipulation 1 avec Matlab :
IMG1=imread('lena_std.tif');
% 1)
figure();
imshow(IMG1);
% % 2)
exemple_verticale(IMG1);
exemple_horizontale(IMG1);


% Manipulation 2 avec Matlab :
% 3)
I = nanmean(IMG1,3);
figure(); imshow(I,[]); daspect([1 1 1 ]); axis off;
colormap jet;
% % 4)
% % Couleur pris de mon colorscheme Neovim
hex = ['#6d00c1';'#8510d8';'#550a8a';'#8601af';'#e9d8f2';'#a167c9';'#363636']
map = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255
I = nanmean(IMG1,3);
figure(); imshow(I,[]); daspect([1 1 1 ]); axis off;
colormap(map);

% Manipulation 3 avec Matlab :
% 5)
I = nuancesDeGris1(IMG1);
figure();
imshow(I);
% % 6) Y-a-t-il une perte d'informations lorsqu’on passe de l'espace RGB aux
% % niveaux de gris? : Oui

% Manipulation 4 avec Matlab :
% 7)
I = rgb2hsv(IMG1)
figure();
imshow(I);

% Manipulation 5 avec Matlab : 
% 8)
I = imread('onion.png');
Linear_Bit(nuancesDeGris1(I))

% Manipulations avancées
% 9)
I = imread('cell.tif');
manipulation1(I);
manipulation2(IMG1);
% 10)
I = imread('cameraman.tif');
imwrite(I,'cameraman.jpeg');
imwrite(I,'cameraman.png');
Ijpeg = imread('cameraman.jpeg');
Ipng = imread('cameraman.png');
K = imabsdiff(Ijpeg,Ipng);
imagesc(K)
%--------------------------------------------------------------------------
function exemple_verticale(z)
    for x=1 : size(z,1)
        for y=1 : size(z,2)/2
        end
    end
    z = imcrop(z,[0 0 y x]);
    figure();
    imshow(z);
end
function exemple_horizontale(z)
    for x=1 : size(z,1)/2
        for y=1 : size(z,2)
        end
    end
    z = imcrop(z,[0 0 y x]);
    figure();
    imshow(z);
end
function [gray_img] = nuancesDeGris1(img)
    R=img(:, :, 1);
    G=img(:, :, 2);
    B=img(:, :, 3);
    [M, N, ~]=size(img);
    gray_img=zeros(M, N, 'uint8');
    for x=1:M
        for y=1:N
            gray_img(x, y)=(R(x, y)*0.2989)+(G(x, y)*0.5870)+(B(x, y)*0.114);
        end
    end
end
function y = Linear_Bit(Img)
    b1 = double(bitget(Img,1));
    b2 = double(bitget(Img,2));
    b3 = double(bitget(Img,3));
    b4 = double(bitget(Img,4));
    b5 = double(bitget(Img,5));
    b6 = double(bitget(Img,6));
    b7 = double(bitget(Img,7));
    b8 = double(bitget(Img,8));
    Img_b = cat(8,b1,b2,b3,b4,b5,b6,b7,b8);
    figure, 
    subplot(4,2,1)
    imshow(b1), title('Bit Plan: 1');
    subplot(4,2,2)
    imshow(b2), title('Bit Plan: 2');
    subplot(4,2,3)
    imshow(b3), title('Bit Plan: 3');
    subplot(4,2,4)
    imshow(b4), title('Bit Plan: 4');
    subplot(4,2,5)
    imshow(b5), title('Bit Plan: 5');
    subplot(4,2,6)
    imshow(b6), title('Bit Plan: 6');
    subplot(4,2,7)
    imshow(b7), title('Bit Plan: 7');
    subplot(4,2,8)
    imshow(b8), title('Bit Plan: 8');
    y = Img_b;
end
function manipulation1(z)
    %carre pour encadrer pixel modifié
    z(100-1:100+1, 20 - 1) = 255; % Left edge
    z(100-1:100+1, 20 + 1) = 255; % Right edge
    z(100-1, 20) = 255; % Top edge
    z(100+1, 20) = 255; % Bottom edge
    for x=1 : size(z,1)/2
        for y=1 : size(z,2)
            if (x == 100) && (y == 20)
                z(x,y)=z(x,y)+25;
            end
        end
    end
    figure();
    imshow(z, 'InitialMagnification', 800);
    for x=1 : size(z,1)/2
        for y=1 : size(z,2)
            if (x == 100) && (y == 20)
                z(x,y)=z(x,y)-25;
            end
        end
    end
    figure();
    imshow(z, 'InitialMagnification', 800);
end
function manipulation2(z)
    %carre pour encadrer pixel modifié
    z(100-1:100+1, 20 - 1) = 255; % Left edge
    z(100-1:100+1, 20 + 1) = 255; % Right edge
    z(100-1, 20) = 255; % Top edge
    z(100+1, 20) = 255; % Bottom edge
    for x=1 : size(z,1)/2
        for y=1 : size(z,2)
            if (x == 100) && (y == 20)
                z(x,y)=z(x,y)+25;
            end
        end
    end
    figure();
    imshow(z, 'InitialMagnification', 200);
    for x=1 : size(z,1)/2
        for y=1 : size(z,2)
            if (x == 100) && (y == 20)
                z(x,y)=z(x,y)-25;
            end
        end
    end
    figure();
    imshow(z, 'InitialMagnification', 200);
end