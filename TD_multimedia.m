%Exercice 1 :
R    = [0 128 255;255 0 0;255 255 0]
G    = [0 128 255;0 255 0;255 0 255]
B    = [0 128 255;0 0 255;0 255 255]
IMG1 = cat(3, R, G, B) / 255;
%figure
%imshow(IMG1, 'InitialMagnification', 10000);
%figure
%imshow(rgb2gray(IMG1), 'InitialMagnification', 10000);
matrice_de_passage = [1/3 1/3 1/3;1/2 0 1/2;-1/4 1/2 -1/4]
IMG2 = cat(3, R * matrice_de_passage, G * matrice_de_passage, B * matrice_de_passage) / 255;
%figure
%imshow(IMG2, 'InitialMagnification', 10000);
uint8Image1 = uint8(255 * mat2gray(IMG2));
uint8Image2 = uint8(255 *  rescale(IMG2));
%figure
%imshow(uint8Image1, 'InitialMagnification', 10000);
%figure
%imshow(uint8Image2, 'InitialMagnification', 10000);

%Exercice 2 :
exercice2 = [0 128 255 255 200 150 ;
    255 0 50 0 255 128 ;
    255 255 50 0 128 150 ;
    150 150 150 200 150 150 ;
    128 200 50 200 150 150 ;
    100 128 128 0 0 0]
IMG3 = mat2gray(exercice2);
%figure
%imshow(IMG3, 'InitialMagnification', 10000);
%imhist(IMG3)

%Exercice 3 :
exercice3 = [25 100 220 220 200 150 ;
    200 30 30 30 220 100 ;
    220 200 30 30 100 150 ;
    150 150 150 200 150 150 ;
    100 200 30 200 150 150 ;
    100 100 100 25 25 25]
IMG4 = mat2gray(exercice3);
%figure
%imshow(IMG4, 'InitialMagnification', 10000);
%IMG5 = expansionDynamique(IMG4)
%figure
%imshow(IMG5, 'InitialMagnification', 10000);

%Exercice 4 :
exercice4 = [1 2 3 6 8 8 10 10 ;
    2 4 5 7 8 11 11 10 ;
    3 5 7 9 12 13 11 8 ;
    6 7 9 14 15 12 9 7 ;
    8 8 12 13 14 9 7 6 ;
    8 11 13 12 9 5 6 5 ;
    10 10 11 9 7 6 4 4 ;
    9 10 8 7 6 5 4 3]
IMG6 = mat2gray(exercice3);
%figure
%imshow(IMG6, 'InitialMagnification', 10000);
IMG7 = histeq(IMG6)
%figure
%imshow(IMG7, 'InitialMagnification', 10000);

%Exercice 5 :
test1 = fspecial('average', [5 5]);
IMG8 = imfilter(IMG3, test1);
%figure
%imshow(IMG8, 'InitialMagnification', 10000);
test2 = fspecial('average', [1 1]);
IMG8 = imfilter(IMG3, test2);
%figure
%imshow(IMG8, 'InitialMagnification', 10000);
test3 = fspecial('average', [3 3]);
IMG8 = imfilter(IMG3, test3);
%figure
%imshow(IMG8, 'InitialMagnification', 10000);

%----------------------------------------------------------------------------------

function [output] = doubleToUint8(input)
% cette fonction prend en entrée une image et renvoie la même image en format uint8
output = uint8(floor(255*input));
end
function [output] = uint8ToDouble(input)
    output = double(input) / 255;
end
function [output] = expansionDynamique(input)
    image = uint8ToDouble(input);
    M = max(max(image));
    m = min(min(image));
    dim = size(image);
    temp = zeros(dim(1),dim(2));
    for j = 1:dim(2)
        for i = 1:dim(1)
            temp(i,j) = (image(i,j)-m)/(M-m);
        end
    end
    output = doubleToUint8(temp);
end