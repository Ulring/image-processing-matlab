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
uint8Image2 = uint8(255 * rescale(IMG2));
%figure
%imshow(uint8Image1, 'InitialMagnification', 10000);
%figure
%imshow(uint8Image2, 'InitialMagnification', 10000);

%Exercice 2 :
exercice2 = [
    0 128 255 255 200 150 ;
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
exercice3 = [
    25 100 220 220 200 150 ;
    200 30 30 30 220 100 ;
    220 200 30 30 100 150 ;
    150 150 150 200 150 150 ;
    100 200 30 200 150 150 ;
    100 100 100 25 25 25]
IMG4 = mat2gray(exercice3);
%figure
%imshow(IMG4, 'InitialMagnification', 10000);
IMG5 = expansionDynamique(IMG4)
%figure
%imshow(IMG5, 'InitialMagnification', 10000);

%Exercice 4 :
exercice4 = [
    1 2 3 6 8 8 10 10 ;
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
IMG9 = imfilter(IMG3, test2);
%figure
%imshow(IMG9, 'InitialMagnification', 10000);
test3 = fspecial('average', [3 3]);
IMG10 = imfilter(IMG3, test3);
%figure
%imshow(IMG10, 'InitialMagnification', 10000);

%Exercice 6 :
exercice6_1 = [
    129 254 255 254 255 150;
    1 0 0 2 128 128;
    255 254 50 48 200 200;
    140 150 148 154 152 152;
    100 200 200 50 148 154;
    126 128 127 2 1 1]
exercice6_2 = [
    128 255 255 255 255 150;
    0 0 0 0 128 128;
    255 255 50 50 200 200;
    150 150 150 150 150 150;
    100 200 200 50 150 150;
    128 128 128 0 0 0]
I1 = mat2gray(exercice6_1);
I2 = mat2gray(exercice6_2);
imwrite(I1,'I1.png');
imwrite(I2,'I2.png');
im1 = dir('I1.png');
im2 = dir('I2.png');
ratio = im1.bytes/im2.bytes;
I1_rle = rle(I1);
I2_rle = rle(I2);
imwrite(I1_rle,'I1_rle.png');
imwrite(I2_rle,'I2_rle.png');
im1_rle = dir('I1_rle.png');
im2_rle = dir('I2_rle.png');
ratio_rle = im1_rle.bytes/im2_rle.bytes;
I1_huff = mat2gray(huffman(I1));
I2_huff = mat2gray(huffman(I2));
imwrite(I1_huff,'I1_huff.png');
imwrite(I2_huff,'I2_huff.png');
im1_huff = dir('I1_huff.png');
im2_huff = dir('I2_huff.png');
ratio_huff = im1_huff.bytes/im2_huff.bytes;
display(ratio);
display(ratio_rle);
display(ratio_huff);

%Exercice 7 :
exercice7 = [
    8 7 4 7 8 10 8 8;
    6 4 5 4 9 8 9 8 ;
    4 5 2 3 4 7 8 10 ;
    9 4 3 4 6 8 6 8 ;
    8 7 4 9 8 7 8 6 ;
    8 8 6 8 12 11 12 13 ;
    8 7 8 8 13 12 10 12 ;
    5 8 8 8 12 14 12 11 ];
serie = [0.1,0.2,0.3,0.4];
% for k = 1:numel(serie)
%     [g,NR,SI,TI] = regiongrow_test(mat2gray(exercice7),1,serie(k));
%     figure, imshow(TI, 'InitialMagnification', 10000);
% end
% Meilleur
% [g,NR,SI,TI] = regiongrow_test(mat2gray(exercice7),1,0.4);
% figure, imshow(TI, 'InitialMagnification', 10000);
% for k = 1:numel(serie)
%     J = regiongrowing_jeu_1(mat2gray(exercice7),1,1,serie(k));
%     figure, imshow(J, 'InitialMagnification', 10000);
% end
% Meilleur
% J = regiongrowing_jeu_1(mat2gray(exercice7),1,1,0.3);
% figure, imshow(J, 'InitialMagnification', 10000);
% for k = 1:numel(serie)
%     J = regiongrowing_jeu_2(mat2gray(exercice7),1,1,serie(k));
%     figure, imshow(J, 'InitialMagnification', 10000);
% end
% Meilleur
% J = regiongrowing_jeu_2(mat2gray(exercice7),1,1,0.3);
% figure, imshow(J, 'InitialMagnification', 10000);
% figure();
% imshow(mat2gray(exercice7), 'InitialMagnification', 10000);

%Exercice 8 :
histo = imhist(mat2gray(exercice7),256);
% figure();
% plot(histo);
% Note(Anass): Add It later

%Exercice 9 :
% Fix(Anass): splitmerge needs a remake 
% for k = 1:numel(serie)
%     g{k} = splitmerge(mat2gray(exercice7),serie(k),1);
%     figure, imshow(g{k})
% end

%Exercice 10 :
exercice10 = [
    8 7 9 8 8 7 8 8 ;
    6 6 7 8 9 8 9 8 ;
    8 8 9 7 7 7 8 6 ;
    1 0 3 1 2 2 0 0 ;
    0 2 1 1 2 3 0 0 ;
    2 2 0 0 1 1 0 0 ];
IMG10 = mat2gray(exercice10);
IMG10_SV = edge(IMG10,'Sobel','vertical');
IMG10_R = edge(IMG10,'Roberts');
IMG10_L = edge(IMG10,'log','vertical');
IMG10_med = medfilt2(IMG10);
IMG10_med_L = edge(IMG10_med,'log','vertical');
% figure();
% imshow(IMG10, 'InitialMagnification', 10000);
% figure();
% imshow(IMG10_SV, 'InitialMagnification', 10000);
% figure();
% imshow(IMG10_R, 'InitialMagnification', 10000);
% figure();
% imshow(IMG10_L, 'InitialMagnification', 10000);
% figure();
% imshow(IMG10_med_L, 'InitialMagnification', 10000);
% Le filtre Laplacien est marche mieux après le lissage.

%Exercice 11 :
R    = [0 0 50 100 200; 200 200 200 0 0; 0 100 150 100 200; 200 0 200 150 250; 200 200 200 200 0]
G    = [255 0 50 100 200; 255 200 0 50 250; 255 100 150 100 200; 0 128 255 200 250; 200 250 200 255 200]
B    = [0 128 50 100 200; 0 128 0 200 200; 128 100 150 100 200; 0 200 0 255 200; 128 200 0 200 200]
IMG11 = cat(3, R, G, B) / 255;
imshow(IMG11, 'InitialMagnification', 10000);
% Note(Anass): A continuer...

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
function [output] = rle(input)
L=length(input);
j=1;
k=1;
i=1;
while i<2*L
    comp=1;
    for j=j:L
        if j==L
            break
        end;
        if input(j)==input(j+1)
            comp=comp+1;
        else
            break
        end;
    end;
    output(k+1)=comp;
    output(k)=input(j);
    if j==L && input(j-1)==input(j)
        break
    end;
    i=i+1;
    k=k+2;
    j=j+1;
    if j==L
        if mod(L,2)==0
            output(k+1)=1;
            output(k)=input(j);
        else
            output(k+1)=1;
            output(k)=input(j);
        end;
        break
    end;
end;
end
function [output] = huffman(input)
A1 = input;
A = A1(:);
[symbols, ~, idx] = unique(A);
p = histcounts(idx, 1:max(idx)+1);
p = p/sum(p);
[dict, avglen] = huffmandict(symbols, p);
output = huffmanenco(A, dict)
end
function J=regiongrowing_jeu_1(I,x,y,reg_maxdist)
if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end
J = zeros(size(I)); 
Isizes = size(I);
reg_mean = I(x,y);
reg_size = 1;
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 
pixdist=0;
%  Ordre d’analyse du voisinage : droite, bas, gauche, haut.
neigb=[1 0; 0 1; -1 0;0 -1];
while(pixdist<reg_maxdist&&reg_size<numel(I))
    for j=1:4,
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    x = neg_list(index,1); y = neg_list(index,2);
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end
end
function J=regiongrowing_jeu_2(I,x,y,reg_maxdist)
if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end
J = zeros(size(I)); 
Isizes = size(I);
reg_mean = I(x,y);
reg_size = 1;
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 
pixdist=0;
% Ordre d’analyse du voisinage : bas, gauche, haut, droite.
neigb=[0 1; -1 0; 0 -1;1 0];
while(pixdist<reg_maxdist&&reg_size<numel(I))
    for j=1:4,
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    x = neg_list(index,1); y = neg_list(index,2);
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end
end
function [g,NR,SI,TI] = regiongrow_test(f,S,T)
f = im2double(f);
if numel(S) == 1
    SI = f == S;
    S1 = S;
else
    SI = bwmorph(S,'shrink',Inf);
    S1 = f(SI);
end
TI = false(size(f));
for K = 1:length(S1)
    seedvalue = S1(K);
    S = abs(f - seedvalue) <= T;
    TI = TI | S;
end
[g,NR] = bwlabel(imreconstruct(SI,TI));
end
function g = splitmerge(f, mindim, fun)
Q=2^nextpow2(max(size(f)));
[M, N]=size(f);
f=padarray(f,[Q-M, Q-N], 'post');
Z=qtdecomp(f,@split_test, mindim, fun);
Lmax=full(max(Z(:)));
g=zeros(size(f));
MARKER=zeros(size(f));
for K=1:Lmax
    [vals, r, c]=qtgetblk(f, Z, K);
    if ~isempty(vals)
        for I=1:length(r)
            xlow=r(I); ylow=c(I);
            xhigh=xlow+K-1;
            yhigh=ylow+K-1;
            region=f(xlow:xhigh, ylow:yhigh);
            flag=fun(region);
            if flag
                g(xlow:xhigh, ylow:yhigh)=1;
                MARKER(xlow, ylow)=1;
            end
        end
    end
end
g=bwlabel(imreconstruct(MARKER,g));
g=g(1:M, 1:N);
end
function v=split_test(B, mindim, fun)
k=size(B,3);
v(1:k)=false;
for I=1:k
    quadregion=B(:,:,I);
    if size(quadregion, 1)<= mindim
        v(I)=false;
        continue
    end
    flag=fun(quadregion);
    if flag
        v(I)=true;
    end
end
end
function flag = predicate(region)
sd = std2(region);
m = mean2(region);
flag = (sd > 0) & (m > 0) & (m < 125);
end
