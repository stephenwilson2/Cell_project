close all;
clear all;
%% Basic, almost function-like, way of manipulating the onecell
c=onecell(1,250,1000,'sc',64,[250 1000,250],1);
imagesc(c);
c

%% An Object-oriented implementation
close all;
clear all;
b=onecell();
b
b.numofmol=400;
b.r=300;
b.pixelsize=10;
b
b=b.refresh_cell();
imshow(b);
close all;
b.pixelsize=64;
imshow(b);