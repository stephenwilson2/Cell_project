close all;
clear all;
c=onecell(10,250,1000,'sc',10,[250 1000],1,0);
c=c.cell_mask();
figure(15);
h1=subplot(5,1,1:4);
hold all;
imagesc(c);axis equal; plot(c); axis equal;
hold off;

% subplot(5,1,5)
figure(13)
% imshow(c);
mask=c.cellmask;
whos mask;
imagesc(c.cellmask')


