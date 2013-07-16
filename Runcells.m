close all;
clear all;
c=onecell(1,250,1000,'sc',64,[250 1000],1,0);

figure(15);
hold all;
imagesc(c);axis equal; plot(c); axis equal;
hold off;
c