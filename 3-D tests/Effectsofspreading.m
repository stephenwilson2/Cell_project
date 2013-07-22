op=8;
numofmol=1:1:op;
o=0;
Vs=cell(op,1);
% hold all;
for num=numofmol
    o=o+1;
    m=zeros(25,100);
    x=zeros(25,1);
    y=zeros(100,1);
    img=cell(20,1);
    me=zeros(20,1);
    v=zeros(20,1);
    for i=1:num
        x(i)=12;
        y(i)=50-i*5;
        m(x(i),y(i))=m(x(i),y(i))+1;
    end
    n=0;
    for sig=.1:.1:50
        n=n+1;
        f=fspecial('gaussian',[100,100],sig);
        f2=fspecial('gaussian',[100,100],sig);
        img{n}=imfilter(m,f);
        figure(2);
        imagesc(m);
        me(n)=mean(mean(img{n}));
        v(n)=var(var((img{n})))*10000;
    end
    
    figure(3);
    hold all;
    subplot(1,3,1)
    plot(.1:.1:50, me/num)
    o=o+1;
    hold off;
    hold all;
    subplot(1,3,2)
    semilogx(.1:.1:50,v/num)
    Vs{round(o/2)}=v;
    hold off;
end
% hold off
subplot(1,3,1)
title(sprintf('The Normalized Mean Pixel Intensity against the Std. Dev. of the Gaussian \n Across All Numbers of Molecules'))
xlabel('\sigma (nm)')
ylabel('Mean Pixel Intensity')
legend('All numbers of molecules')
legend(sprintf('%d molecule',numofmol(1)),sprintf('%d molecules',numofmol(2)),sprintf('%d molecules',numofmol(3))...
    ,sprintf('%d molecules',numofmol(4)),sprintf('%d molecules',numofmol(5)),sprintf('%d molecules',numofmol(6)),...
    sprintf('%d molecules',numofmol(7)), sprintf('%d molecules',numofmol(8)))
title(sprintf('The Variance against the Std. Dev. of the Gaussian \n and Varying amounts of Molecules'))
xlabel('\sigma (nm)')

subplot(1,3,2)
legend(sprintf('%d molecule',numofmol(1)),sprintf('%d molecules',numofmol(2)),sprintf('%d molecules',numofmol(3))...
    ,sprintf('%d molecules',numofmol(4)),sprintf('%d molecules',numofmol(5)),sprintf('%d molecules',numofmol(6)),...
    sprintf('%d molecules',numofmol(7)), sprintf('%d molecules',numofmol(8)))
title(sprintf('The Variance against the Std. Dev. of the Gaussian \n and Varying amounts of Molecules'))
xlabel('\sigma (nm)')
ylabel('Variance (1/nm^2)')
subplot(1,3,3)
hold all;
xs=zeros(8);
ys=zeros(8);
s=cell(5,1);
in=0;
for o=5:1:10
    in=in+1;
    for i=1:op
        ys(i)=Vs{i}(o);
        xs(i)=i;
    end
    semilogy(xs,ys)
    s{in}=sprintf('%g nm of spread',o/10);
end
hold off;
title(sprintf('The Number of Molecules versus the Variance'))
xlabel('Number of Molecules')
ylabel('Variance (1/nm^2)')
legend(s)