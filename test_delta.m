function test_delta()
    close all;
    datapts=30;
    numofmol=linspace(1,1000,10);
    imgs=cell(1,datapts);
    n=0;
    for o=1:length(numofmol)
        for i=1:datapts        
            n=n+1;
            imgs{n}=MAKE(numofmol(o));
            figure(1)
            subplot(1,4,1)
            imagesc(imgs{n}{1}')
            title('Images of cells. Resolution 100 nm^2/pixel')
            subplot(1,4,2)
            imagesc(imgs{n}{2}')
            title('Images of cells. Resolution 100 nm^2/pixel')
            subplot(1,4,3)
            imagesc(imgs{n}{3}')
            title('Images of cells. Resolution 100 nm^2/pixel')
            subplot(1,4,4)
            imagesc(imgs{n}{4}')
            title('Images of cells. Resolution 100 nm^2/pixel')
        end
    end
    analyze(imgs,datapts);
end
function k=calctheo(sig)
    k=1/(4*pi*sig*sig)-1/(25*100)
end
function analyze(imgs,pts)
    pairspsf(length(imgs)/pts+1,3)=0;
    e=0;
    for i=1:length(imgs)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            e=e+1;
            V(o)=var(imgs{e}{1}(:));
            n(o)=mean(imgs{e}{1}(:));
        end
        pairspsf(i,1)=mean(n);
        pairspsf(i,2)=mean(V);
        pairspsf(i,3)=std(V)/pts^0.5;
    end
    %with no PSF    
    figure(74);
    [p,s]=polyfit(pairspsf(:,1), pairspsf(:,2),1);
    y=pairspsf(:,2);
    yfit = polyval(p,pairspsf(:,1));
    R2psf = corrcoef(pairspsf(:,2), yfit);
    theox=0:.001:max(sort(pairspsf(:,1)));
    hold all;
    errorbar(pairspsf(:,1),pairspsf(:,2),pairspsf(:,3),'ob');
    plot(pairspsf(:,1),yfit,'color', 'red');
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    legend('Simulation with PSF',n1)
    
    e=0;
    pairlpsf(length(imgs)/pts+1,3)=0;
    for i=1:length(imgs)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            e=e+1;
            V(o)=var(imgs{e}{2}(:));
            n(o)=mean(imgs{e}{2}(:));
        end
        pairlpsf(i,1)=mean(n);
        pairlpsf(i,2)=mean(V);
        pairlpsf(i,3)=std(V)/pts^0.5;
    end
    %with small PSF    
    figure(76);
    [p,s]=polyfit(pairlpsf(:,1), pairlpsf(:,2),1);
    y=pairlpsf(:,2);
    yfit = polyval(p,pairlpsf(:,1));
    R2psf = corrcoef(pairlpsf(:,2), yfit);
    theox=0:.001:max(sort(pairlpsf(:,1)));
    theoy=calctheo(1)*theox;
    hold all;
    errorbar(pairlpsf(:,1),pairlpsf(:,2),pairlpsf(:,3),'ob');
    plot(pairlpsf(:,1),yfit,'color', 'red');
    plot(theox,theoy,'color','blue')
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    n2=sprintf('Theoretical: sigma=%i nm',10);
    legend('Simulation with PSF',n1,n2)
     
    e=0;
    pairlpsf(length(imgs)/pts+1,3)=0;
    for i=1:length(imgs)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            e=e+1;
            V(o)=var(imgs{e}{4}(:));
            n(o)=mean(imgs{e}{4}(:));
        end
        pairlpsf(i,1)=mean(n);
        pairlpsf(i,2)=mean(V);
        pairlpsf(i,3)=std(V)/pts^0.5;
    end
   %with large PSF    
    figure(78);
    [p,s]=polyfit(pairlpsf(:,1), pairlpsf(:,2),1);
    y=pairlpsf(:,2);
    yfit = polyval(p,pairlpsf(:,1));
    R2psf = corrcoef(pairlpsf(:,2), yfit);
    theox=0:.001:max(sort(pairlpsf(:,1)));
    theoy=calctheo(10)*theox;
    hold all;
    errorbar(pairlpsf(:,1),pairlpsf(:,2),pairlpsf(:,3),'ob');
    plot(pairlpsf(:,1),yfit,'color', 'red');
    plot(theox,theoy,'color','blue')
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    n2=sprintf('Theoretical: sigma=%i nm',100);
    legend('Simulation with PSF',n1,n2)
    
    e=0;
        
    pairlpsf(length(imgs)/pts+1,3)=0;
    for i=1:length(imgs)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            e=e+1;
            V(o)=var(imgs{e}{3}(:));
            n(o)=mean(imgs{e}{3}(:));
        end
        pairlpsf(i,1)=mean(n);
        pairlpsf(i,2)=mean(V);
        pairlpsf(i,3)=std(V)/pts^0.5;
    end
    %with real PSF
    figure(80);
    [p,s]=polyfit(pairlpsf(:,1), pairlpsf(:,2),1);
    y=pairlpsf(:,2);
    yfit = polyval(p,pairlpsf(:,1));
    R2psf = corrcoef(pairlpsf(:,2), yfit);
    theox=0:.001:max(sort(pairlpsf(:,1)));
    theoy=calctheo(8.0)*theox;
    hold all;
    errorbar(pairlpsf(:,1),pairlpsf(:,2),pairlpsf(:,3),'ob');
    plot(pairlpsf(:,1),yfit,'color', 'red');
    plot(theox,theoy,'color','blue')
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    n2=sprintf('Theoretical: sigma=%i nm',80);
    legend('Simulation with PSF',n1,n2)
    
end


function imgs=MAKE(numofmol)
    m=zeros(25,100);
    x=zeros(25,1);
    y=zeros(100,1);
    for i=1:numofmol
        x(i)=randi([1 25]);
        y(i)=randi([1 100]);
        if x(i)>=1 && x(i)<=25 && y(i)>=1 && y(i)<=100
            m(x(i),y(i))=m(x(i),y(i))+1;
        end
    end
    f=fspecial('gaussian',[100,100],1);
    f2=fspecial('gaussian',[100,100],10);
    f3=fspecial('gaussian',[100,100],8.4);
    img=imfilter(m,f);
    img2=imfilter(m,f2);
    img3=imfilter(m,f3);
    imgs={m,img, img3, img2};
end