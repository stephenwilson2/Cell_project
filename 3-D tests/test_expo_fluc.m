function test_expo_fluc()
    if ~isequal(exist('test_expo.mat','file'),2)
        datapts=20;
        imgnum=100;
        numofmol=100:100:100*imgnum;
        imgs=cell(imgnum*datapts,1);
        i=0;
    %     figure(2);
        for u=1:datapts        
            for n=1:imgnum
                i=i+1;
                imgs{i}=MAKE(numofmol(n));
    %             subplot(4,imgnum*datapts,i)
    %             imagesc(imgs{i}{1})
    %             subplot(4,imgnum*datapts,i+imgnum*datapts)
    %             imagesc(imgs{i}{2})
            end
        end
        save('test_expo')
    else
        load('test_expo')
    end
    analyze(imgs,datapts,numofmol);
end

function analyze(imgs,pts,molpcell)
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
    %with small PSF    
    figure(74);
    [p,s]=polyfit(pairspsf(:,1), pairspsf(:,2),1);
    y=pairspsf(:,2);
    yfit = polyval(p,pairspsf(:,1));
    R2psf = corrcoef(pairspsf(:,2), yfit);    
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
    
%     saveas(gcf, sprintf('testexpo3D_%i_%i.fig',pts,max(molpcell)))
    
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
    hold all;
    errorbar(pairlpsf(:,1),pairlpsf(:,2),pairlpsf(:,3),'ob');
    plot(pairlpsf(:,1),yfit,'color', 'red');
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    legend('Simulation with PSF',n1)
    
%     saveas(gcf, sprintf('testexpo3D_%i_%i.fig',pts,max(molpcell)))
    
end



function imgs=MAKE(numofmol)
    m=zeros(25,100);
    x=zeros(75,1);
    y=zeros(300,1);
    for i=1:numofmol
        x(i)=randi([-25 50]);
        y(i)=randi([-100 200]);
        if x(i)>=1 && x(i)<=25 && y(i)>=1 && y(i)<=100
            m(x(i),y(i))=m(x(i),y(i))+1;            

        end
    end
    f=fspecial('gaussian',[100,100],1);
    f2=fspecial('gaussian',[100,100],10);
    img=imfilter(m,f);
    img2=imfilter(m,f2);
    imgs={img, img2};
end