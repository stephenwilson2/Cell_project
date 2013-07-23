function testIntensity3D()        
    datapts=5;
    molpcell=round(1:1:400);
    if ~isequal(exist(sprintf('testIntensity3D_%i_%i.mat',datapts,max(molpcell)),'file'),2)
        r=250;
        w=1000;
        c=cell(datapts,1);
        i=0;
        for numofmol=molpcell
            for o=1:datapts
                i=i+1
                c{i}=onecell(numofmol,r,w,'sc',64);
            end
        end
        save(sprintf('testIntensity3D_%i_%i',datapts,max(molpcell)))
    else
        load(sprintf('testIntensity3D_%i_%i',datapts,max(molpcell)))
    end
    analyze(c,datapts,molpcell);
end

function analyze(cells,pts,molpcell)
    e=0;
    pairpsf(length(cells)/pts+1,3)=0;
    for i=1:length(cells)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            e=e+1;
            V(o)=var(cells{e}.img{1}(:)); %Should this be 
            n(o)=mean(cells{e}.img{1}(:));
        end
        pairpsf(i,1)=mean(n);
        pairpsf(i,2)=mean(V);
        pairpsf(i,3)=std(V)/pts^0.5;
    end
    %with PSF
    theok=9.5567e-004;% photons per voxel
    figure(74);
    [p,s]=polyfit(pairpsf(:,1), pairpsf(:,2),1);
    y=pairpsf(:,2);
    yfit = polyval(p,pairpsf(:,1));
    R2psf = corrcoef(pairpsf(:,2), yfit);
    theox=0:.01:max(pairpsf(:,1));
    theoy=theok*theox;
    hold all;
    errorbar(pairpsf(:,1),pairpsf(:,2),pairpsf(:,3),'ob');
    plot(pairpsf(:,1),yfit,'color', 'red');
    plot(theox,theoy,'color','blue')
    hold off;
    title('Variance compared to Mean Pixel Intensity',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')

    n1=sprintf('Fit-Slope: %d photons/pixel \n Intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    n2=sprintf('Theoretical Slope: %d photons/pixel',theok);
    legend('Simulation with PSF',n1,n2)
    
    saveas(gcf, sprintf('testIntensity3D_%i_%i.fig',pts,max(molpcell)))
    
    save(sprintf('testIntensity3D_%i_%i',pts,max(molpcell)),'-append','pairpsf')
end
