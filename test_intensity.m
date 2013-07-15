function test_intensity()
    clear all;
    close all;
    if ~isequal(exist('test_intensity.mat','file'),2)
        datapts=20;
        molpcell=round(1:1:200);
        r=250;
        w=1000;
        c=cell(datapts,1);
        i=0;
        for numofmol=molpcell
            for o=1:datapts
                i=i+1
                c{i}=onecell(numofmol,r,w,'sc',64,[250 1000],1);
            end
        end
        save('test_intensity')
    else
        load('test_intensity')
    end
    analyze(c,datapts);
end

function analyze(cells,pts)
    pairpsf(length(cells)/pts+1,3)=0;
    for i=1:length(cells)/pts
        V=zeros(pts,1);
        n=zeros(pts,1);
        for o=1:pts
            V(o)=var(cells{i+o}.img(:));
            n(o)=mean(cells{i+o}.img(:));
        end
        pairpsf(i,1)=mean(n);
        pairpsf(i,2)=mean(V);
        pairpsf(i,3)=std(V)/pts^0.5;
    end
    %with PSF    
    figure(74);
    [p,s]=polyfit(pairpsf(:,1), pairpsf(:,2),1);
    y=pairpsf(:,2);
    yfit = polyval(p,pairpsf(:,1));
    R2psf = corrcoef(pairpsf(:,2), yfit);    
    hold all;
    errorbar(pairpsf(:,1),pairpsf(:,2),pairpsf(:,3),'ob');
    plot(pairpsf(:,1),yfit,'color', 'red');
    hold off;
    title('Variance compared to number of molecules',...
        'FontWeight','bold')
    xlabel('Mean Pixel Intensity')
    ylabel('Variance')

    n1=sprintf('Fit-Slope: %d, intercept %d R^2: %d', p(1),p(2), R2psf(1,2));
    legend('Simulation with PSF',n1)
    
    saveas(gcf, 'TestIntensity.fig')
    
    save('test_intensity','-append','pairpsf')
end