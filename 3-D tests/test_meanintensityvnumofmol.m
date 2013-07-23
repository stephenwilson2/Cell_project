function test_meanintensityvnumofmol(ROI)
    %works with spherocylinders only right now
%     ROI=0;
    if ~isequal(exist(sprintf('mvn_ROI_%i.mat',ROI),'file'),2)        
        datapts=10;
        numofmol=round(linspace(1,1000,10)); %min, max, number of pts b/w
        c=cell(length(numofmol),1);
        m=zeros(length(numofmol),1);
        v=zeros(length(numofmol),1);
        sem=zeros(length(numofmol),1);
        sev=zeros(length(numofmol),1);
        tmpm=zeros(datapts,1);
        tmpv=zeros(datapts,1);
        
        n=0;
        for i=1:length(numofmol)
            numofmol(i)
            for o=1:datapts
                n=n+1;
                c{n}=onecell(numofmol(i));
                if ROI
                    f=round(c{n}.PSF.sigma/c{n}.pixelsize/2);  %half a PSF on each side
                    img=c{n}.img{1};
                    tmpm(o)=mean(mean(img(f:size(img,1)-f,f:size(img,2)-f))); 
                    tmpv(o)=mean(var(img(f:size(img,1)-f,f:size(img,2)-f)));
                else
                    tmpm(o)=mean(mean(c{n}.img{1}(:)));
                    tmpv(o)=var(c{n}.img{1}(:));
                end
            end
            m(i)=mean(tmpm);
            v(i)=mean(tmpv);
            sem(i)=std(tmpm)/datapts^0.5;
            sev(i)=std(tmpv)/datapts^0.5;
        end
        save(sprintf('mvn_ROI_%i.mat',ROI))
    else
        load(sprintf('mvn_ROI_%i.mat',ROI))
    end
    if ROI
        s='The ROI is 1 PSF on both axes smaller';
    else
        s='THE ROI is the entire cell';
    end
        
    figure('name',s);
    subplot(1,2,1)
    errorbar(numofmol,m,sem);
       title('Mean Compared to Number of Molecules',...
            'FontWeight','bold')
    ylabel('Mean Pixel Intensity (1/pixel)')
    xlabel('Number of Molecules')

    if ROI
%         l=c{n}.l-f-f*c{n}.pixelsize;
%         r=c{n}.r-f-f*c{n}.pixelsize/2;
%         a=2*pi*r+(l-2*r)*r;
        theok=3.6456e-009*c{n}.pixelsize^2;
    else
%         l=c{n}.l;
%         r=c{n}.r;
%         a=2*pi*r+(l-2*r)*r;
        theok=3.6456e-009*c{n}.pixelsize^2;% photons per area of z-slice. The area is the area of the cell not the camera
    end
    subplot(1,2,2)
    hold all;
    errorbar(m,v,sev)
    
    theox=0:.001:max(m);
    theoy=theok*theox;
    plot(theox,theoy,'color','red');
    hold off;
    title('Variance Compared to Mean Pixel Intensity',...
            'FontWeight','bold')
    xlabel('Mean Pixel Intensity (1/pixel)')
    ylabel('Variance (1/pixel^2)')
    legend('Simulation','Theory')
    
end