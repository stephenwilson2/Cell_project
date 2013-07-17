function testPSF()
%TESTPSF Tests the PSF to ensure that it is working
    clear all;
    close all;
    c=onecell(1,250,1000,'sc',10,[250 1000 250],1,50);
    s=0;
    syms a;
    d=double(int(c.PSF.gaussDistribution(a,0,c.PSF.sigmaz), -inf,inf));
    disp(sprintf('Integration of the zPSF gives a value of %d. \n',d));

    for i=1:length(c.img)
        s=s+sum(sum(c.img{i}));
    end
    
    atz=sum(sum(c.img{round(c.fl(1,3)/c.pixelsize)}));
    disp(sprintf('A z-slice of the simulation has a value of %d at the origin of the molecule.\n',atz));
    atz2=c.PSF.gaussDistribution(1,1,c.PSF.sigmaz/c.PSF.pxsz);
    disp(sprintf('That z-slice should have a value of %d in order to have the whole cell sum to 1.\n',atz2));
    disp(sprintf('The apparent sum of the cell is %d. \n',s));
    er=(1-s)*100;
    disp(sprintf('The percent error from reality is %f%%. \n',er));
end

