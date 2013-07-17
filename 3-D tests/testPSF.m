function testPSF()
%TESTPSF Tests the PSF to ensure that it is working
    clear all;
    close all;
    c=onecell(1,250,1000,'sc',1,[250 1000 250],1,500);
    s=0;
    for i=1:length(c.img)
        s=s+sum(sum(c.img{i}));
    end
    s

end

