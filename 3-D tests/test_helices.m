function test_helices()
    numrange=25:25:125;
    l=1000;
    r=250;
    pxsz=64;
    shape='sc';
    %normal cell
    if ~isequal(exist(sprintf('test_helices%s_%i_%i_%i.mat',shape,pxsz,r,l),'file'),2)
        graph(shape,numrange,pxsz,r,l);
    end
    
    %no cap cell
    shape='b';
    l=2000;
    if ~isequal(exist(sprintf('test_helices%s_%i_%i_%i.mat',shape,pxsz,r,l),'file'),2)
        graph(shape,numrange,pxsz,r,l)
    end
    l=1000;
    shape='sc';
    
    %High res
    pxsz=5;
    if ~isequal(exist(sprintf('test_helices%s_%i_%i_%i.mat',shape,pxsz,r,l),'file'),2)
        graph(shape,numrange,pxsz,r,l)
    end
    pxsz=64;
    
    %Fission yeast
    r=1750;
    l=5000;
    if ~isequal(exist(sprintf('test_helices%s_%i_%i_%i.mat',shape,pxsz,r,l),'file'),2)
        graph(shape,numrange,pxsz,r,l) % mod pxsz or numrange (the density is diff!)????????????????????????????????
    end
    r=250;
    l=1000;
    
    %long
    l=10000;
    if ~isequal(exist(sprintf('test_helices%s_%i_%i_%i.mat',shape,pxsz,r,l),'file'),2)
        graph(shape,numrange,pxsz,r,l)
    end
end

function graph(shape,numrange,pxsz,r,l)
    figure(1);
    pts=3;
    n=0;
    for i=1:length(numrange)
        for o=1:pts
            n=n+1;
            subplot(length(numrange)*pts/3,3,n);
            c=onecell(numrange(i),r,l,shape,pxsz);
            c.showslice(4);
            c=c.refresh_all();
            c.showslice(4);
            c=c.refresh_all();
            c.showslice(4);
        end
    end
    save(sprintf('test_helices%s_%i_%i_%i',shape,pxsz,r,l))
    saveas(gcf,sprintf('test_helices%s_%i_%i_%i.fig',shape,pxsz,r,l))
%     saveas(gcf,sprintf('test_helices%s_%i_%i_%i.jpg',shape,pxsz,r,l))
end