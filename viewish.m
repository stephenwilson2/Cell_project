function viewish()
    c=onecell(64);
    threed=zeros(size(c.img{1},1),size(c.img{1},2),length(c.img));
    for s=1:length(c.img)
        threed(:,:,s)=c.img{s};
    end
%     size(threed)
%     figure(2)
%     [x,y,z] = meshgrid(1:1:size(threed,2),1:1:size(threed,1),1:1:size(threed,3));
%     size(x)
%     slice(x,y,z,threed,[],[],[1,2,3,4,5,6,7,8])
    colormap(gray);
    colorbar;
    [x,y] = meshgrid(1:1:size(c.img{1},2),1:1:size(c.img{1},1));
    z=c.img{1};
    surf(x,y,z)
    view([118,64]);axis tight
end