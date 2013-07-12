function test_helices()
    clear all;
    close all;
%     if ~isequal(exist('test_helices.mat','file'),2)
        datapts=3;
        molpcell=round(25:25:125);
        l=length(molpcell);
        r=250;
        w=1000;
        c=cell(datapts,1);
        i=0;
        for numofmol=molpcell
            for o=1:datapts
                i=i+1;
                c{i}=onecell(numofmol,r,w,'sc',64,[r w],1,45);
            end
        end
        save('test_helices_rotated')
%     else
%         load('test_helices')
%     end
    analyze(c,datapts,l);
end

function analyze(cel,pts,l)
figure(1);
    for i=1:length(cel)
            i
            subplot(l,pts,i)
            imagesc(cel{i}'); axis equal;
            axis tight; axis off;
    end
    saveas(gcf,'test_helices_rotated.fig')
end