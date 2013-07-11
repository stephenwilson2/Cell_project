function test_helices()
    if ~isequal(exist('test_helices.mat','file'),2)
        datapts=3;
        molpcell=round(10:10:100);
        l=length(molpcell);
        r=250;
        w=1000;
        c=cell(datapts,1);
        i=0;
        for numofmol=molpcell
            for o=1:datapts
                i=i+1;
                c{i}=onecell(numofmol,r,w,'sc',64,[250 1000],1);
            end
        end
%         save('test_helices')
    else
        load('test_helices')
    end
    analyze(c,datapts,l);
end

function analyze(cel,pts,l)
figure(1);
    for i=1:length(cel)
            i
            subplot(l,pts,i)
            imagesc(cel{i})
            axis tight; axis off;
            
            s=num2str(cel{i}.numofmol);
            s=strcat('\color{white}',s,' Molecules');
            text(.9,1.1,s);
    end
end