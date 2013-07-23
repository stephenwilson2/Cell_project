function test_delta()
    datapts=3;
    imgs=cell(1,datapts);
    for i=1:datapts
        imgs{i}=MAKE(5);
        figure(1)
        subplot(1,3,1)
        imagesc(imgs{i}{1})
        subplot(1,3,2)
        imagesc(imgs{i}{1})
        subplot(1,3,3)
        imagesc(imgs{i}{1})
        pause
    end

end

function imgs=MAKE(numofmol)
    m=zeros(25,100);
    x=zeros(25,1);
    y=zeros(100,1);
    for i=1:numofmol
        x(i)=randi([1 25]);
        y(i)=randi([1 100]);
        if x(i)>=1 && x(i)<=25 && y(i)>=1 && y(i)<=100
            m(x(i),y(i))=m(x(i),y(i))+1;
        end
    end
    f=fspecial('gaussian',[100,100],1);
    f2=fspecial('gaussian',[100,100],10);
    img=imfilter(m,f);
    img2=imfilter(m,f2);
    imgs={m,img, img2};
end