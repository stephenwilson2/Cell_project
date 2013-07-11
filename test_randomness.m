function test_randomness()
    if ~isequal(exist('test_randomness.mat'),2)
        datapts=20;
        molpcell=100;
        r=250;
        w=1000;
        c=cell(datapts,1);
        x=[];
        y=[];
        for i=1:datapts
            c{i}=onecell(molpcell,r,w,'sc',64,[250 1000],1,0);
            x=[x; c{i}.pts(:,1)];
            y=[y; c{i}.pts(:,2)];
        end
        save('test_randomness')
    else
        load('test_randomness')
    end
    analyze(datapts,x,y,c)
end


function analyze(datapts,molx,moly,cells)
bin1=cells{1}.l;
bin2=cells{1}.r*2;
expect=[];
figure(1);
[f2,x2]=hist(molx,1:1:bin1);
bar(x2*10,f2*10/trapz(x2,f2))
hold on;
molx=sort(molx);
% expect(1:datapts)=0;
% for w=1:datapts
%     cel=cells{w};
%     expect(w)=(sum(cel.img(molx(w),:)/sum(sum(cel.img))));
% end
% 
% expect=[molx *10;expect*10];
% plot(expect(1,:),expect(2,:),'color','red')
hold off;
title('Distribution of a Large Number of Randomly Chosen X-values',... 
  'FontWeight','bold')
xlabel('X (nm), Resolution 10nm/pixel')
ylabel('Probability')
saveas(gcf, 'testhistmolx.fig')

expect=[];
figure(2);
[f3,x3]=hist(moly,1:1:bin2);
bar(x3*10,f3*10/trapz(x3,f3))
hold on;
moly=sort(moly);
% expect(1:datapts,1)=0;
% for w=1:datapts
%     cel=cells{w};
%     expect(w)=(sum(cel.img(:,moly(w))/sum(sum(cel.img))));
% end
% expect=expect';
% expect=[moly *10;expect*10];
% plot(expect(1,:),expect(2,:),'color','red')
hold off;
title('Distribution of a Large Number of Randomly Chosen Y-values',... 
  'FontWeight','bold')
xlabel('Y (nm), Resolution 10nm/pixel')
ylabel('Probability')
saveas(gcf, 'testhistmoly.fig')


end