clear all;
close all;
x = 100;

for i = 1:100
    data = [zeros(1,x+i)  ones(1,x-i)];
    m(i) = mean(data);
    v(i) = var(data);
end

figure;
plot(m,v,'o');

xlabel('Mean Pixel Intensity')
ylabel('Variance')