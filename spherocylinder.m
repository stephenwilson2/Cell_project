function [X Y]=spherocylinder(x,y,a,b,angle,steps)
    beta = -angle * (pi / 180);
    if a>b
        tmp=a;
        a=b;
        b=tmp;
    end   
     
    alpha = linspace(0, 180, steps/4)' .* (pi / 180); %creates a linearly spaced vector from 0 to 180 of length steps.
    X = x + (a * cos(alpha) * cos(beta) - a * sin(alpha) * sin(beta))-(b-a)*sin(beta);
    Y = y + (a * cos(alpha) * sin(beta) + a * sin(alpha) * cos(beta))+(b-a)*cos(beta);
    for n=1:(b-a)*2
        X=[X;X(length(X))+1*sin(beta)];
        Y=[Y;Y(length(Y))-1*cos(beta)];
    end
    alpha = linspace(180, 360, steps/4)' .* (pi / 180); %creates a linearly spaced vector from 180 to 360 of length steps.  
    X2 = x + (a * cos(alpha) * cos(beta) - a * sin(alpha) * sin(beta))+(b-a)*sin(beta);
    Y2 = y + (a * cos(alpha) * sin(beta) + a * sin(alpha) * cos(beta))-(b-a)*cos(beta);
    for n=1:(b-a)*2
        X2=[X2;X2(length(X2))-1*sin(beta)];
        Y2=[Y2;Y2(length(Y2))+1*cos(beta)];
    end
    X=[X;X2];
    Y=[Y;Y2];    

    if nargout==1, X = [X Y]; end
end