classdef psf
    %PSF Summary of this class goes here
    %   Detailed explanation goes here
    
    
    properties (SetAccess=private)
        emwave=520; %nm
        sigma
        img=[]
        fl
        pxsz=10;
    end
    
    methods (Static)
        function f = gaussDistribution(x, mu, s)
            p1 = -.5 * ((x - mu)/s) .^ 2;
            p2 = (s * sqrt(2*pi));
            f = exp(p1) ./ p2;
        end
    end %Static methods
    
    methods
        function obj = psf(cell,varargin)
            c=varargin;
            if isa(c,'psf')
                obj.emwave=c.emwave;
                obj.sigma=c.sigma;
                obj.img=c.img;
                obj.fl=c.fl;
                obj.pxsz=c.pxsz;
            else
                optargs = {obj.pxsz,obj.emwave};
                % Checks to ensure 2 optional inputs at most
                numvarargs = length(c);
                if numvarargs > 2
                    error('Takes at most 2 optional inputs');
                end
                % Overwrites defaults if optional input exists
                optargs(1:numvarargs) = c;
                obj.pxsz= cell2mat(optargs(1));
                obj.emwave= cell2mat(optargs(2));
                
                obj=calcsig(obj);
                obj.fl=cell.fl;
                obj=obj.applyPSF(cell);
            end
        end %constructor
        
        function obj = calcsig(obj)
            n=1.515; %refractive index for immersion oil
            NA=1.4; %numerical apperature
            a=asin(NA/n);
            k=(2*pi/obj.emwave);
            
            num=4-7*power(cos(a),3/2)+3*power(cos(a),7/2);
            de=7*(1-power(cos(a),3/2));
            obj.sigma=1/n/k*power(num/de,-0.5);
        end %calcsig
        
        function obj = applyPSF(obj,cel)
            obj.img=zeros(cel.l*2/obj.pxsz,cel.r*2/obj.pxsz);%creates blank matrix of the camera pixel size
            obj.fl=obj.fl/obj.pxsz; %Scales the pts of fluorescence

            syms a;
            for i=1:size(obj.fl,1)
                mu=obj.fl(1);
                f=gaussDistribution(a,mu,obj.sigma);
%                 d=zeros(size(obj.img,1));
                for x=1:size(obj.img,1)
                    x1=obj.pxsz*gaussDistribution(x,mu,obj.sigma);
                    myx=x:obj.pxsz/100:x+obj.pxsz;
                    x2=trapz(gaussDistribution(myx,mu,obj.sigma))/10;
                    d(x)=(x1-x2)/x2;
                end
                max(d)
            end
           
        end %applyPSF
        

        
    end %methods
    
end

