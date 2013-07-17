classdef psf
    %PSF Summary of this class goes here
    %   Takes the onecell object and two optional arguments: 
	%     - the pixel size of the camera 
    %     - the z-plane to be simulated
	%     - the wavelength of the emission light
    
    properties (SetAccess=private)
        emwave=520; %nm
        sigma
        sigmaz
        numofplanes
        img={}
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
        function obj = psf(cel,varargin)
            c=varargin;
            if isa(c,'psf')
                obj.emwave=c.emwave;
                obj.sigma=c.sigma;
                obj.img=c.img;
                obj.fl=c.fl;
                obj.pxsz=c.pxsz;
            else
                optargs = {obj.pxsz,obj.numofplanes,obj.emwave};
                % Checks to ensure 3 optional inputs at most
                numvarargs = length(c);
                if numvarargs > 3
                    error('Takes at most 3 optional inputs');
                end
                % Overwrites defaults if optional input exists
                optargs(1:numvarargs) = c;
                obj.pxsz= cell2mat(optargs(1));
                obj.numofplanes= cell2mat(optargs(2));
                obj.emwave= cell2mat(optargs(3));
                
                obj=calcsig(obj);
                obj.fl=cel.fl;
                obj.fl=obj.fl/obj.pxsz; %Scales the pts of fluorescence
                obj.fl(:,1)=obj.fl(:,1)+cel.l/obj.pxsz*.3;
                obj.fl(:,2)=obj.fl(:,2)+cel.r/obj.pxsz*.3;
                obj.fl(:,3)=obj.fl(:,3);
                planes=linspace(1,cel.r*2/obj.pxsz,obj.numofplanes); %rounding prob might come up
                obj.img=cell(round(cel.r*2/obj.pxsz),1);
                if length(planes)>length(obj.img)
                    error('Too many sections for that cell length and pixel size. Reduce the number of sections for that cell');
                end
                for plane=planes
                    obj=obj.applyPSF(cel,plane);
                end
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
            
            num=5*(7^.5)*(1-power(cos(a),3/2));
            de=(6^.5)*n*k*(4*power(cos(a),5)-25*power(cos(a),7/2)+42*power(cos(a),5/2)...
                -25*power(cos(a),3/2)+4)^.5;
            obj.sigmaz=num/de;
        end %calcsig
        
        function obj = applyPSF(obj,cel,currentplane)
            %applyPSF Takes the PSF object, the onecells object, and the
            % z-plane (in nm) plane to be simulated.
            im=zeros(round(cel.l*2/obj.pxsz*1.3),round(cel.r*2/obj.pxsz*1.3));%creates blank matrix of the camera pixel size
            for i=1:size(obj.fl,1)
                tmp=zeros(round(cel.l*2/obj.pxsz*1.3),round(cel.r*2/obj.pxsz*1.3));
                mu=obj.fl(i,1);% value in X
                y=ceil(obj.fl(i,2)); % value in Y
                if y>size(im,2)
                    y=y-1;
                end
                z=ceil(obj.fl(i,3)); % value in Z
                if z>size(im,2)
                    z=z-1;
                end
                norm=1/gaussDistribution(z,z,obj.sigmaz/obj.pxsz);
                %Spread contributions for that pt in the plane
                for x=1:size(im,1)
                    tmp(x,y)=tmp(x,y)+...
                        (gaussDistribution(x,mu,obj.sigma/obj.pxsz))*...
                        gaussDistribution(z,currentplane,obj.sigmaz/obj.pxsz);
                end
                for h=1:size(im,1)
                    o=1:size(im,2);
                    im(h,:)=im(h,:)+...
                            (gaussDistribution(o,obj.fl(i,2),obj.sigma/obj.pxsz)*sum(tmp(h,:)));
                end
                obj.img{round(currentplane)}=im;
            end

        end %applyPSF

    end %methods
    
end

