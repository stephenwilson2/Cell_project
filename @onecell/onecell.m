classdef onecell
    %onecell Creates one cell that can be manipulated.
    %   WILL ADDD DETAILS HERE

    properties
        r=1000;
        l=1000;
        numofmol=1;
        angle=0; %in degrees for some reason
        algo='c';
        ori
        pixelsize=64;
        
    end
    
    properties (SetAccess=private)
        mol
        img=[];
        pts=[];
        fl=[];
        current=1;
        gopsf=1;
    end
    
    
    % Class methods
    methods
        function obj = onecell(varargin)
            % Sets defaults for optional inputs in order: numofmol,r,l,algo,pixelsize,ori,angle
            optargs = {obj.numofmol,obj.r,obj.l,obj.algo,obj.pixelsize,obj.ori,obj.gopsf,obj.angle};
            
            % Checks to ensure 8 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 8
                error('Takes at most 8 optional inputs');
            end
            
            % Overwrites defaults if optional input exists
            optargs(1:numvarargs) = varargin;
            obj.numofmol= cell2mat(optargs(1));
            obj.r = cell2mat(optargs(2));
            obj.l = cell2mat(optargs(3));
            obj.algo = cell2mat(optargs(4));
            obj.pixelsize = cell2mat(optargs(5));
            if isempty(cell2mat(optargs(6)))
                obj.ori = [obj.r,obj.l];
            else
                obj.ori = cell2mat(optargs(6));
            end
            obj.gopsf = cell2mat(optargs(7));
            obj.angle = cell2mat(optargs(8));
            
            % Construct a onccell object
            if obj.algo=='c'
                obj.l=obj.r;
            elseif strcmp(obj.algo,'sc')
                obj.ori=obj.ori-[obj.r,obj.l];
                obj.l=obj.l*2;
                
            end
            
            if strcmp(obj.algo,'sc') && obj.l<obj.r*3
                error('Spherocylinders need to be long... Increase the length of the cell to at least 3 times the length.');
            end
            obj=refresh_all(obj);
        end %onecell
        
        function obj=refresh_all(obj)
            obj=obj.addMolecules(obj.numofmol);
            obj=label(obj);
            if obj.gopsf==1
                if strcmp(obj.algo,'sc')
                    obj.l=obj.l/2;
                    obj=applyPSF(obj);
                    obj.l=obj.l*2;
                else
                    obj=applyPSF(obj);
                end                
            end
            obj=rotate(obj);
            obj.current=1;
        end %refresh_all
        
        function obj=refresh_psf(obj)
            if obj.gopsf==1
                if strcmp(obj.algo,'sc')
                    obj.l=obj.l/2;
                    obj=applyPSF(obj);
                    obj.l=obj.l*2;
                else
                    obj=applyPSF(obj);
                end                
            end
            obj=rotate(obj);
            obj.current=1;
        end %refresh_psf
        
        %set information about the cell
        function obj = set.ori(obj,val)
            if ~isa(val,'double')
                error('Origin must be of class double')
            end
            obj.ori(1)=val(1);
            obj.ori(2)=val(2);
            obj.current=0; %#ok<*MCSUP>
        end % set.ori
        function obj = set.r(obj,val)
            if ~isa(val,'double')
                error('Radius must be of class double')
            end
            obj.r = val;
            obj.current=0;
        end % set.r
        function obj = set.pixelsize(obj,val)
            if ~isa(val,'double')
                error('Radius must be of class double')
            end
            obj.pixelsize = val;
            obj.current=0;
        end % set.pixelsize
        function obj = set.numofmol(obj,val)
            if ~isa(val,'double')
                error('Radius must be of class double')
            end
            obj.numofmol = val;
            obj.current=0;
        end % set.r
        function obj = set.l(obj,val)
            if ~isa(val,'double')
                error('Length must be of class double')
            end
            obj.l = val;
            obj.current=0;
        end % set.l
        function obj = set.angle(obj,val)
            if ~isa(val,'double')
                error('Angle must be of class double')
            end
            obj.angle = val;
            obj.current=0;
        end % set.angle
        
        %get information about the cell
        function val = get.l(obj)
            val=obj.l;
        end % get.l
        function val = get.r(obj)
            val=obj.r;
        end % get.r
        function val = get.angle(obj)
            val=obj.angle;
        end % get.angle
        function val = get.ori(obj)
            val=obj.ori;
        end % get.ori
        
        function val = incell(obj,x,y)
            %INCELL incell takes an x,y position and determines if it is in
            %the cell or not.
            %It returns a 0 if it is not and a 1 if it is.
            if strcmp(obj.algo,'s')
                X=asin((x)/obj.r);
                Y=acos((y)/obj.l);
                if isreal(X) || isreal(Y)
                    val=1;
                else
                    val=0;
                end
            elseif strcmp(obj.algo,'c')
                 X=asin((abs(x-obj.r)^2+abs(y-obj.r)^2)^0.5/obj.r);
                 if isreal(X)
                    val=1;
                else
                    val=0;
                 end
            elseif strcmp(obj.algo,'sc')
                X=asin((abs(x-obj.r)^2+abs(y-obj.r)^2)^0.5/obj.r);%the circle closer at r
                Y=asin((abs(x-obj.l+obj.r)^2+abs(y-obj.r)^2)^0.5/obj.r);%the cirlce closer l-r
                X1=(asin((x-obj.r)/(obj.l-obj.r*2)))^0.5;
                Y2=(acos((y-obj.l-obj.r)/(2*obj.r))^0.5);
                if (isreal(X) || isreal(Y))||(isreal(X1) || isreal(Y2))
                    val=1;
                else
                    val=0;
                end
            end
            
        end %incell(x,y)
        
        %Adjusts or adds to the cell
        function obj = addMolecules(obj,val)
            %addMolecules Adds molecules randomly into the cell according to
            %the number given
            obj.numofmol=val;
            obj.mol=molecules(obj,val);
            obj.pts=[obj.mol.x obj.mol.y];
            obj.current=0;
        end
        function obj = label(obj)
            tmp=labels(obj,obj.mol);
            obj.fl=tmp.flpts;
            obj.current=0;
        end
        function obj = applyPSF(obj)
            obj.img=psf(obj,obj.pixelsize).img;
            obj.current=0;
%             obj.img
        end
        function obj = rotate(obj)
            obj.img=imrotate(obj.img, obj.angle);
        end
        
        %turn back on after diagnostics
        
        %       function disp(obj)
        %          % DISP Display object in MATLAB syntax
        %  
        %       end % disp
        %Show the cell
        function imshow(obj)
            % IMSHOW Shows the cell
            if obj.current==0;
                reply = input('Cell not current, refresh cell? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'Y';
                end
                if strcmp(reply,'Y')
                    obj=obj.refresh();
                end
            end
            if isempty(obj.img)
                plot(obj.fl(:,1),obj.fl(:,2),'o');
            else
               imshow(mat2gray(flipud(obj.img')));
            end
        end %imshow
        function imagesc(obj)
            % IMAGESC Shows the cell
            if obj.current==0;
                reply = input('Cell not current, refresh cell? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'Y';
                end
                if strcmp(reply,'Y')
                    obj=obj.refresh();
                end
            end
            if isempty(obj.img)
                plot(obj.fl(:,1),obj.fl(:,2),'o');
            else
               colormap('Jet');imagesc(obj.img');
            end
        end %imagesc
        function plot(obj)
            % PLOT  PLOT(obj) plots the cell
            if obj.current==0;
                reply = input('Cell not current, refresh cell? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'Y';
                end
                if strcmp(reply,'Y')
                    obj=obj.refresh();
                end
            end
%             (obj.fl(:,1)-obj.ori(1))-(1-cos(obj.angle*pi/180))
            if strcmp(obj.algo,'sc')
                plot(obj.fl(:,1)/obj.pixelsize+obj.l/2/obj.pixelsize*.3,obj.fl(:,2)/obj.pixelsize+obj.r/2/obj.pixelsize*.3,'o');
            else
                plot(obj.fl(:,1)/obj.pixelsize+obj.l/obj.pixelsize*.3,obj.fl(:,2)/obj.pixelsize+obj.r/obj.pixelsize*.3,'o');
            end
%             if strcmp(obj.algo,'sc')
%                 plot(obj.fl(:,1)/obj.pixelsize,obj.fl(:,2)/obj.pixelsize,'o');
%             else
%                 plot(obj.fl(:,1)/obj.pixelsize,obj.fl(:,2)/obj.pixelsize,'o');
%             end
        end % plot
        
    end % methods
end % classdef