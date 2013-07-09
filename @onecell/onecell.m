classdef onecell
    %onecell Creates one cell that can be manipulated.
    %   WILL ADDD DETAILS HERE

    properties
        r=1000;
        l=1000;
        numofmol=1;
        angle=0;
        algo='c';
        ori
        pixelsize=64;
    end
    
    properties (SetAccess=private)
        mol
        weird
        img=[];
        pts=[];
        fl=[];
        current=1;
    end
    
    
    % Class methods
    methods
        function obj = onecell(varargin)
            % Sets defaults for optional inputs
            optargs = {obj.numofmol,obj.r,obj.l,obj.algo,obj.pixelsize,obj.ori,obj.angle};
            
            % Checks to ensure 7 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 7
                error('Takes at most 7 optional inputs');
            end
            
            % Overwrites defaults if optional input exists
            optargs(1:numvarargs) = varargin;
            obj.numofmol= cell2mat(optargs(1));
            obj.r = cell2mat(optargs(2));
            obj.l = cell2mat(optargs(3));
            obj.algo = cell2mat(optargs(4));
            obj.pixelsize = cell2mat(optargs(5));
            obj.ori = cell2mat(optargs(6));
            obj.angle = cell2mat(optargs(7));
            
            % Construct a onccell object
            if isempty(obj.ori)
                obj.ori = [obj.r/2,obj.l/2];
            end
            obj=obj.addMolecules(obj.numofmol);
            obj=label(obj);
            obj=applyPSF(obj);
            obj.current=1;
        end %onecell
        
        function obj=refresh(obj)
            obj=obj.addMolecules(obj.numofmol);
            obj=label(obj);
            obj=applyPSF(obj);
            obj.current=1;
        end %Refresh needs to be updated for angles and origins
        
        %set information about the cell
        function obj = set.ori(obj,val)
            if ~isa(val,'double')
                error('Origin must be of class double')
            end
            if ~isempty(val)
                if val(1)~=0 && val(1)<obj.r
                    obj.ori(1)=val(1);
                else
                    error('origins start at 1,1 and end at %s,%s',...
                        num2str(obj.r),num2str(obj.l))
                end
                
                if val(2)~=0 && val(2)<obj.l
                    obj.ori(2)=val(2);
                else
                    error('origins start at 1,1 and end at %s,%s',...
                        num2str(obj.r),num2str(obj.l))
                end
            end
            obj.current=0;
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
%         function val = get.mol(obj)
%             val=[obj.mol.x obj.mol.y];
%         end
        
        function val = incell(obj,x,y)
            %INCELL incell takes an x,y position and determines if it is in
            %the cell or not.
            %It returns a 0 if it is not and a 1 if it is.
            if strcmp(obj.algo,'s')
                X=asin((x-obj.ori(1))/obj.r);
                Y=acos((y-obj.ori(2))/obj.l);
                if isreal(X) && isreal(Y)
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
            
        
        %turn back on after diagnostics
        
        %       function disp(obj)
        %          % DISP Display object in MATLAB syntax
        %          imshow(obj);
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
               imshow(mat2gray(obj.img));
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
               imagesc(obj.img);
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
            if isempty(obj.img)
                plot(obj.fl(:,1),obj.fl(:,2),'o');
            else
                imshow(mat2gray(obj.img));
            end
        end % plot
        
    end % methods
end % classdef