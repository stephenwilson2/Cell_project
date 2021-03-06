classdef onecell
    %%onecell Creates one cell that can be manipulated.
    %Constructor:
    %Takes input for any of the following variables:
    %    - numofmol: the number of molecules to be added to the cell
    %    - r: The radius of the cell
    %    - l: The length of the cell (this is half the actual length)
    %    - algo: Selects for the type of cell to be used. Inputs are:
    %        - 's': sphere
    %        - 'b': box (or rectangle)
    %        - 'sc': spherocylinder
    %    - pixelsize: Gives the size in nm^2 of the camerapixel
    %    - ori: States the origin of the cell in nm. 
    %        - 'b' and 's' cells are set to [r,l] by default
    %        - 'sc' cells are set to [0,0] by default
    %    - gopsf: If set to 1, the PSF is applied to the fluorophores. If set to 0, the PSF is not set
    %    - sections: This is the number of the z-sections that will be automatically calculated
    %    - angle: Gives the angle to rotate the cell in degrees
    %
    %Variables:
    %User changeable variables are the same as those in the constructor and can be changed in the following manner:
    %    c=onecell();
    %    c.r=20;
    %The radius would then be 20 nm.
    %
    %The following variables are accessible to the user but cannot be directly changed:
    %    - mol: This stores the molecules object. There should be functions in onecell to manipulate the molecules object.
    %    - img: Stores cells that contain the image matricies
    %    - pts: Stores the locations of the molecules as pairs of points
    %    - fl: Stores the locations of the fluorophores as pairs of points
    %    - cellmask: Stores an image of the cell's mask. Must be constructed with the cell_mask function
    %    - current: If 1, then the img, pts, and fl variables are up-to-date. If 0,  one of onecell's refresh functions need to be used
    %
    %Functions:
    %addMolecules (int): Adds a given integer number molecules to the cell
    %applyPSF ():  Applies the PSF and creates onecell's img
    %cell_mask(): Creates a cell mask stored under onecell as cellmask
    %imagesc(onecell): Specifies the way that imagesc displays the onecell object. Uses a color bar and labels
    %imshow(onecell): Specifies the way that imshow displays the onecell object
    %incell(double x, double y):  Checks to see if point [x,y] is in the cell
    %label(): Labels the molecules with flurophores
    %Onecell(varargin): onecell's constructor 
    %plot(onecell): Specifies the way that plot displays the onecell object. When plotted the fluorophores are 
    %             plotted so that they could be overlaid onto a cell
    %refresh_all(): Refreshes the entire cell; Used if the molecules or fluorophores were changed
    %refresh_cell(): Refreshes everything but does not get new positions for the molecules or change the 
    %                number of molecules. Used if anything but the molecules was changed
    %rotate(): Rotates the cell according to onecell's angle

    properties
        r=250; % The radius of the cell
        l=1000; % The length of the cell (this is half the actual length)
        numofmol=10; %The number of molecules to be added to the cell
        angle=0; %Gives the angle to rotate the cell in degrees
        algo='sc'; %Selects for the type of cell to be used. Inputs are: - 's': sphere - 'b': box (or rectangle)- 'sc': spherocylinder
        ori %States the origin of the cell in nm. - 'b' and 's' cells are set to [r,l,r] by default - 'sc' cells are set to [r,l,r] by default
        pixelsize=64;%Gives the size in nm^2 of the camerapixel
        gopsf=1;%If set to 1, the PSF is applied to the fluorophores. If set to 0, the PSF is not set
        sections=8;%This is the number of the z-sections that will be automatically calculated
        PSF % Stores the PSF
        moltype %describes the type of molecule to be built
    end
    
    properties (SetAccess=private)
        mol %This stores the molecules object. There should be functions in onecell to manipulate the molecules object.
        img; %Stores cells that contain the image matricies
        pts=[]; %Stores the locations of the molecules as pairs of points
        fl=[]; %Stores the locations of the fluorophores as pairs of points
        cellmask=[]; %Stores an image of the cell's mask. Must be constructed with the cell_mask function
        current=1; %If 1, then the img, pts, and fl variables are up-to-date. If 0,  one of onecell's refresh functions need to be used
    end
    
    properties(SetAccess=private,GetAccess=private)
       oldr
       oldl
       oldpixelsize
       oldori
    end
    
    methods(Static)
        function bool=sphere(obj,x,y,z,varargin)
            if isempty(varargin)
                r=obj.r;
            else
                r=varargin{1};
            end
            X=asin((abs(x-r)^2+abs(y-obj.r)^2)^0.5/obj.r);
            Z=asin((abs(x-r)^2+abs(z-obj.r)^2)^0.5/obj.r);
            Z2=asin((abs(y-obj.r)^2+abs(z-obj.r)^2)^0.5/obj.r);
            if isreal(X) && isreal(Z) && isreal(Z2)
                bool=1;
            else
                bool=0;
            end                
        end
        function bool=cylinder(obj,x,y,z,varargin)
            if isempty(varargin)
                r=obj.r;
            else
                r=varargin{1};
            end
            X=asin((abs(x-r)^2+abs(y-obj.r)^2)^0.5/(obj.l/2-obj.r));
            Z=asin((abs(x-r)^2+abs(z-obj.r)^2)^0.5/(obj.l/2-obj.r));
            Z2=asin((abs(y-obj.r)^2+abs(z-obj.r)^2)^0.5/obj.r);
            if isreal(X) && isreal(Z) && isreal(Z2)
                bool=1;
            else
                bool=0;
            end      
        end
        function bool=box(obj,x,y,z,varargin)
            if isempty(varargin)
                shift=0;
            else
                shift=varargin{1};
            end
            X=(asin((x-shift)/(obj.l)))^0.5;
            X2=(asin((x+shift)/(obj.l)))^0.5;
            Y=(asin((y)/(obj.r*2)))^0.5;
            Z=(asin((z)/(obj.r*2)))^0.5;
            if isreal(X) && isreal(Y) && isreal(Z) && isreal(X2)
                bool=1;
            else
                bool=0;
            end
        end
    end
    
    %% Class methods
    methods
        function obj = onecell(varargin)
            % Sets defaults for optional inputs in order: numofmol,r,l,algo,pixelsize,ori,gopsf,angle
            optargs = {obj.numofmol,obj.r,obj.l,obj.algo,obj.pixelsize,obj.ori,obj.gopsf,obj.sections,obj.angle};
            
            % Checks to ensure 9 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 9
                error('Takes at most 9 optional inputs');
            end
            
            % Overwrites defaults if optional input exists
            optargs(1:numvarargs) = varargin;
            obj.numofmol= cell2mat(optargs(1));
            obj.r = cell2mat(optargs(2));
            obj.l = cell2mat(optargs(3));
            obj.algo = cell2mat(optargs(4));
            obj.pixelsize = cell2mat(optargs(5));
            if isempty(cell2mat(optargs(6)))
                obj.ori = [obj.r,obj.l,obj.r];
            else
                obj.ori = cell2mat(optargs(6));                
            end
            obj.gopsf = cell2mat(optargs(7));
            obj.sections = cell2mat(optargs(8));
            obj.angle = cell2mat(optargs(9));
            
            % Construct a onccell object
            if obj.algo=='s'
                obj.l=obj.r;
            elseif strcmp(obj.algo,'sc')
                obj.ori=obj.ori-[obj.r,obj.l,obj.r];
                obj.l=obj.l*2;
            end
            
            if strcmp(obj.algo,'sc') && obj.l<obj.r*3
                error('Spherocylinders need to be long... Increase the length of the cell to at least 3 times the length.');
            end
            obj=obj.addMolecules(obj.numofmol);
        end %onecell
        
        function obj=refresh_all(obj)
            %refresh_all Refreshes the entire cell; Used if the molecules
            %were changed.
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
            obj.oldr=obj.r;
            obj.oldl=obj.l;
            obj.oldori=obj.ori;
            obj.current=1;
        end %refresh_all
        
        function obj=refresh_cell(obj)
            %refresh_all Refreshes everything but does not get new positions for the molecules or change the 
            %number of molecules. Used if anything but the molecules was
            %changed. Not automated because refreshing takes a
            %significant amount of time
            f1=(obj.r-obj.oldr)+(obj.ori(1)-obj.oldori(1));
            f2=(obj.l-obj.oldl)+(obj.ori(2)-obj.oldori(2));
            obj.pts(:,1)=obj.pts(:,1)+f1;
            obj.pts(:,2)=obj.pts(:,2)+f2;
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
            obj.oldr=obj.r;
            obj.oldl=obj.l;
            obj.oldori=obj.ori;
            obj.current=1;
        end %refresh_cell
        
        %%Set information about the cell
        function obj = set.ori(obj,val)
            obj.oldori=obj.ori;
            if ~isa(val,'double')
                error('Origin must be of class double')
            end
            obj.ori(1)=val(1);
            obj.ori(2)=val(2);
            obj.ori(3)=val(3);
            obj.current=0; %#ok<*MCSUP>
        end % set.ori
        function obj = set.r(obj,val)
            obj.oldr=obj.r;
            if ~isa(val,'double')
                error('Radius must be of class double')
            end
            obj.r = val;
            obj.current=0;

        end % set.r
        function obj = set.pixelsize(obj,val)
            obj.oldpixelsize=obj.pixelsize;
            if ~isa(val,'double')
                error('pixelsize must be of class double')
            end
            obj.pixelsize = val;
            obj.current=0;
        end % set.pixelsize
        function obj = set.numofmol(obj,val)
            if ~isa(val,'double')
                error('numofmol must be of class double')
            end
            obj.numofmol = val;
            obj.current=0;
        end % set.numofmol
        function obj = set.l(obj,val)
            obj.oldl=obj.l;
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
        
        function obj = set.moltype(obj,val)
            if ~isa(val,'char')
                error('moltype must be of class char')
            end
            obj.moltype = val;
            obj=obj.addMolecules(obj.numofmol);
            obj=refresh_cell(obj);
        end % set.moltype
        
        function obj = set.sections(obj,val)
            if ~isa(val,'double')
                error('Sections must be of class double')
            end
            obj.sections = val;
            obj.current=0;
        end % set.sections
        
        function obj = check(obj)
            % check Checks to see if the cell needs to be refreshed
            if obj.current==0;
                reply = input('Cell not current, refresh cell? Y/N [Y]: ', 's');
                if isempty(reply)
                    reply = 'Y';
                end
                if strcmpi(reply,'Y')
                    obj=obj.refresh_cell();
                    obj.current
                end
            end
        end
        
        %%get information about the cell
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
        
        function val = incell(obj,x,y,z)
            %INCELL Checks to see if point [x,y,z] is in the cell
            %the cell or not.
            %It returns a 0 if it is not and a 1 if it is.
            if strcmp(obj.algo,'b')
                val= obj.box(obj,x,y,z);
            elseif strcmp(obj.algo,'s')
                val=obj.sphere(obj,x,y,z);
            elseif strcmp(obj.algo,'sc')
                val1=obj.sphere(obj,x,y,z,obj.r);%the sphere closer at r
                val2=obj.sphere(obj,x,y,z,obj.l-obj.r);
                val3=obj.cylinder(obj,x,y,z,obj.l/2);
                if  val1 || val2 ||val3
                    val=1;
                else
                    val=0;
                end
            end
        end %incell(x,y)
        
        
        %%Adjusts or adds to the cell
        function obj = addMolecules(obj,val)
            %addMolecules Adds a given integer number molecules to the cell
            obj.numofmol=val;
            if ~isempty(obj.moltype)
                obj.mol=molecules(obj,val,obj.moltype);
            else
                obj.mol=molecules(obj,val);
            end
            obj.pts=[obj.mol.x obj.mol.y obj.mol.z];
            obj=refresh_all(obj);
        end
        function obj = label(obj)
            % label Labels the molecules with flurophores
            tmp=labels(obj,obj.mol);
            obj.fl=tmp.flpts;
            obj.current=0;
        end
        function obj = applyPSF(obj)
            %applyPSF Applies the PSF and creates onecell's img
            obj.PSF=psf(obj,obj.pixelsize,obj.sections);
            obj.img=obj.PSF.img;
            obj.current=0;
        end
        function obj = rotate(obj)
            for i=1:length(obj.img)
                if ~isempty(obj.img{i})
                    obj.img{i}=imrotate(obj.img{i}, obj.angle);
                end
            end
        end
        function obj = cell_mask(obj)
            %cell_mask Creates a cell mask stored under onecell as cellmask
            if strcmp(obj.algo,'sc')
                obj.cellmask=zeros(round(obj.l/obj.pixelsize),round(obj.r*2/obj.pixelsize)/2);
            else
                error('Only Spherocylinders currently supported')
            end
            pt=[];
            for x=1:1:obj.r
                for y=1:1:obj.r
                    X=asin((abs(x-obj.r)^2+abs(y-obj.r)^2)^0.5/obj.r);%the sphere closer at r
                    Y=asin((abs(x-obj.l+obj.r)^2+abs(y-obj.r)^2)^0.5/obj.r);%the cirlce at l-r
                    X1=(asin((x-obj.r)/(obj.l-obj.r*2)))^0.5;
                    Y2=(acos((y-obj.l-obj.r)/(2*obj.r))^0.5);
                    if (isreal(X) || isreal(Y))||(isreal(X1) || isreal(Y2))
                        pt=[pt;[x y]];
                    end
                end
            end
            x=pt(:,1);y=pt(:,2);
            mx=obj.r+1:obj.l-obj.r;
            x=[x;(x+obj.l-obj.r);mx';mx'];
            my=zeros(obj.l-2*obj.r,1);
            my(:)=obj.r;
            my2=ones(obj.l-2*obj.r,1);

            y=[y;flipud(y);my2;my];
            
            for i=1:length(x)
                obj.cellmask(ceil(x(i)/obj.pixelsize),ceil(y(i)/obj.pixelsize))=1;
            end
            tmp=obj.cellmask;
            
            obj.cellmask=zeros(round(obj.l/obj.pixelsize),round(obj.r*2/obj.pixelsize));
            obj.cellmask(1:size(tmp,1),1:size(tmp,2))=tmp;
            obj.cellmask(1:size(tmp,1),size(tmp,2)+1:2*size(tmp,2))=fliplr(tmp);
            obj.cellmask=imfill(obj.cellmask,8,'holes');
            
        end

        %%Show the cell
        function obj=imshow(obj)
            % IMSHOW Specifies the way that imshow displays the onecell object
            obj=check(obj);
            if isempty(obj.img)
                plot3(obj.fl(:,1),obj.fl(:,2),obj.fl(:,3),'o');axis equal;axis tight;
            else
                n=0;                
                if obj.current==1
                    r=obj.r; %#ok<*PROP>
                    l=obj.l;
                    px=obj.pixelsize;
                else
                    px=obj.oldpixelsize;
                    r=obj.oldr;
                    l=obj.oldl;
                end
                if strcmp(obj.algo,'s')
                    l=l*2;
                end
                mi=1;
                ma=0;
                for i=1:length(obj.img)
                    if ~isempty(obj.img{i})
                        tmin=min(min(obj.img{i}));
                        tmax=max(max(obj.img{i}));
                        if tmin<mi
                            mi=tmin;
                        end
                        if tmax>ma
                            ma=tmax;
                        end
                    end
                end
               planes=linspace(1,obj.r*2,obj.sections);
               for i=1:length(obj.img)
                    if ~isempty(obj.img{i})
                        n=n+1;
                        figure(1);
                        subplot(obj.sections,1,n)
                        imshow(flipud(obj.img{i}'),[mi ma]);colormap(gray);colorbar;axis equal;axis tight;
                        title(sprintf('%i molecules in a %i nm by %i nm by %i nm cell\n Depth: %i nm Resolution: %i^2 nm^2 / pixel',obj.numofmol,r*2,l,r*2,round(planes(n)),px),'FontWeight','bold');
                    end
               end
            end
        end %imshow
        function obj=showslice(obj,num,varargin)
            if ~isempty(varargin)
                m=varargin{1};
            end
            %showslice This function displays a single z-slice at the given
            %depth.
            if obj.current==1
                r=obj.r; %#ok<*PROP>
                l=obj.l;
                px=obj.pixelsize;
            else
                px=obj.oldpixelsize;
                r=obj.oldr;
                l=obj.oldl;
            end
            if strcmp(obj.algo,'s')
                    l=l*2;
            end
            planes=linspace(1,obj.r*2,obj.sections);
            n=0;
            for i=1:length(obj.img)
                if ~isempty(obj.img{i})
                    n=n+1;
                    if n==num
                        if isempty(varargin)
                            m=max(max(obj.img{num}));
                        end
                        imshow(flipud(obj.img{num}'),[0 m]);colormap(gray);colorbar;axis equal;axis tight;
                        title(sprintf('%i molecules in a %i nm by %i nm by %i nm cell\n Depth: %i nm Resolution: %i^2 nm^2 / pixel',obj.numofmol,r*2,l,r*2,round(planes(n)),px),'FontWeight','bold');
                    end
                end
            end
        end
        function obj=surf(obj)
            %surf Displays the onecell object's slices in 3-D
            obj=check(obj);
            if isempty(obj.img)
                plot3(obj.fl(:,1),obj.fl(:,2),obj.fl(:,3),'o');axis equal;axis tight;
            else
                n=0;                
                if obj.current==1
                    r=obj.r; %#ok<*PROP>
                    l=obj.l;
                    px=obj.pixelsize;
                else
                    px=obj.oldpixelsize;
                    r=obj.oldr;
                    l=obj.oldl;
                end
                if strcmp(obj.algo,'s')
                    l=l*2;
                end
                mi=1;
                ma=0;
                for i=1:length(obj.img)
                    if ~isempty(obj.img{i})
                        tmin=min(min(obj.img{i}));
                        tmax=max(max(obj.img{i}));
                        if tmin<mi
                            mi=tmin;
                        end
                        if tmax>ma
                            ma=tmax;
                        end
                    end
                end
                planes=linspace(1,obj.r*2,obj.sections);
                for i=1:length(obj.img)
                    if ~isempty(obj.img{i})
                        n=n+1;
                        figure(1);
                        subplot(4,obj.sections/4,n)
                        colormap(gray);
                        [x,y] = meshgrid(1:1:size(obj.img{i},2),1:1:size(obj.img{i},1));
                        z=obj.img{i};
                        surf(x,y,z)
                        view([100,64]);axis tight;axis([1 size(obj.img{1},2) 1 size(obj.img{1},1) mi ma]);
                        title(sprintf('%i molecules in a %i nm by %i nm by %i nm cell\n Depth: %i nm Resolution: %i^2 nm^2 / pixel',obj.numofmol,r*2,l,r*2,round(planes(n)),px),'FontWeight','bold');
                    end
               end
            end          
        end
        function obj=imagesc(obj)
            % IMAGESC Specifies the way that imagesc displays the onecell object. Uses a color bar and labels
            obj=check(obj);
            if isempty(obj.img)
                plot3(obj.fl(:,1),obj.fl(:,2),obj.fl(:,3),'o');
            else
               if obj.current==1
                   r=obj.r; %#ok<*PROP>
                   l=obj.l;
                   px=obj.pixelsize;
               else
                   px=obj.oldpixelsize;
                   r=obj.oldr;
                   l=obj.oldl;
               end
               if strcmp(obj.algo,'s')
                    l=l*2;
               end
               n=0;
               planes=linspace(1,obj.r*2,obj.sections);
               for i=1:length(obj.img)
                    if ~isempty(obj.img{i})
                        n=n+1;
                        figure(1);
                        subplot(obj.sections,1,n)
                        imagesc(obj.img{i}');colormap(gray);colorbar;axis equal;axis tight;
                        title(sprintf('%i molecules in a %i nm by %i nm by %i nm cell\n Depth: %i nm Resolution: %i^2 nm^2 / pixel',obj.numofmol,r*2,l,r*2,round(planes(n)),px),'FontWeight','bold');
                    end
               end
 
            end
        end %imagesc
        function obj=plot(obj)
            % PLOT  Specifies the way that plot displays the onecell object. When plotted the fluorophores are plotted so that they could be overlaid onto a cell
 
            %             (obj.fl(:,1)-obj.ori(1))-(1-cos(obj.angle*pi/180))
            obj=check(obj);
            if obj.current==1
                r=obj.r; %#ok<*PROP>
                l=obj.l;
                px=obj.pixelsize;
            else
                px=obj.oldpixelsize;
                r=obj.oldr;
                l=obj.oldl;
            end
            if isempty(obj.img)
                if strcmp(obj.algo,'sc')
                    plot3(obj.fl(:,1)/px,obj.fl(:,2)/px,obj.fl(:,3)/px,'o');axis equal;axis tight;
                else
                    plot3(obj.fl(:,1)/px,obj.fl(:,2)/px,obj.fl(:,3)/px,'o');axis equal;axis tight;
                end
            else
                if strcmp(obj.algo,'sc')
                    plot3(obj.fl(:,1)/px+l/2/px*.3,obj.fl(:,2)/px+r/2/px*.3,obj.fl(:,3)/px+r/2/px*.3,'o');axis equal;axis tight;
                else
                    plot3(obj.fl(:,1)/px+l/px*.3,obj.fl(:,2)/px+r/px*.3,obj.fl(:,3)/px+r/px*.3,'o');axis equal;axis tight;
                end
            end
        end % plot
        
    end % methods
end % classdef