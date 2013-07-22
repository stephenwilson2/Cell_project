classdef molecules
    %MOLECULES Constructs molecules for the onecell class
    %   Takes the onecell object and one of these:
	%      - An old molecules object
	%      - Number of molecules to add

    
    properties (SetAccess=private)
        numofmol=1 %Number of molecules
        x %X-coordinates of the molecules
        y %Y-coordinates of the molecules
        z %Z-coordinates of the molecules
        type
        siz=10; %default length of the molecule
        kuhnlength=10;%length of each segment in nm
    end %properties
    
    
    methods
       function obj = molecules(cell,c,varargin)
           %Constructs molecules for the onecell class
            if isa(c,'molecules')
                obj.x=c.x;
                obj.y=c.y;
                obj.z=c.z;
            else
                if isempty(varargin)
                    obj.type='pts';
                else
                    obj.type=varargin{1};
                end
                obj.numofmol=c;
                obj.x=zeros(obj.numofmol,1);
                obj.y=zeros(obj.numofmol,1);
                obj.z=zeros(obj.numofmol,1);
                obj=obj.addmolecules(cell);
                if strcmpi(obj.type,'dna')||strcmpi(obj.type,'protein')||strcmpi(obj.type,'rna')
                    obj=obj.randomwalk(cell);
                end
            end
       end %constructor
       
       function obj = addmolecules(obj,cel)
           %Calculates random origins for the molecules
           a=1;
           b=cel.l*2;
           c=cel.r*2;
%            obj.x=zeros(obj.numofmol,1);
%            obj.y=zeros(obj.numofmol,1);
           for i=1:obj.numofmol
               obj.x(i)=a + (b-a).*rand(1);
               obj.y(i)=a + (c-a).*rand(1);
               obj.z(i)=a + (c-a).*rand(1);
               
               while ~cel.incell(obj.x(i),obj.y(i),obj.z(i))
                   obj.x(i)=a + (b-a).*rand(1);
                   obj.y(i)=a + (c-a).*rand(1);
                   obj.z(i)=a + (c-a).*rand(1);
               end
               if strcmp(cel.algo,'sc')
                   obj.x(i)=obj.x(i)-cel.ori(2);
                   obj.y(i)=obj.y(i)-cel.ori(1);
                   obj.z(i)=obj.z(i)-cel.ori(3);
               else
                   obj.x(i)=obj.x(i)+cel.ori(2)-cel.l;
                   obj.y(i)=obj.y(i)+cel.ori(1)-cel.r;
                   obj.z(i)=obj.z(i)+cel.ori(3)-cel.r;
               end
               
           end
           
       end  %addmolecules
       
      
       function obj=randomwalk(obj,cel)
           % randomwalk (NOT WORKING)Creates a random walk from the start point of each molecule and uses
           % the same shape for the molecules so they are identical copies
           % at different locations.
           molx=zeros(length(obj.x),obj.siz);
           moly=zeros(length(obj.y),obj.siz);
           molz=zeros(length(obj.z),obj.siz);
           for i=1:length(obj.x)
               i
               while ~cel.incell(molx(i,1),moly(i,1),molz(i,1))
                   zang=randn(1,length(molx))*2*pi;
                   ang=randn(1,length(molx))*pi;
                   obj.z
                   molz
                   sin(zang)*obj.kuhnlength
                   obj.z(i)+sin(zang)*obj.kuhnlength
                   hey=obj.z(i)+sin(zang)*obj.kuhnlength;

                   molz(i,:)=hey;
                   t=cos(zang)*obj.kuhnlength;
                   molx(i,:)=obj.x(i)+t.*sin(ang);
                   moly(i,:)=obj.y(i)+t.*sin(ang);
                   [molx moly molz]
               end
           end
           obj.x=reshape(molx,length(obj.x)*obj.siz,1);
           obj.y=reshape(moly,length(obj.y)*obj.siz,1);
           obj.z=reshape(molz,length(obj.z)*obj.siz,1);
       end
       
       
    end %methods
    
end %classdef

