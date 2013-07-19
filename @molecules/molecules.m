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
        siz=10; %default length of the molecule
        kuhnlength=10;%length of each segment in nm
    end %properties
    
    
    methods
       function obj = molecules(cell,c)
           %Constructs molecules for the onecell class
            if isa(c,'molecules')
                obj.x=c.x;
                obj.y=c.y;
                obj.z=c.z;
            else
                obj.numofmol=c;
                obj.x=zeros(obj.numofmol,1);
                obj.y=zeros(obj.numofmol,1);
                obj.z=zeros(obj.numofmol,1);
                obj=obj.addmolecules(cell);
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
       
      
       function obj=rnaodna(obj)
           molx=zeros(length(obj.x),obj.siz);
           moly=zeros(length(obj.y),obj.siz);
           molz=zeros(length(obj.z),obj.siz);
           for i=1:obj.siz
               ang=rand(1,length(obj.x))*2*pi;
               ang2=rand(1,length(obj.x))*2*pi;

               molx(:,i)=obj.x+cos(ang')*obj.kuhnlength;
               moly(:,i)=obj.y+sin(ang')*obj.kuhnlength;
               molz(:,i)=obj.z+sin(ang2');
           end
           obj.x=reshape(molx,1,length(obj.x)*obj.siz);
           obj.y=reshape(moly,1,length(obj.y)*obj.siz);
           obj.z=reshape(molz,1,length(obj.z)*obj.siz);
       end
       
       
    end %methods
    
end %classdef

