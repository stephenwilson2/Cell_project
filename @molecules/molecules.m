classdef molecules
    %MOLECULES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        numofmol=1
        x
        y
    end %properties
    
    
    methods
       function obj = molecules(cell,c)
            if isa(c,'molecules')
                obj.x=c.x;
                obj.y=c.y;
            else
                obj.numofmol=c;
                obj.x=zeros(obj.numofmol,1);
                obj.y=zeros(obj.numofmol,1);
                obj=obj.addmolecules(cell);
            end
       end %constructor
%        function obj=set.x(obj,val)
%           obj.x=val; 
%        end
       
       function obj = addmolecules(obj,cel)  
           a=1;
           b=cel.l*2;
           c=cel.r*2;
           obj.x=zeros(obj.numofmol,1);
           obj.y=zeros(obj.numofmol,1);
           for i=1:obj.numofmol
               obj.x(i)=a + (b-a).*rand(1);
               obj.y(i)=a + (c-a).*rand(1);
               
               while ~cel.incell(obj.x(i),obj.y(i))
                   obj.x(i)=a + (b-a).*rand(1);
                   obj.y(i)=a + (c-a).*rand(1);
               end
               obj.x(i)=obj.x(i)+cel.ori(2)-cel.l;
               obj.y(i)=obj.y(i)+cel.ori(1)-cel.r;
           end
       end  %addmolecules
    end %methods
    
end %classdef

