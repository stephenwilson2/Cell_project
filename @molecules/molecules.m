classdef molecules
    %MOLECULES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (SetAccess=private)
        numofmol=1
        x=zeros(1,1);
        y=zeros(1,1);
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
       
       function obj = addmolecules(obj,cel)
            for i=1:obj.numofmol
                a=1;
                b=cel.r*2;
                c=cel.l*2;
                
                obj.x(i)=a + (b-a).*rand(1);
                obj.y(i)=a + (c-a).*rand(1);

                while ~cel.incell(obj.x(i),obj.y(i))
                    obj.x(i)=a + (b-a).*rand(1);
                    obj.y(i)=a + (c-a).*rand(1);
                end
            end
       end  %addmolecules
    end %methods
    
end %classdef

