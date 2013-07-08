classdef molecules
    %MOLECULES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
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
       end
       function obj = addmolecules(obj,cel)
            for i=1:obj.numofmol
                obj.x(i)=randi([1 cel.r]);
                obj.y(i)=randi([1 cel.l]);

                while ~cel.incell(obj.x(i),obj.y(i))
                    obj.x(i)=randi([1 cel.r]);
                    obj.y(i)=randi([1 cel.l]);
                end
            end
       end   
       function val = get.x(obj)
           val=obj.x;
       end
       function val = get.y(obj)
           val=obj.y;
       end
    end %methods
    
end %classdef

