classdef Elastic
    %Elastic 3D elastic class
    %   Detailed explanation goes here
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
        finiteDisp = 0;
        
        Young;
        Poisson;
        linear = true;
    end
    
    properties(SetAccess = private)
        strss
    end
    %%
    methods
        %% Construct
        function obj = Elastic(num, props)
            obj.ndm   = num.ndm;
            obj.ndof  = num.ndof;
            obj.strss = zeros(num.str, num.gp, num.el, num.steps+1);
            
            if strcmp(props{1,1} ,'E') && strcmp(props{2,1} ,'v')
                obj.Young  = props{1,2};
                obj.Poisson= props{2,2};
            elseif strcmp(props{1,1} ,'v') && strcmp(props{2,1} ,'E')
                obj.Poisson= props{1,2};
                obj.Young  = props{2,2};
            else
                error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
            end
        end
        %% Epsilon
        function [deps, obj] = computeStrain(obj, gp, el, ~)
            deps = gp.B * el.Uvc;
        end
        %% Sigma
        function [sigma_voigt, obj] = computeCauchy(obj, gp, step)
            n.strss = obj.strss(:, gp.i, gp.iel, step);
            
            sigma_voigt = gp.D*gp.deps;
            
            obj.strss(:, gp.i, gp.iel, step+1) = sigma_voigt;
        end
        %% Tangential stiffness
        function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~)
            E = obj.Young;
            v = obj.Poisson;

            Eh= E/(1-2*v)/(1+v);
            G = 0.5*E/(1+v);
            D =[Eh*(1-v) Eh*v     Eh*v     0 0 0
                Eh*v     Eh*(1-v) Eh*v     0 0 0
                Eh*v     Eh*v     Eh*(1-v) 0 0 0
                0        0        0        G 0 0
                0        0        0        0 G 0
                0        0        0        0 0 G];
            
            ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
            
        end
        %% Element K
        function Kel = computeK_el(~, Kel, gp, ~)
            Kel = Kel + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
        end
        %% Element Fint
        function Fint = computeFint(~, gp, el)
            Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
        end
    end
end

