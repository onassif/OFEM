classdef PlaneStrain
    %PlaneStrain Computes Cauchy stress based on 2D plane strain
    %   Currently only if ndm=2, ndof=2
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
    end
    %%
    methods
        function obj = PlaneStrain(ndm, ndof)
            obj.ndm  = ndm;
            obj.ndof = ndof;
        end
        %% B
        function B = Compute_B(~,gp)
            dx = gp.dNdX(:,1);
            dy = gp.dNdX(:,2);
            B  =[ dx(1),     0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0
                    0.0,   dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4)
                  dy(1),   dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4)];
        end
        %% Sigma
        function sigma_voigt = Compute_cauchy(obj, gp, props)
            [D, ~] = obj.Compute_tangentstiffness( gp, props);
            if (obj.ndm==2 && obj.ndof ==2)
                sigma_voigt = D *gp.eps;
            end
        end
        %% Tangential stiffness
        function [D, ctan] = Compute_tangentstiffness(~, ~, props)
            
            if strcmp(props{1,1} ,'E') && strcmp(props{2,1} ,'v')
                E = props{1,2};
                v = props{2,2};
            elseif strcmp(props{1,1} ,'v') && strcmp(props{2,1} ,'E')
                v = props{1,2};
                E = props{2,2};
            else
                error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
            end
            
            Eh= E/(1-2*v)/(1+v);
            G = 0.5*E/(1+v);
            c =[Eh*(1-v)    Eh*v        Eh*v        0 0 0
                Eh*v        Eh*(1-v)    Eh*v        0 0 0
                Eh*v        Eh*v        Eh*(1-v)    0 0 0
                0           0           0           G 0 0
                0           0           0           0 G 0
                0           0           0           0 0 G];
            
            
            D = [ c(1,1) c(1,2) c(1,4)
                c(2,1) c(2,2) c(2,4)
                c(4,1) c(4,2) c(4,4)];
            ctan = reshape(c([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);

        end
        %% Element K
        function Kel = Compute_Kel(~, Kel, gp, ~)
            % Definitions
            B=gp.B;
            D=gp.D;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            Kel = Kel + (B'*D*B) *gp.J *gp.w;
        end
        %% Element Fint
        function Fint = Compute_Fint(~, Fint, gp)
            Fint = Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
        end
    end
end

