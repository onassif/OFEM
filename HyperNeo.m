classdef HyperNeo
    %HyperNeo Computes Cauchy stress based on NeoHookean model
    %   Currently only if ndm=2, ndof=2
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        finiteDisp = 1;
    end
    %%
    methods
        function obj = HyperNeo(ndm, ndof)
            obj.ndm  = ndm;
            obj.ndof = ndof;
        end
        %% Sigma
        function sigma_voigt = Compute_cauchy(obj, gp, props)
            
            if strcmp(props{1,1} ,'mu') && strcmp(props{2,1} ,'lambda')
                mu      = props{1,2};
                lambda  = props{2,2};
            elseif strcmp(props{1,1} ,'lambda') && strcmp(props{2,1} ,'mu')
                mu      = props{1,2};
                lambda  = props{2,2};
            else
                error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
            end
            
            if (obj.ndm==2 && obj.ndof ==2)
                I = eye(obj.ndm);
                b = gp.b;
                J = gp.det_F;
                
                sigma  = mu/J* (b-I) + lambda*(J-1)*I;
                sigma_voigt = [sigma(1,1); sigma(2,2); sigma(1,2)];
            end
        end
        %% Tangential stiffness
        function [D, ctan] = Compute_tangentstiffness(~, gp, props)
            
            if strcmp(props{1,1} ,'mu') && strcmp(props{2,1} ,'lambda')
                mu      = props{1,2};
                lambda  = props{2,2};
            elseif strcmp(props{1,1} ,'lambda') && strcmp(props{2,1} ,'mu')
                mu      = props{1,2};
                lambda  = props{2,2};
            else
                error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
            end
            
            c = zeros(3,3,3,3);
            I = eye(3);
            J = gp.det_F;
            
            for i=1:3
                for j=1:3
                    for k=1:3
                        for l=1:3
                            c(i,j,k,l) = mu/J*( I(i,k)*I(j,l) + I(i,l)*I(j,k) )+...
                                lambda*(  (2*J-1)*I(i,j)*I(k,l) - (J-1)*( I(i,k)*I(j,l) + I(i,l)*I(j,k) )  );
                        end
                    end
                end
            end
            
            D = [ c(1,1,1,1) c(1,1,2,2) c(1,1,1,2)
                c(2,2,1,1) c(2,2,2,2) c(2,2,1,2)
                c(1,2,1,1) c(1,2,2,2) c(1,2,1,2)];
            ctan = c;
        end
        %% Element K
        function Kel = Compute_Kel(obj, Kel, gp, ngp)
            % Definitions
            B=gp.B;
            D=gp.D;
            ndf = obj.ndof;
            sigma=[gp.sigma(1) gp.sigma(3)
                gp.sigma(3) gp.sigma(2)];
            dNdx=gp.dNdx;
            K_geo = zeros(ndf*ngp,ndf*ngp);
            I = eye(ndf);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for a =1:ngp
                for b=1:ngp
                    for k=1:ndf
                        for i=1:ndf
                            for l=1:ndf
                                for j=1:ndf
                                    aa=(a-1)*ndf; bb=(b-1)*ndf;
                                    
                                    K_geo(aa+i, bb+k) = K_geo(aa+i, bb+k) +...
                                        I(l,k)*dNdx(a,l) *sigma(i,j)*dNdx(b,j);
                                end
                            end
                        end
                    end
                end
            end
            Kel = Kel + ( (B'*D*B) + K_geo ) *gp.j*gp.w;
        end
        %% Element Fint
        function Fint = Compute_Fint(~, Fint, gp)
            
            Fint = Fint + (gp.B'*gp.sigma) *gp.j *gp.w;
        end
    end
end

