classdef ClassicPlasticityRI
    %Elastic 3D elastic class
    %   Detailed explanation goes here
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
        finiteDisp = 0;
        
        I2 = eye(3);
        IFour
        linear = false;
    end
    
    properties (SetAccess = private)
        
        Young
        Poisson
        iso_modulus
        kin_modulus
        yield_strss
        Shear
        Bulk
        
        ep
        beta
        alpha
        
        plastic = false;
        sTrial
        dGamma
        dir;
    end
    %%
    methods
        %% Construct
        function obj = ClassicPlasticityRI(num, props, identity)
            obj.ndm  = num.ndm;
            obj.ndof = num.ndof;
            obj.ep   = zeros(3,3, num.gp, num.el, num.steps+1);
            obj.beta = zeros(3,3, num.gp, num.el, num.steps+1);
            obj.alpha= zeros(     num.gp, num.el, num.steps+1);
            
            for i=1:length(props)
                switch props{i,1}
                    case 'E'
                        obj.Young       = props{i,2};
                    case 'v'
                        obj.Poisson     = props{i,2};
                    case 'K'
                        obj.iso_modulus = props{i,2};
                    case 'H'
                        obj.kin_modulus = props{i,2};
                    case 'Y'
                        obj.yield_strss = props{i,2};
                    otherwise
                        error("You've chosen Rate Independent Plastic material but specified incompatible material properties, I'm disapponted");
                end
            end
            obj.Shear =  0.5*obj.Young/(1+  obj.Poisson);
            obj.Bulk  =(1/3)*obj.Young/(1-2*obj.Poisson);
            
            obj.IFour = identity;
            
            
        end
        %% Epsilon
        function [eps, obj] = computeStrain(obj, gp, el, ~)
            eps = gp.B * el.Uvc;
        end
        %% Sigma
        function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
            if obj.plastic
                G    = obj.Shear;
                kappa= obj.Bulk;
                sigma= kappa*sum(gp.eps(1:obj.ndm))*obj.I2 + obj.sTrial - 2*G*obj.dGamma*obj.dir;
                sigma_voigt = [...
                    sigma(1,1)
                    sigma(2,2)
                    sigma(3,3)
                    sigma(1,2)
                    sigma(2,3)
                    sigma(1,3)];
            else
                sigma_voigt = gp.D *gp.eps;
            end
        end
        %% Tangential stiffness
        function [D, ctan, obj] = computeTangentStiffness(obj, gp, step)
            % Identities
            I    = obj.I2;
            I4   = obj.IFour;
            % Material properties
            E    = obj.Young;
            v    = obj.Poisson;
            G    = obj.Shear;
            kappa= obj.Bulk;
            Kp   = obj.iso_modulus;
            Hp   = obj.kin_modulus;
            Y    = obj.yield_strss;
            % Values from previous step
            n    = struct(...
                'ep'   , obj.ep   (:,:, gp.i, gp.iel, step), ...
                'beta' , obj.beta (:,:, gp.i, gp.iel, step), ...
                'alpha', obj.alpha(     gp.i, gp.iel, step));
            % Strain full matrix
            np1.eps  = [...
                1.0*gp.eps(1) 0.5*gp.eps(4) 0.5*gp.eps(6)
                0.5*gp.eps(4) 1.0*gp.eps(2) 0.5*gp.eps(5)
                0.5*gp.eps(6) 0.5*gp.eps(5) 1.0*gp.eps(3)];
            
            np1.e = np1.eps - 1/3.*trace(np1.eps)*I;
            %             [K_n,~,~,~] = obj.hardening(alpha_n);
            n.K = ( Y + Kp*n.alpha );
            
            np1.sTrial      = 2*G*(np1.e - n.ep);
            np1.xiTrial     = np1.sTrial - n.beta;
            np1.normXiTrial = norm(eig(np1.xiTrial));
            np1.fTrial      = np1.normXiTrial - sqrt(2/3) * n.K;
            
            if np1.fTrial<=0
                obj.plastic = false;
                
                Eh= E/(1-2*v)/(1+v);
                
                D =[Eh*(1-v)    Eh*v        Eh*v        0 0 0
                    Eh*v        Eh*(1-v)    Eh*v        0 0 0
                    Eh*v        Eh*v        Eh*(1-v)    0 0 0
                    0           0           0           G 0 0
                    0           0           0           0 G 0
                    0           0           0           0 0 G];
                
                ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
                
                obj.ep   (:,:, gp.i, gp.iel, step+1)= obj.ep   (:,:, gp.i, gp.iel, step);
                obj.beta (:,:, gp.i, gp.iel, step+1)= obj.beta (:,:, gp.i, gp.iel, step);
                obj.alpha(     gp.i, gp.iel, step+1)= obj.alpha(     gp.i, gp.iel, step);
                
            else
                obj.plastic = true;
                
                np1.N      = np1.xiTrial ./ np1.normXiTrial;
                np1.dGamma = np1.fTrial / (2*G + (2/3)*(Kp+Hp));
                obj.dir    = np1.N;
                obj.dGamma = np1.dGamma;
                obj.sTrial = np1.sTrial;
                
                obj.ep   (:,:, gp.i, gp.iel, step+1)= n.ep    + np1.dGamma.*np1.N;
                obj.beta (:,:, gp.i, gp.iel, step+1)= n.beta  + np1.dGamma.*np1.N*(2/3)*Hp;
                obj.alpha(     gp.i, gp.iel, step+1)= n.alpha + sqrt(2/3)*np1.dGamma;
                
                np1.theta   = 1 - (2*G*np1.dGamma)/np1.normXiTrial;
                %                 thetaBar= deltaGamma + theta - 1;
                np1.thetaBar= 1 / (1 + (Kp+Hp)/(3*G)) + np1.theta - 1;
                
                ctan = zeros(3,3,3,3);
                for i=1:3
                    for j=1:3
                        for k=1:3
                            for l=1:3
                                ctan(i,j,k,l) = ...
                                    kappa*I(i,j)*I(k,l) + ...
                                    2*G*np1.theta*( I4(i,j,k,l) - (1/3)*I(i,j)*I(k,l) ) - ...
                                    2*G*np1.thetaBar * ( np1.N(i,j)*np1.N(k,l) );
                            end
                        end
                    end
                end
                D =[...
                    ctan(1,1,1,1) ctan(1,1,2,2) ctan(1,1,3,3) ctan(1,1,1,2) ctan(1,1,2,3) ctan(1,1,1,3)
                    ctan(2,2,1,1) ctan(2,2,2,2) ctan(2,2,3,3) ctan(2,2,1,2) ctan(2,2,2,3) ctan(2,2,1,3)
                    ctan(3,3,1,1) ctan(3,3,2,2) ctan(3,3,3,3) ctan(3,3,1,2) ctan(3,3,2,3) ctan(3,3,1,3)
                    ctan(1,2,1,1) ctan(1,2,2,2) ctan(1,2,3,3) ctan(1,2,1,2) ctan(1,2,2,3) ctan(1,2,1,3)
                    ctan(2,3,1,1) ctan(2,3,2,2) ctan(2,3,3,3) ctan(2,3,1,2) ctan(2,3,2,3) ctan(2,3,1,3)
                    ctan(1,3,1,1) ctan(1,3,2,2) ctan(1,3,3,3) ctan(1,3,1,2) ctan(1,3,2,3) ctan(1,3,1,3)];
            end
        end
        %% Element K
        function Kel = computeK_el(~, Kel, gp, ~)
            % Definitions
            B=gp.B;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Kel = Kel + (B'*gp.D*B) *gp.J *gp.w;
        end
        %% Element Fint
        function Fint = computeFint(~, gp, el)
            Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
        end
    end
    
    methods (Static)
        
        %         function [deltaGamma, alpha_np1] = box3d1(alpha, normXiTrial, G)
        %             tol = 1e-9;
        %             g = inf;
        %
        %             deltaGamma = 0;
        %             alpha_k = alpha;
        %             [~,~,H_n,~] = ClassicPlasticityRI.hardening(alpha);
        %
        %             while(abs(g)>tol)
        %                 [K,Kp,H,Hp] = ClassicPlasticityRI.hardening(alpha_k);
        %
        %                 g = -sqrt(2/3)*K + normXiTrial...
        %                     -( 2*G*deltaGamma + sqrt(2/3)*(H - H_n) );
        %
        %                 Dg= -2*G*( 1 + (Hp+Kp)/(3*G) );
        %
        %                 deltaGamma = deltaGamma - g/Dg;
        %                 alpha_k      = alpha + sqrt(2/3)*deltaGamma;
        %             end
        %             alpha_np1 = alpha_k;
        %         end
        
        function [K,Kp,H,Hp] = hardening(alpha)
            K = (60 + 10*alpha);
            Kp= 10;
            H = 40*alpha;
            Hp= 40;
        end
    end
end

