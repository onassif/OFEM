classdef ClassicPlasticityRI
    %Elastic 3D elastic class
    %   Detailed explanation goes here
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
        finiteDisp = 0;
        epsp;
        threshold;
        q;
        sigma_dev;
        
        Young
        Poisson
        iso_modulus
        kin_modulus
        yield_strss
        Shear
        Bulk
        e0
        edot
        
        
        dt
        de_p
        
        
        I2 = eye(3);
        IFour
        
        e_p
        Beta
        Alpha
    end
    %%
    methods
        %% Construct
        function obj = ClassicPlasticityRI(num, props, dt)
            obj.ndm  = num.ndm;
            obj.ndof = num.ndof;
            obj.epsp = zeros(num.str,1);
            obj.e_p  = zeros(3,3, num.gp,num.el);
            obj.Beta = zeros(3,3, num.gp,num.el);
            obj.Alpha= zeros(     num.gp,num.el);
            obj.dt   = dt;
            
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
                    case 'e0'
                        obj.e0          = props{i,2};
                    case 'edot'
                        obj.edot        = props{i,2};
                    otherwise
                        error("You've chosen Rate Independent Plastic material but specified incompatible material properties, I'm disapponted");
                end
            end
            obj.Shear =  0.5*obj.Young/(1+  obj.Poisson);
            obj.Bulk  =(1/3)*obj.Young/(1-2*obj.Poisson);
            
            obj.IFour = zeros(3,3,3,3);
            for i=1:3
                for j=1:3
                    obj.IFour(i,j,i,j)=1;
                end
            end
            
            
        end
        %% Sigma
        function sigma_voigt = Compute_cauchy(~, gp)
            sigma_voigt = gp.D *gp.eps;
%             sigma_voigt = gp.D([1,2,4],[1,2,4]) *gp.eps;
%             %             sigma_voigt_trial = gp.D *(gp.eps - obj.epsp);
%             %
%             %             if (sigma_voigt_trial <= obj.threshold)
%             %                 sigma_voigt = sigma_voigt_trial;
%             %             else
%             %
%             %             end
        end
        %% Tangential stiffness
        function [D, ctan, obj] = Compute_tangentstiffness(obj, gp)
            E = obj.Young;
            v = obj.Poisson;
            G    = obj.Shear;
            kappa= obj.Bulk;
            I    = obj.I2;
            I4   = obj.IFour;
            eps  = [...
                1.0*gp.eps(1) 0.5*gp.eps(4) 0.5*gp.eps(6)
                0.5*gp.eps(4) 1.0*gp.eps(2) 0.5*gp.eps(5)
                0.5*gp.eps(6) 0.5*gp.eps(5) 1.0*gp.eps(3)];
            
            ep     = obj.e_p  (:,:,gp.i, gp.iel);
            beta   = obj.Beta (:,:,gp.i, gp.iel);
            alpha_n= obj.Alpha(    gp.i, gp.iel);
            
            e = eps - 1/3.*trace(eps)*I;
            [K_n,~,H_n,~] = obj.hardening(alpha_n);
            
            sTrial      = 2*G*(e - ep);
            xiTrial     = sTrial - beta;
            normXiTrial = norm(xiTrial);
            ftrial      = normXiTrial - sqrt(2/3) * K_n;
            
            if ftrial<=0
%                 decision.s = sTrial;
%                 decision.xi=xiTrial;

            Eh= E/(1-2*v)/(1+v);
            
            D =[Eh*(1-v)    Eh*v        Eh*v        0 0 0
                Eh*v        Eh*(1-v)    Eh*v        0 0 0
                Eh*v        Eh*v        Eh*(1-v)    0 0 0
                0           0           0           G 0 0
                0           0           0           0 G 0
                0           0           0           0 0 G];
            
            ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);

            else
                n = xiTrial ./ normXiTrial;
                [deltaGamma, alpha_np1] = obj.box3d1(alpha_n, normXiTrial, G);
                [~,Kp_np1,H_np1,Hp_np1] = obj.hardening(alpha_np1);
                
                
                obj.e_p  (:,:,gp.i, gp.iel)= ep + deltaGamma.*n;
                obj.Beta (:,:,gp.i, gp.iel)= beta + sqrt(2/3).*(H_np1 - H_n).*n;
                obj.Alpha(    gp.i, gp.iel)= alpha_np1;
%                 ep   = ep + deltaGamma.*n;
%                 beta = beta + sqrt(2/3).*(H_np1 - H_n).*n;
%                 sigma=kappa*trace(eps)*I + sTrial - 2*G*deltaGamma*n;
                
                theta   =1-2*G*deltaGamma/normXiTrial;
                thetaBar=1/(1+(Kp_np1+Hp_np1)/(3*G)) - (1-theta);
                
                ctan = zeros(3,3,3,3);
                for i=1:3
                    for j=1:3
                        for k=1:3
                            for l=1:3
                                ctan(i,j,k,l) = ...
                                    ctan(i,j,k,l) + ...
                                    kappa*I(i,j)*I(k,l) + ...
                                    2*G*theta* (I4(i,j,k,l) - (1/3)*I(i,j)*I(k,l)) - ...
                                    2*G*thetaBar*n(i,j)*n(k,l);
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
        function Kel = Compute_Kel(~, Kel, gp, ~)
            % Definitions
            B=gp.B;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Kel = Kel + (B'*gp.D*B) *gp.J *gp.w;
        end
        %% Element Fint
        function Fint = Compute_Fint(~, Fint, gp)
            Fint = Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
        end
        %% Yielding function
        function value = get.threshold(obj)
            K     = obj.iso_modulus;
            H     = obj.kin_modulus;
            sig_y = obj.yield_strss;
            
            value = sig_y + K*alpha;
            
        end
        %% Compute plastic strain increment
        function value = get.de_p(~)
            value=1;
            
        end
        
%         function decision = box3d2(obj, gp)
%             E = obj.Young;
%             v = obj.Poisson;
%             G    = obj.Shear;
%             kappa= obj.Bulk;
%             I    = obj.I2;
%             I4   = obj.obj.IFour;
%             eps  = [...
%                 1.0*gp.eps(1) 0.5*gp.eps(4) 0.5*gp.eps(6)
%                 0.5*gp.eps(4) 1.0*gp.eps(2) 0.5*gp.eps(5)
%                 0.5*gp.eps(6) 0.5*gp.eps(5) 1.0*gp.eps(3)];
%             
%             ep     = obj.e_p  (:,:,gp.i, gp.iel);
%             beta   = obj.Beta (:,:,gp.i, gp.iel);
%             alpha_n= obj.Alpha(    gp.i, gp.iel);
%             
%             e = eps - 1/3.*trace(eps)*I;
%             [K_n,~,H_n,~] = hardening(alpha_n);
%             
%             sTrial      = 2*G*(e - ep);
%             xiTrial     = sTrial - beta;
%             normXiTrial = norm(xiTrial);
%             ftrial      = normXiTrial - sqrt(2/3) * K_n;
%             
%             if ftrial<=0
% %                 decision.s = sTrial;
% %                 decision.xi=xiTrial;
% 
%             Eh= E/(1-2*v)/(1+v);
%             
%             D =[Eh*(1-v)    Eh*v        Eh*v        0 0 0
%                 Eh*v        Eh*(1-v)    Eh*v        0 0 0
%                 Eh*v        Eh*v        Eh*(1-v)    0 0 0
%                 0           0           0           G 0 0
%                 0           0           0           0 G 0
%                 0           0           0           0 0 G];
%             
%             ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
% 
%             else
%                 n = xiTrial ./ normXiTrial;
%                 [deltaGamma, alpha_np1] = box3d1(alpha_n, normXiTrial, G);
%                 [~,Kp_np1,H_np1,Hp_np1] = hardening(alpha_np1);
%                 
%                 beta = beta + sqrt(2/3).*(H_np1 - H_n).*n;
%                 ep   = ep + deltaGamma.*n;
%                 sigma=kappa*trace(eps)*I + sTrial - 2*G*deltaGamma*n;
%                 
%                 theta   =1-2*G*deltaGamma/normXiTrial;
%                 thetaBar=1/(1+(Kp_np1+Hp_np1)/(3*G)) - (1-theta);
%                 
%                 ctan = zeros(3,3,3,3);
%                 for i=1:3
%                     for j=1:3
%                         for k=1:3
%                             for l=1:3
%                                 ctan(i,j,k,l) = ...
%                                     ctan(i,j,k,l) + ...
%                                     kappa*I(i,j)*I(k,l) + ...
%                                     2*G*theta* (I4(i,j,k,l) - (1/3)*I(i,j)*I(k,l)) - ...
%                                     2*G*thetaBar*n(i,j)*n(k,l);
%                             end
%                         end
%                     end
%                 end
%                 D =[...
%                     ctan(1,1,1,1) ctan(1,1,2,2) ctan(1,1,3,3) ctan(1,1,1,2) ctan(1,1,2,3) ctan(1,1,1,3)
%                     ctan(2,2,1,1) ctan(2,2,2,2) ctan(2,2,3,3) ctan(2,2,1,2) ctan(2,2,2,3) ctan(2,2,1,3)
%                     ctan(3,3,1,1) ctan(3,3,2,2) ctan(3,3,3,3) ctan(3,3,1,2) ctan(3,3,2,3) ctan(3,3,1,3)
%                     ctan(1,2,1,1) ctan(1,2,2,2) ctan(1,2,3,3) ctan(1,2,1,2) ctan(1,2,2,3) ctan(1,2,1,3)
%                     ctan(2,3,1,1) ctan(2,3,2,2) ctan(2,3,3,3) ctan(2,3,1,2) ctan(2,3,2,3) ctan(2,3,1,3)
%                     ctan(1,3,1,1) ctan(1,3,2,2) ctan(1,3,3,3) ctan(1,3,1,2) ctan(1,3,2,3) ctan(1,3,1,3)];
%             end
%         end
    end
    
    methods (Static)
        
        function [deltaGamma, alpha_np1] = box3d1(alpha, normXiTrial, G)
            tol = 1e-9;
            g = inf;
            
            deltaGamma = 0;
            alpha_k = alpha;
            [~,~,H_n,~] = hardening(alpha);
            
            while(abs(g)<tol)
                [K,Kp,H,Hp] = hardening(alpha_k);
                
                g = -sqrt(2/3)*K + normXiTrial...
                    -( 2*G*deltaGamma + sqrt(2/3)*(H - H_n) );
                
                Dg= -2*G*( 1 + (Hp+Kp)/(3*G) );
                
                deltaGamma = deltaGamma - g/Dg;
                alpha_k      = alpha + sqrt(2/3)*deltaGamma;
            end
            alpha_np1 = alpha_k;
        end
        
        function [K,Kp,H,Hp] = hardening(alpha)
            K = 60 + 10*alpha;
            Kp= 60;
            H = 40*alpha;
            Hp= 40;
        end
    end
end

