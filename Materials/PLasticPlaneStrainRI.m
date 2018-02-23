classdef PLasticPlaneStrainRI
    %Elastic 3D elastic class
    %   Detailed explanation goes here
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
        finiteDisp = 0;
        epsp;
        threshold;
        alpha;
        q;
        sigma_dev;
        
        Young
        Poisson
        mixedHardParam
        initYieldStress
        
        iso_modulus
        kin_modulus
        yield_strss
        e0
        edot
        
        dt
        de_p
        
        
        I = eye(3);
        
        time
        hard
        linear = false;
    end
    %%
    methods
        %% Construct
        function obj = PLasticPlaneStrainRI(num, props, time, hard)
            obj.ndm  = num.ndm;
            obj.ndof = num.ndof;
            obj.epsp = zeros(num.str,1);
            obj.time = time;
            obj.hard = hard;
            for i=1:length(props)
                switch props{i,1}
                    case 'E'
                        obj.Young          = props{i,2};
                    case 'v'
                        obj.Poisson        = props{i,2};
                    case 'M'
                        obj.mixedHardParam = props{i,2};
                    case 'sig0'
                        obj.initYieldStress= props{i,2};
                    otherwise
                        error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
                end
            end
            
        end
        %% Epsilon
        function [eps, obj] = computeStrain(obj, gp, el, ~)
            eps = gp.B * el.Uvc;
        end
        %% Sigma
        function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
            sigma_voigt = gp.D([1,2,4],[1,2,4]) *gp.eps;
%             sigma_voigt_trial = gp.D *(gp.eps - obj.epsp);
%             
%             if (sigma_voigt_trial <= obj.threshold)
%                 sigma_voigt = sigma_voigt_trial; 
%             else
%                 
%             end
        end
        %% Tangential stiffness
        function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~)
            E = obj.Young;
            v = obj.Poisson;
            
            Eh= E/(1-2*v)/(1+v);
            G = 0.5*E/(1+v);
            D =[Eh*(1-v)    Eh*v        Eh*v        0 0 0
                Eh*v        Eh*(1-v)    Eh*v        0 0 0
                Eh*v        Eh*v        Eh*(1-v)    0 0 0
                0           0           0           G 0 0
                0           0           0           0 G 0
                0           0           0           0 0 G];
            
            ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
            
        end
        %% Element K
        function Kel = computeK_el(~, Kel, gp, ~)
            % Definitions
            B=gp.B;
            D=gp.D([1,2,4],[1,2,4]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Kel = Kel + (B'*D*B) *gp.J *gp.w;
        end
        %% Element Fint
        function Fint = computeFint(~, gp, el)
            Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
        end
        %% Yielding function
        function value = get.threshold(obj)
            K     = obj.iso_modulus;
            H     = obj.kin_modulus;
            sig_y = obj.yield_strss;
            
            value = sig_y + K*alpha;
            
        end
        %% Compute plastic strain increment
        function value = get.de_p(obj)
            value=1;
            
        end
    end
end

