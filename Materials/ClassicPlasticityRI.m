classdef ClassicPlasticityRI
    %Elastic 3D elastic class
    %   Detailed explanation goes here
    
    properties (Hidden, SetAccess = private)
        ndm;
        ndof;
        dNdX;
        finiteDisp = 0;
        
        I2 = eye(3);
        I4_dev
        I4_bulk
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
        
        eEff
        de
        dep
        theta
        bkStrss
        strss
    end
    %%
    methods
        %% Construct
        function obj = ClassicPlasticityRI(num, props, identity)
            obj.ndm  = num.ndm;
            obj.ndof = num.ndof;
            obj.eEff    = zeros(     num.gp, num.el, num.steps+1);
            obj.bkStrss = zeros(3,3, num.gp, num.el, num.steps+1);
            obj.strss   = zeros(3,3, num.gp, num.el, num.steps+1);
            
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
            
            obj.I4_dev= identity.I4_dev;
            obj.I4_bulk= identity.I4_bulk;
            
        end
        %% Epsilon
        function [eps, obj] = computeStrain(obj, gp, el, ~)
            eps = gp.B * el.Uvc;
            de=[1.0*eps(1) 0.5*eps(4) 0.5*eps(6)
                0.5*eps(4) 1.0*eps(2) 0.5*eps(5)
                0.5*eps(6) 0.5*eps(5) 1.0*eps(3)];
            de = de - (1/3)*trace(de)*eye(3);
            obj.de = de;
        end
        %% Sigma
        function [sigma_voigt, obj] = computeCauchy(obj, gp, step)
            I = obj.I2;
            G = obj.Shear;
            K = obj.Bulk;
            
            n.strss     = obj.strss(:,:, gp.i, gp.iel, step);
            np1.de      = obj.de;
            np1.dep     = obj.dep;
           
            sigma      = n.strss + 2*G* (np1.de - np1.dep) + K*sum(gp.eps(1:obj.ndm))*I;
            
            obj.strss(:,:, gp.i, gp.iel, step+1) = sigma;
            
            sigma_voigt = obj.voigtize(sigma,'col');
        end
        %% Tangential stiffness
        function [D, ctan, obj] = computeTangentStiffness(obj, gp, step)
            % Identities
            I    = obj.I2;
            I4_dev = obj.I4_dev;
            I4_bulk= obj.I4_bulk;
            % Material properties
            G    = obj.Shear;
            kappa= obj.Bulk;
            Kp   = obj.iso_modulus;
            Hp   = obj.kin_modulus;
            Y    = obj.yield_strss;
            % Values from previous step
            n.bkStrss = obj.bkStrss(:,:,  gp.i, gp.iel, step);
            n.strss   = obj.strss  (:,:,  gp.i, gp.iel, step);
            n.eEff    = obj.eEff   (gp.i, gp.iel, step);
            
            np1.de  = obj.de;
            np1.S   = n.strss - (1/3)*trace(n.strss)*I + 2*G*np1.de;
            
            np1.S_ht= np1.S - n.bkStrss;
            np1.sEqv= sqrt( (3/2)*sum(dot(np1.S_ht, np1.S_ht)) );
            
            
            np1.deEff    = obj.getEffPlasStrnInc(n.eEff, G, Y, Kp, Hp, np1.sEqv);
            obj.eEff(gp.i, gp.iel, step+1) =  n.eEff + np1.deEff;
            
            if (np1.sEqv*np1.deEff>0)
                np1.dir      = np1.S_ht/np1.sEqv;
                
                np1.dbkStrss = np1.deEff * Hp * np1.dir;
                np1.dep      = (3/2)*np1.deEff* np1.dir;
                
                np1.theta = 1 - 3*G* (np1.deEff/np1.sEqv);
                c2        = (1-np1.theta) - 1/(1+ (Kp+Hp)/(3*G)) ;
                c3        = (3/2)*c2/(np1.sEqv)^2;
            else
                np1.dbkStrss = zeros(size(np1.S_ht));
                np1.dep      = 0;
                
                np1.theta    = 1;
                c3 = 0;
            end
            % Save states
            obj.bkStrss(:,:, gp.i, gp.iel, step+1) = n.bkStrss + np1.dbkStrss;
            obj.dep = np1.dep;
            
            % Compute D
            S = obj.voigtize(np1.S_ht, 'row');
            
            D = kappa*I4_bulk            +... % Elastic bulk
                np1.theta*(2*G)*I4_dev   +... % Elastic deviatoric
                c3*(2*G)*(S'*S);              % Plastic
            
            ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
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
        
        function deEff = getEffPlasStrnInc(eEff_n, G, Y, Kp,Hp, sEquiv)
            deEff = (sEquiv - Y - eEff_n*Kp) / (Kp + Hp + 3*G);
            if deEff < 0
                deEff = 0;
            end
        end
        
        function vec = voigtize(mat, orientation)
            switch orientation
                case {'Col', 'col', 'Column', 'column'}
                    vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)]';
                case {'Row', 'row'}
                    vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)];
                otherwise
                    error('Unknown input');
            end
        end
    end
end

