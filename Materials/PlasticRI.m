classdef PlasticRI
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      dNdX;
      finiteDisp = 0;
      
      I
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
      G
      Bulk
      
      eEff
      e
      ep
      theta
      bkStrss
   end
   %%
   methods
      %% Construct
      function obj = PlasticRI(num, props, identity)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.eEff    = zeros(     num.gp, num.el, num.steps+1);
         obj.bkStrss = zeros(num.ndm,num.ndm, num.gp, num.el, num.steps+1);
         obj.ep      = zeros(num.ndm,num.ndm, num.gp, num.el, num.steps+1);
         
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
         obj.G    =  0.5*obj.Young/(1+  obj.Poisson);
         obj.Bulk =(1/3)*obj.Young/(1-2*obj.Poisson);
         
         obj.I      = eye(obj.ndm);
         obj.I4_dev = identity.I4_dev;
         obj.I4_bulk= identity.I4_bulk;
         
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         I = obj.I;
         
         eps= gp.B * el.Uvc;
         if obj.ndm == 2
            e=[1.0*eps(1) 0.5*eps(3)
               0.5*eps(3) 1.0*eps(2)];
         elseif obj.ndm == 3
            e=[1.0*eps(1) 0.5*eps(4) 0.5*eps(6)
               0.5*eps(4) 1.0*eps(2) 0.5*eps(5)
               0.5*eps(6) 0.5*eps(5) 1.0*eps(3)];
         end
         e = e - (1/3)*trace(e)*I;
         obj.e = e;
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, step)
         I = obj.I;
         G = obj.G;
         K = obj.Bulk;
         
         np1.e  = obj.e;
         np1.ep = obj.ep(:,:, gp.i, gp.iel, step+1);
         
         sigma       = 2*G* (np1.e - np1.ep) + K*sum(gp.eps(1:obj.ndm))*I;
         sigma_voigt = obj.voigtize(sigma,'col', obj.ndm);
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, step)
         
         % Identities
         I4_dev = obj.I4_dev;
         I4_bulk= obj.I4_bulk;
         % Material properties
         G    = obj.G;
         kappa= obj.Bulk;
         Kp   = obj.iso_modulus;
         Hp   = obj.kin_modulus;
         Y    = obj.yield_strss;
         % Values from previous step
         n.ep      = obj.ep     (:,:, gp.i, gp.iel, step);
         n.eEff    = obj.eEff   (     gp.i, gp.iel, step);
         n.bkStrss = obj.bkStrss(:,:, gp.i, gp.iel, step);
         
         % Trial equivalent stress
         np1.e    = obj.e;
         np1.S_ht = 2*G*(np1.e - n.ep) - n.bkStrss;
         np1.sEqv = sqrt( (3/2)*sum(dot(np1.S_ht, np1.S_ht)) );
         
         
         np1.deEff    = obj.getEffPlasStrnInc(n.eEff, G, Y, Kp, Hp, np1.sEqv);
         obj.eEff(gp.i, gp.iel, step+1) =  n.eEff + np1.deEff;
         
         if (np1.sEqv*np1.deEff>0)
            np1.dbkStrss = np1.deEff * Hp * (np1.S_ht/np1.sEqv);
            np1.dep      = (3/2)*np1.deEff* (np1.S_ht/np1.sEqv);
            
            np1.theta = 1 - 3*G* (np1.deEff/np1.sEqv);
            c2        = (1-np1.theta) - 1/(1+ (Kp+Hp)/(3*G)) ;
            c3        = (3/2)*c2/(np1.sEqv)^2;
         else
            np1.dbkStrss = zeros(size(np1.S_ht));
            np1.dep      = 0;
            
            np1.theta    = 1;
            c3 = 0;
         end
         np1.ep = n.ep + np1.dep;
         
         % Save states
         obj.bkStrss(:,:, gp.i, gp.iel, step+1) = n.bkStrss + np1.dbkStrss;
         obj.ep     (:,:, gp.i, gp.iel, step+1) = np1.ep;
         
         % Compute D
         if     obj.ndm == 2
            S=[np1.S_ht(1,1) np1.S_ht(2,2) 0 ...
               np1.S_ht(1,2) 0             0];
         elseif obj.ndm == 3
            S=[np1.S_ht(1,1) np1.S_ht(2,2) np1.S_ht(3,3)...
               np1.S_ht(1,2) np1.S_ht(2,3) np1.S_ht(1,3)];
         end
         
         D =kappa*I4_bulk            +... % Elastic bulk
            np1.theta*(2*G)*I4_dev   +... % Elastic deviatoric
            c3*(2*G)*(S'*S);              % Plastic
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         if obj.ndm == 2
            D =D([1,2,4],[1,2,4]);
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
      
      function deEff = getEffPlasStrnInc(eEff_n, G, Y, Kp,Hp, sEquiv)
         deEff = (sEquiv - Y - eEff_n*Kp) / (Kp + Hp + 3*G);
         if deEff < 0
            deEff = 0;
         end
      end
      
      function vec = voigtize(mat, orientation, ndm)
         switch orientation
            case {'Col', 'col', 'Column', 'column'}
               if     ndm == 2
                  vec = [mat(1,1) mat(2,2) mat(1,2)]';
               elseif ndm == 3
                  vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)]';
               end
            case {'Row', 'row'}
               if     ndm == 2
                  vec = [mat(1,1) mat(2,2) mat(1,2)];
               elseif ndm == 3
                  vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)];
               end
            otherwise
               error('Unknown input');
         end
      end
   end
end

