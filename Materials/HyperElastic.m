classdef HyperElastic
   %FiniteElastic Computes Cauchy stress based on NeoHookean model
   properties (SetAccess = private)
      ndm;
      ndof;
      ngp;
      finiteDisp = 1;
      
      lam;
      mu;
      
      I
      I4_dev
      I4_bulk
   end
   %%
   methods
      %% Construct
      function ob = HyperElastic(num, props, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.ngp  = num.gp;
         
         [ob.lam, ob.mu] = ob.getProps(props);
         
         ob.I      = eye(ob.ndm);
         ob.I4_dev = identity.I4_dev;
         ob.I4_bulk= identity.I4_bulk;
      end
      %% Epsilon
      function [eps, ob] = Strain(ob, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, ~, ~)
         matE = diag([2,2,2,1,1,1]);
         JxX = det(gp.F);
         
         D = ob.mu*matE + JxX*ob.lam*( (2*JxX-1)*ob.I4_bulk - (JxX-1)*matE );
         
         sigma   = ob.mu*(gp.b - ob.I) + JxX*ob.lam*(JxX-1)*ob.I;
         sigma_v = T2T1(sigma,1);
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         B=gp.Bf;
         
         if ob.ndm == 2
            gp.D = gp.D([1,2,4],[1,2,4]);
         end
         D = formCombD(gp.sigma, gp.D, gp.finiteDisp);
         
         Kel = el.K + gp.J*gp.w* (B'*D*B);
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, ~)
         Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
      end
   end
   %% Static Methods
   methods (Static)
      function [lam, mu] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E  = props{j,2};
               case 'nu'
                  nu = props{j,2};
            end
         end
         lam = nu*E/((1+nu)*(1-2*nu));
         mu  = E/(2*(1+nu));
      end
   end
end