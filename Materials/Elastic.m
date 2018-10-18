classdef Elastic
   %Elastic 3D class
   properties (SetAccess = private)
      ndm;
      ndof;
      ngp;
      finiteDisp = 0;
      linear = true;
      C0;
   end
   %%
   methods
      %% Construct
      function ob = Elastic(num, props, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.ngp  = num.gp;
         
         [E, nu] = ob.getProps(props);
         
         G =   0.5*E/(1+  nu);
         K = (1/3)*E/(1-2*nu);
         
         ob.C0 = K*identity.I4_bulk + 2*G*identity.I4_dev;
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma and Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, ~, ~)
         D = ob.C0;
         if ob.ndm == 2
            sigma_v = D([1,2,4],[1,2,4])*gp.eps;
         elseif ob.ndm == 3
            sigma_v = D*gp.eps;
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         if (ob.ndm == 2); gp.D = gp.D([1,2,4],[1,2,4]); end
         
         Kel = el.K + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, ~)
         Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
      end
   end
   methods (Static)
      function [E, nu] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E = props{j,2};
               case 'nu'
                  nu = props{j,2};
            end
         end
      end
   end
end