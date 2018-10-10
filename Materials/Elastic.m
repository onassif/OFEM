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
      %% Sigma
      function [sigma_voigt, ob] = computeCauchy(ob, gp, ~, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, ob] = computeTangentStiffness(ob, ~, ~, ~)
         D = ob.C0;
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if ob.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(~, gp, el, ~)
         if gp.i == 1
            Kel = (gp.B'*gp.D*gp.B) *gp.J *gp.w;
         else
            Kel = el.K + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
         end
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, ~)
         if gp.i == 1
            Fint = (gp.B'*gp.sigma) *gp.J *gp.w;
         else
            Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
         end
      end
   end
   methods (Static)
      function [E, v] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E = props{j,2};
               case 'nu'
                  v = props{j,2};
            end
         end
      end
   end
end