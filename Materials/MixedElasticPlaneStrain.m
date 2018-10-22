classdef MixedElasticPlaneStrain
   %PlaneStrain Computes Cauchy stress based on 2D plane strain
   %   Currently only if ndm=2, ndof=2
   properties (SetAccess = private)
      ndm;
      ndof;
      nen;
      ngp;
      dNdX;
      diff;
      finiteDisp = 0;
      
      C0;
      Bulk;
      
      npress = 1;
      M = struct('i',[],'N',[]);
   end
   %%
   methods
      %% Construct
      function ob = MixedElasticPlaneStrain(num, props, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.nen  = num.nen;
         ob.ngp  = num.gp;
         ob.diff = ob.ndof - ob.ndm;
         switch num.nen
            case {3,4} % T3 and Q4
               ob.M.N = 1;
            case 6 % T6
               xi = 1/6*[...
                  4 1 1
                  1 1 4]';
               ob.M = T3(0, 0, xi);
            case 9 % Q9
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1  0 +1  0 -1 0
                  -1 -1 +1 +1 -1  0 +1  0 0]';
               ob.M = Q4(0, 0, xi);
         end
         
         [E,v] = ob.getProps(props);
         G =   0.5*E/(1+  v);
         K = (1/3)*E/(1-2*v);
         
         ob.C0 = K*identity.I4_bulk + 2*G*identity.I4_dev;
         
         ob.Bulk = K;
      end
      %% Epsilon
      function [eps, ob] = Strain(ob, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = computeTangentStiffness(ob, gp, ~, ~)
         D       = ob.C0;
         sigma_v = ob.C0([1,2,4],[1,2,4])*gp.eps;
      end
      %% Element K
      function kel = computeK_el(ob, gp, el, ~)
         if (ob.ndm == 2); gp.D = gp.D([1,2,4],[1,2,4]); end
         dN      = reshape(gp.dNdx',ob.nen*ob.ndm,1);
         ob.M.i = gp.i;
         K = ob.Bulk;
         kel = el.K + [...
            (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN*ob.M.N
            ob.M.N'*dN'                          (1/K)*ob.M.N'*ob.M.N] .*gp.j*gp.w;
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, ~)
         numU = numel(el.nodes);
         dN = reshape(gp.dNdx',ob.nen*ob.ndm,1);
         K = ob.Bulk;
         
         Sdev = [gp.sigma(1:2) - 1/3*sum(ob.C0(1:3,1:3)*[gp.eps(1:2);0]); gp.sigma(3)];
         p = el.U_global(numU+el.i);
         
         Fint = el.Fint + [gp.B'*Sdev; 0] *gp.j *gp.w;
      end
   end
   methods (Static)
      function [E, v] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E = props{j,2};
               case 'v'
                  v = props{j,2};
            end
         end
      end
   end
end