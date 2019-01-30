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
      Cdev
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
         
         ob.C0   = K*identity.I4_bulk + 2*G*identity.I4_dev;
         
         ob.Cdev = 2*G*identity.I4_dev;
         
         ob.Bulk = K;
      end
      %% Compute gp K, F and in case of dynamics: M
      function [gp, el, ob] = computeKFM(ob, gp, el, ~)
         gp.eps = gp.B * el.Uvc;

         gp.D = ob.Cdev;
         
         if (ob.ndm == 2)
            D = gp.D([1,2,4],[1,2,4]);
            Iv = [1 1 0];
         else
            D = gp.D;
            Iv = [1 1 1 0 0 0];
         end
         dN = Iv*gp.B;
         gp.sigma = D*gp.eps;
         
         ob.M.i = gp.i;
         K = ob.Bulk;

         el.K = el.K + [...
            gp.B'*D*gp.B    dN'*ob.M.N
            ob.M.N'*dN    -(1/K)*ob.M.N'*ob.M.N] .*gp.j*gp.w;

%       function Fint = computeFint(ob, gp, el, ~)
         numU = numel(el.nodes);

         p = el.U_global(numU+el.i);
         
         el.Fint = el.Fint + [...
            gp.B'*gp.sigma    + dN'*ob.M.N*p 
            ob.M.N'*Iv*gp.eps - (1/K)*ob.M.N'*ob.M.N*p] *gp.j *gp.w;
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