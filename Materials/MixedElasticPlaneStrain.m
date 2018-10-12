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
      function obj = MixedElasticPlaneStrain(num, props, identity)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.nen  = num.nen;
         obj.ngp  = num.gp;
         obj.diff = obj.ndof - obj.ndm;
         switch num.nen
            case {3,4} % T3 and Q4
               obj.M.N = 1;
            case 6 % T6
               xi = 1/6*[...
                  4 1 1
                  1 1 4]';
               obj.M = T3(0, 0, xi);
            case 9 % Q9
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1  0 +1  0 -1 0
                  -1 -1 +1 +1 -1  0 +1  0 0]';
               obj.M = Q4(0, 0, xi);
         end
         
         [E,v] = obj.getProps(props);
         G =   0.5*E/(1+  v);
         K = (1/3)*E/(1-2*v);
         
         obj.C0 = K*identity.I4_bulk + 2*G*identity.I4_dev;
         
         obj.Bulk = K;
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~, ~)
         sigma_voigt = gp.D *gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~, ~)
         c = obj.C0;
         
         D = c([1,2,4],[1,2,4]);
         ctan = reshape(c([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
      end
      %% Element K
      function kel = computeK_el(obj, gp, el, ~)
         dN      = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         obj.M.i = gp.i;
         K = obj.Bulk;
         if gp.i == 1
            kel = [...
               (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN*obj.M.N
               obj.M.N'*dN'                          (1/K)*obj.M.N'*obj.M.N] .*gp.j*gp.w;
         else
            kel = el.K + [...
               (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN*obj.M.N
               obj.M.N'*dN'                          (1/K)*obj.M.N'*obj.M.N] .*gp.j*gp.w;
         end
      end
      %% Element Fint
      function Fint = computeFint(obj, gp, el, ~)
         numU = numel(el.nodes);
         dN = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         K = obj.Bulk;
         
         Sdev = [gp.sigma(1:2) - 1/3*sum(obj.C0(1:3,1:3)*[gp.eps(1:2);0]); gp.sigma(3)];
         p = el.U_global(numU+el.i);
         if gp.i == 1
            Fint = [gp.B'*Sdev; 0] *gp.j *gp.w;
         else
            Fint = el.Fint + [gp.B'*Sdev; 0] *gp.j *gp.w;
         end
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