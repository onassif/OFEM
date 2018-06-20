classdef MixedElasticPlaneStrain
   %PlaneStrain Computes Cauchy stress based on 2D plane strain
   %   Currently only if ndm=2, ndof=2
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      nen;
      ngp;
      dNdX;
      diff;
      finiteDisp = 0;
      
      Young;
      Poisson;
      Bulk;
      
      npress = 1;
      M = struct('i',[],'N',[]);
      linear = true;
      
      name = 'MixedElasticPlaneStrain';
   end
   %%
   methods
      %% Construct
      function obj = MixedElasticPlaneStrain(num, props)
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
         
         if strcmp(props{1,1} ,'E') && strcmp(props{2,1} ,'v')
            obj.Young  = props{1,2};
            obj.Poisson= props{2,2};
         elseif strcmp(props{1,1} ,'v') && strcmp(props{2,1} ,'E')
            obj.Poisson= props{1,2};
            obj.Young  = props{2,2};
         else
            error(['You''ve chosen Plane Strain material but specified',...
               'incompatible material properties, I''m disapponted']);
         end
         
         obj.Bulk = (1/3)*obj.Young/(1-2*obj.Poisson);
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
         sigma_voigt = gp.D *gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~, ~)
         E = obj.Young;
         v = obj.Poisson;
         
         Eh= E/(1-2*v)/(1+v);
         G = 0.5*E/(1+v);
         c =[Eh*(1-v)	Eh*v        Eh*v        0 0 0
            Eh*v        Eh*(1-v)    Eh*v        0 0 0
            Eh*v        Eh*v        Eh*(1-v)    0 0 0
            0           0           0           G 0 0
            0           0           0           0 G 0
            0           0           0           0 0 G];
         
         
         D=[c(1,1) c(1,2) c(1,4)
            c(2,1) c(2,2) c(2,4)
            c(4,1) c(4,2) c(4,4)];
         ctan = reshape(c([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
      end
      %% Element K
      function kel = computeK_el(obj, kel, gp, ~)
         dN      = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         obj.M.i = gp.i;
         K = obj.Bulk;
         if gp.i == 1
         kel = [...
            (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN*obj.M.N
            obj.M.N'*dN'                          (1/K)*obj.M.N'*obj.M.N] .*gp.j*gp.w;
         else
         kel = kel + [...
            (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN*obj.M.N
            obj.M.N'*dN'                          (1/K)*obj.M.N'*obj.M.N] .*gp.j*gp.w;
         end
      end
      %% Element Fint
      function Fint = computeFint(obj, gp, el)
         dN = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         K = obj.Bulk;
         if gp.i == 1
         Fint = [(gp.B'*gp.sigma - (K*(dN*dN'))*el.Uvc ); zeros(size(obj.M.N,2),1)] *gp.j *gp.w;
         else
         Fint = el.Fint +...
            [(gp.B'*gp.sigma - (K*(dN*dN'))*el.Uvc ); zeros(size(obj.M.N,2),1)] *gp.j *gp.w;  
         end
      end
   end
end