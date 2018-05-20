classdef MixedElasticPlaneStrain
   %PlaneStrain Computes Cauchy stress based on 2D plane strain
   %   Currently only if ndm=2, ndof=2
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      nen;
      dNdX;
      diff;
      finiteDisp = 0;
      
      Young;
      Poisson;
      Bulk;
      
      npress = 1;
      M = 1;
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
         obj.diff = obj.ndof - obj.ndm;
         
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
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~)
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
         dN = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         K = obj.Bulk;
         kel = kel + [...
            (gp.B'*gp.D*gp.B) - (K*(dN*dN'))    dN
            dN'                          (1/K)*obj.M*obj.M] .*gp.j*gp.w;
      end
      %% Element Fint
      function Fint = computeFint(obj, gp, el)
         dN = reshape(gp.dNdx',obj.nen*obj.ndm,1);
         K = obj.Bulk;
         
         Fint = el.Fint +...
            [(gp.B'*gp.sigma - (K*(dN*dN'))*el.Uvc ); 0] *gp.j *gp.w;
      end
   end
end