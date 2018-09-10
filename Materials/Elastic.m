classdef Elastic
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      ngp;
      finiteDisp = 0;
      
      Young;
      Poisson;
      linear = true;
      
      name = 'Elastic';
   end
   %%
   methods
      %% Construct
      function obj = Elastic(num, props)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.ngp  = num.gp;
         
         if strcmp(props{1,1} ,'E') && strcmp(props{2,1} ,'v')
            obj.Young  = props{1,2};
            obj.Poisson= props{2,2};
         elseif strcmp(props{1,1} ,'v') && strcmp(props{2,1} ,'E')
            obj.Poisson= props{1,2};
            obj.Young  = props{2,2};
         else
            error(['You''ve chosen Elastic material but specified ',...
               'incompatible material properties, I''m disapponted']);
         end
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~, ~)
         E = obj.Young;
         v = obj.Poisson;
         
         Eh= E/(1-2*v)/(1+v);
         G = 0.5*E/(1+v);
         D =[...
            Eh*(1-v) Eh*v     Eh*v     0 0 0
            Eh*v     Eh*(1-v) Eh*v     0 0 0
            Eh*v     Eh*v     Eh*(1-v) 0 0 0
            0        0        0        G 0 0
            0        0        0        0 G 0
            0        0        0        0 0 G];
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if obj.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(~, gp, ~, ~)
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
end

