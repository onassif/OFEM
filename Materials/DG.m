classdef DG
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      finiteDisp = 0;
      
      EL;
      ER;
      vL;
      vR;
      lamdaR;
      lamdaL;
      muR;
      muL;
      linear = true;
      
      sGP;
      name = 'DG';
   end
   %%
   methods
      %% Construct
      function obj = DG(num, props, propsList)
         obj.ndm  = num.ndm - 1;
         obj.ndof = num.ndof;
         for i = 1:2
            switch props{i,1}
               case 'L'
                  [obj.EL, obj.vL] = obj.getProps(propsList{props{i,2}});
               case 'R'
                  [obj.ER, obj.vR] = obj.getProps(propsList{props{i,2}});
            end
         end
         
      obj.lamdaR = obj.vR*obj.ER/((1+obj.vR)*(1-2*obj.vR));
      obj.lamdaL = obj.vL*obj.EL/((1+obj.vL)*(1-2*obj.vL));
      obj.muR    = obj.ER/(2*(1+obj.vR));
      obj.muL    = obj.EL/(2*(1+obj.vL));
      
      obj.sGP = L2(num,0
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, ~, ~, ~)
         if obj.ndm == 1
            eps = zeros(3,1);
         elseif obj.ndm == 2
            eps = zeros(6,1);
         end
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~)
         lamdaR = obj.lamdaR;
         lamdaL = obj.lamdaL;
         muR = obj.muR;
         muL = obj.muL;

         DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
         DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
         
         
         
      end
      %% Element K
      function Kel = computeK_el(~, Kel, gp, ~)
         Kel = Kel + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el)
         Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
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
   
