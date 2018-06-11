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
      
      pencoeff = 4
      sGP;
      bGP;
      el;
      name = 'DG';
   end
   %%
   methods
      %% Construct
      function obj = DG(num, props, propsList, el)
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
         
         xi = [...
            -sqrt(0.6)  0  sqrt(0.6)]';
         w = (1/9).*[5 8 5];
         obj.sGP = L2(0,0,xi,w);
         
         xi = sqrt(0.6)*[...
            -1 +1 +1 -1 0 +1 0 -1 0
            -1 -1 +1 +1 -1 0 +1 0 0]';
         w = (1/81).*[25 25 25 25 40 40 40 40 64];
         obj.bGP = Q4(0, 0, xi, w);
         
         obj.el = el.elements;
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
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, el, ~)
         lamdaR = obj.lamdaR;
         lamdaL = obj.lamdaL;
         muR = obj.muR;
         muL = obj.muL;
         
         DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
         DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
         
         ulresL = [0 0 0 0 0.1 0 0.1 0]';
         ulresR = [0 0 0 0 0.0 0 0.0 0]';
         
         obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,:));
         [tauL, intbL] = obj.computeTau(obj.bGP, DmatL, obj.ndm);
         obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i+1,:));
         [tauR, intbR] = obj.computeTau(obj.bGP, DmatR, obj.ndm);
         
         obj.sGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,1:2));
         for i = 1:size(obj.sGP.xi,1)
            r = obj.sGP.xi(i);
            s = -1;
            bL = 1/2*(1-s)*(1-r^2);
            bR = bL;
         end
            
         
         
         
         obj.bGP;
         %          obj.sGP.dXdxi = [1; -1]'*obj.sGP.dNdxi*;
         %          obj.sGP.det_dXdxi_list(1) = det(obj.sGP.dXdxi);
         
         
         num = struct('el', 2, 'nen', 9, 'gp', 9, 'ndm', 2);
         eGP = Q9(num, 0);
         
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
      %% Compute tau
      function [tau, intb] = computeTau(bGP, D, ndm)
         bGP.iel=1; tau = zeros(ndm, ndm); intb = 0;
         for i = 1:size(bGP.xi,1)
            r = bGP.xi(i,1);
            s = bGP.xi(i,2);
            bGP.i = i;
            dbdxi = [r*(s-1), 1/2*(r^2-1)];
            dbdX  = dbdxi / bGP.dXdxi;
            b  = 1/2*(1-s)*(1-r^2);
            B = [...
               dbdX(1) 0
               0       dbdX(2)
               dbdX(2) dbdX(1)];
            tau = tau   + bGP.J*bGP.w* (B'*D*B);
            intb = intb + bGP.J*bGP.w* b;
         end
         tau = inv(tau);
      end
   end
end