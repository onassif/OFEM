classdef DG
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      numeq;
      ngp
      finiteDisp = 0;
      
      EL;
      ER;
      vL;
      vR;
      lamdaR;
      lamdaL;
      muR;
      muL;
      NmatL;
      NmatR;
      dNdxi
      c1;
      bnAdN1
      bnAdN2
      tvtr
      jumpu
      ep
      toggle = false;
      linear = true;
      
      pencoeff = 4
      eGPL;
      eGPR;
      sGP;
      bGP;
      el;
      name = 'DG';
   end
   %%
   methods
      %% Construct
      function obj = DG(num, props, propsList, el)
         obj.ndm   = num.ndm;
         obj.ndof  = num.ndof;
         obj.numeq = num.nen*num.ndm;
         if num.gp == 1 || num.gp == 4
               obj.ngp = 3;
               xiL = [-sqrt(0.6)  0  sqrt(0.6); -1 -1 -1]';
               xiR = [sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1]';
               w = (1/18).*[5 8 5];
         end
         switch num.gp
            case 1
               xiL = (1+xiL)/2;
               xiR = (1+xiR)/2;
               obj.eGPL = T3(0,0,xiL,w);     obj.eGPR = T3(0,0,xiR,w);
               
               xi = [...
                  1/3 1/3
                  0.05971587179 0.47014206410
                  0.47014206410 0.05971587179
                  0.47014206410 0.47014206410
                  0.79742698540 0.10128650730
                  0.10128650730 0.79742698540
                  0.10128650730 0.10128650730];
               w = [0.1125 0.0662 0.0662 0.0662 0.0630 0.0630 0.0630];
               obj.bGP = T3(0, 0, xi, w);
               obj.dNdxi = [-1;1;0];
            case 4
               obj.eGPL = Q4(0,0,xiL,w);     obj.eGPR = Q4(0,0,xiR,w);
               obj.sGP = L3(0);
               
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1 0 +1 0 -1 0
                  -1 -1 +1 +1 -1 0 +1 0 0]';
               w = (1/81).*[25 25 25 25 40 40 40 40 64];
               obj.bGP = Q4(0, 0, xi, w);
               obj.dNdxi = [-0.394337567;0.394337567;0.105662433;-0.105662433];
            otherwise
               error("unimplemented shape");
         end
         
         for i = 1:2
            switch props{i,1}
               case 'L'
                  [obj.EL, obj.vL] = obj.getProps(propsList{props{i,2}});
               case 'R'
                  [obj.ER, obj.vR] = obj.getProps(propsList{props{i,2}});
            end
         end
         obj.lamdaL = obj.vL*obj.EL/((1+obj.vL)*(1-2*obj.vL));
         obj.lamdaR = obj.vR*obj.ER/((1+obj.vR)*(1-2*obj.vR));
         obj.muL    = obj.EL/(2*(1+obj.vL));
         obj.muR    = obj.ER/(2*(1+obj.vR));
         
         obj.el = el.elements;
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, ~, ~)
         if gp.i == 1
            obj.toggle = ~obj.toggle;
         end
         numstr = (obj.ndm*obj.ndm + obj.ndm) /2;
         
         eps = zeros(numstr,1);
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, el, ~)
         if obj.toggle
            ndm  = obj.ndm;
            
            lamdaL = obj.lamdaL;    lamdaR = obj.lamdaR;
            muL = obj.muL;          muR = obj.muR;
            
            DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
            DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
            
            ulresL = el.ulres(el.i,:)';
            ulresR = el.ulres(el.i+1,:)';
            
            obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,:));
            [tauL, ~] = obj.computeTau(obj.bGP, DmatL, obj.ndm-1, class(obj.bGP));
            
            obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i+1,:));
            [tauR, ~] = obj.computeTau(obj.bGP, DmatR, obj.ndm-1, class(obj.bGP));
            
            obj.eGPL.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,:));   obj.eGPL.iel = 1;
            obj.eGPR.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i+1,:)); obj.eGPR.iel = 1;
            
            surfTan = el.nodes(el.conn(el.i,:),:)'*obj.dNdxi;
            
            [eb, intedge, obj.c1, nvect] = obj.edgeInt(obj.sGP, surfTan);
            
            edgeK  = (tauL + tauR)*eb^2;
            gamL   = eb^2*(edgeK\tauL);
            gamR   = eb^2*(edgeK\tauR);
            obj.ep = obj.pencoeff*intedge*inv(edgeK);
            
            nen = size(el.conn,2);
            obj.eGPL.i = gp.i;   obj.eGPR.i = gp.i;
            NL = obj.eGPL.N;  NR = obj.eGPR.N;
            pad = zeros(ndm, nen);
            
            obj.NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
            obj.NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
            
            BmatL = obj.eGPL.B;
            BmatR = obj.eGPR.B;
            
            obj.bnAdN1 = gamL*nvect*DmatL*BmatL;
            obj.bnAdN2 = gamR*nvect*DmatR*BmatR;
            
            obj.tvtr  = (obj.bnAdN1*ulresL + obj.bnAdN2*ulresR);
            obj.jumpu = obj.NmatR*ulresR - obj.NmatL*ulresL;
         end
         Eh= obj.EL/(1-2*obj.vL)/(1+obj.vL);
         G = 0.5*obj.EL/(1+obj.vL);
         v = obj.vL;
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
      function Kel = computeK_el(obj, Kel, gp, ~)
         i = gp.i;
         
         if obj.toggle
            NL     = obj.NmatL;         NR = obj.NmatR;
            bnAdN1 = obj.bnAdN1;    bnAdN2 = obj.bnAdN2;
            c = obj.c1(i);

            if gp.i == 1
               ElemKLL = c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*obj.ep*NL));
               ElemKLR = c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*obj.ep*NR));
               ElemKRL = c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*obj.ep*NL));
               ElemKRR = c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*obj.ep*NR));
            else
               mid = size(Kel,1)/2;
               ElemKLL = Kel(    1:mid,     1:mid) + c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*obj.ep*NL));
               ElemKLR = Kel(    1:mid, mid+1:end) + c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*obj.ep*NR));
               ElemKRL = Kel(mid+1:end,     1:mid) + c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*obj.ep*NL));
               ElemKRR = Kel(mid+1:end, mid+1:end) + c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*obj.ep*NR));
            end
            Kel = [...
               ElemKLL ElemKLR
               ElemKRL ElemKRR];
         else
            Kel = zeros(obj.numeq, obj.numeq);
         end
         
      end
      %% Element Fint
      function Fint = computeFint(obj, gp, el)
         i = gp.i;
         
         if obj.toggle
            NL     = obj.NmatL;         NR = obj.NmatR;
            bnAdN1 = obj.bnAdN1;    bnAdN2 = obj.bnAdN2;
            c = obj.c1(i);
            
            tvtr  = obj.tvtr;
            jumpu = obj.jumpu;
            if gp.i == 1
               ElemFL = - c*( - NL'*(tvtr + obj.ep*jumpu) + bnAdN1'*jumpu);
               ElemFR = - c*( + NR'*(tvtr + obj.ep*jumpu) + bnAdN2'*jumpu);
            else
               mid = size(el.Fint,1)/2;
               ElemFL = el.Fint(    1:mid) - c*( - NL'*(tvtr + obj.ep*jumpu) + bnAdN1'*jumpu);
               ElemFR = el.Fint(mid+1:end) - c*( + NR'*(tvtr + obj.ep*jumpu) + bnAdN2'*jumpu);
            end
            Fint = [ElemFL; ElemFR];
         else
            Fint = zeros(obj.numeq,1);
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
      %% Compute tau
      function [tau, intb] = computeTau(bGP, D, ndm, elType)
         ngp = size(bGP.xi,1);
         bGP.iel=1; tau = zeros(ndm, ndm); intb = 0;
         for i = 1:ngp
            bGP.i = i;
            
            [b,dbdX] = DG.edgeBubble(bGP.xi(i,1), bGP.xi(i,2), bGP.dXdxi, elType);
            B = [...
               dbdX(1) 0
               0       dbdX(2)
               dbdX(2) dbdX(1)];
            tau = tau   + bGP.J*bGP.w* (B'*D*B);
            intb = intb + bGP.J*bGP.w* b;
         end
         tau = inv(tau);
      end
      %% Integrate normal vectors
      function [eb, intedge, c1, nvect] = edgeInt(oGP, surfTan)
         oGP.iel = 1;
         J = norm(surfTan);
         n = cross([surfTan;0],[0;0;1] )/J;
         nvect = [...
            n(1) 0    n(2)
            0    n(2) n(1)];
        c1 = J * oGP.weights;
        eb = sum(J * oGP.weights .* oGP.Nmat(:,2));
        intedge = sum(J * oGP.weights);         
      end
      %% Compute Edge Bubble shape function
      function [b, dbdX] = edgeBubble(r,s,dXdxi,elType)
         switch elType
            case 'T3'
               dbdX = 4*[(1-2*r-s), -r];
               b = 4*(1-r-s)*r;
            case 'Q4'
               dbdxi = [r*(s-1), 1/2*(r^2-1)];
               dbdX  = dbdxi / dXdxi;
               b     = 1/2*(1-s)*(1-r^2);
         end
      end
   end
end