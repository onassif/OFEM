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
      NmatL;
      NmatR;
      c1;
      bnAdN1
      bnAdN2
      tvtr
      jumpu
      ep
      linear = true;
      
      pencoeff = 4
      sGPL;
      sGPR;
      bGP;
      el;
      name = 'DG';
   end
   %%
   methods
      %% Construct
      function obj = DG(num, props, propsList, el)
         obj.ndm  = num.ndm;
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
            -sqrt(0.6)  0  sqrt(0.6)
            -1 -1 -1]';
         w = (1/9).*[5 8 5];
         obj.sGPL = Q4(0,0,xi,w);
         
         xi = [...
            sqrt(0.6)  0  -sqrt(0.6)
            -1 -1 -1]';
         obj.sGPR = Q4(0,0,xi,w);
         
         xi = sqrt(0.6)*[...
            -1 +1 +1 -1 0 +1 0 -1 0
            -1 -1 +1 +1 -1 0 +1 0 0]';
         w = (1/81).*[25 25 25 25 40 40 40 40 64];
         obj.bGP = Q4(0, 0, xi, w);
         
         obj.el = el.elements;
         
         ngpS = 3;      ndm = 2; nen = 4;
         
         obj.NmatL = zeros(size(el.nodes,2), size(el.nodes,1), ngpS);
         obj.NmatR = zeros(size(el.nodes,2), size(el.nodes,1), ngpS);
         
         obj.c1 = zeros(ngpS, 1);
         
         obj.bnAdN1 = zeros(ndm, ndm*nen, ngpS);
         obj.bnAdN2 = zeros(ndm, ndm*nen, ngpS);
         
         obj.tvtr  = zeros(ndm, ngpS);
         obj.jumpu = zeros(ndm, ngpS);
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, ~, ~, ~)
         if obj.ndm == 2
            eps = zeros(3,1);
         elseif obj.ndm == 3
            eps = zeros(6,1);
         end
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, el, ~)
         ndm  = obj.ndm;
         ngpS = size(obj.sGPL.xi,1);
         
         lamdaL = obj.lamdaL;    lamdaR = obj.lamdaR;
         muL = obj.muL;          muR = obj.muR;
         
         DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
         DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
         
         ulresL = [0 0 0 0 0.1 0 0.1 0]';
         ulresR = [0 0 0 0 0.0 0 0.0 0]';
         
         obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,:));
         [tauL, ~] = obj.computeTau(obj.bGP, DmatL, obj.ndm-1);
         
         obj.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i+1,:));
         [tauR, ~] = obj.computeTau(obj.bGP, DmatR, obj.ndm-1);
         
         obj.sGPL.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,:));
         obj.sGPL.iel = 1;
         [ebL, intedge, obj.c1, nvect] = obj.edgeInt(obj.sGPL);
         
         obj.sGPR.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i+1,:));
         obj.sGPR.iel = 1;
         [ebR, ~, ~, ~] = obj.edgeInt(obj.sGPR);
         
         edgeK  = tauL*ebL^2 + tauR*ebR^2;
         gamL   = ebL^2*(edgeK\tauL);
         gamR   = ebR^2*(edgeK\tauR);
         obj.ep = obj.pencoeff*intedge*inv(edgeK);
         
         nen = size(el.conn,2);
         for i = 1:ngpS
            obj.sGPL.i = i;   obj.sGPR.i = i;
            NL = obj.sGPL.N;  NR = obj.sGPR.N;
            pad = zeros(ndm, nen);
            
            NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
            NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
            
            BmatL = obj.sGPL.B;
            BmatR = obj.sGPR.B;
            
            bnAdN1 = gamL*nvect*DmatL*BmatL;
            bnAdN2 = gamR*nvect*DmatR*BmatR;
            
            obj.tvtr(:,i)  = (bnAdN1*ulresL + bnAdN2*ulresR);
            obj.jumpu(:,i) = NmatR*ulresR - NmatL*ulresL;
            
            obj.NmatL(:,:,i)  = NmatL;    obj.NmatR(:,:,i)  = NmatR;
            obj.bnAdN1(:,:,i) = bnAdN1;   obj.bnAdN2(:,:,i) = bnAdN2;
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
      function Kel = computeK_el(obj, ~, ~, ~)
         ngpS = size(obj.sGPL.xi,1);
         
         ElemKLL = zeros(8,8);
         ElemKLR = zeros(8,8);
         ElemKRL = zeros(8,8);
         ElemKRR = zeros(8,8);
         for i = 1:ngpS
            NL = obj.NmatL(:,:,i);
            NR = obj.NmatR(:,:,i);
            
            c = obj.c1(i);
            
            bnAdN1 = obj.bnAdN1(:,:,i);
            bnAdN2 = obj.bnAdN2(:,:,i);
            
            ElemKLL = ElemKLL + c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*obj.ep*NL));
            ElemKLR = ElemKLR + c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*obj.ep*NR));
            ElemKRL = ElemKRL + c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*obj.ep*NL));
            ElemKRR = ElemKRR + c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*obj.ep*NR));
         end
         Kel = [...
            ElemKLL ElemKLR
            ElemKRL ElemKRR];
      end
      %% Element Fint
      function Fint = computeFint(obj, ~, ~)
         ngpS = size(obj.sGPL.xi,1);
         
         ElemFL = zeros(8,1);
         ElemFR = zeros(8,1);
         
         for i = 1:ngpS
            NL = obj.NmatL(:,:,i);
            NR = obj.NmatR(:,:,i);
            
            c = obj.c1(i);
            
            bnAdN1 = obj.bnAdN1(:,:,i);
            bnAdN2 = obj.bnAdN2(:,:,i);
            
            tvtr  = obj.tvtr(:,i);
            jumpu = obj.jumpu(:,i);
            
            ElemFL = ElemFL - c*( - NL'*(tvtr + obj.ep*jumpu) + bnAdN1'*jumpu);
            ElemFR = ElemFR - c*( + NR'*(tvtr + obj.ep*jumpu) + bnAdN2'*jumpu);
         end
         
         Fint = [ElemFL; ElemFR];
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
         ngp = size(bGP.xi,1);
         bGP.iel=1; tau = zeros(ndm, ndm); intb = 0;
         for i = 1:ngp
            bGP.i = i;
            
            [b,dbdX] = DG.edgeBubble(bGP.xi(i,1), bGP.xi(i,2), bGP.dXdxi);
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
      function [eb, intedge, c1, nvect] = edgeInt(sGP)
         ngp = size(sGP.xi,1);
         c1 = zeros(ngp,1);
         sGP.iel  = 1; eb = 0; intedge = 0;
         for i =1:ngp
            sGP.i = i;
            [ebe, ~] = DG.edgeBubble(sGP.xi(i,1), sGP.xi(i,2), sGP.dXdxi);
            t1 = [sGP.dXdxi(:,1); 0];
            t2 = [0 0 1];
            t3m = norm(cross(t1,t2));
            t3u = cross(t1,t2)/t3m;
            c1(i) = sGP.w*t3m;
            eb      = eb + c1(i)*ebe;
            intedge = intedge + c1(i);
         end
         nvect = [...
            t3u(1)   0        t3u(2)
            0        t3u(2)   t3u(1)];
      end
      %% Compute Edge Bubble shape function
      function [b, dbdX] = edgeBubble(r,s,dXdxi)
         dbdxi = [r*(s-1), 1/2*(r^2-1)];
         
         dbdX  = dbdxi / dXdxi;
         b     = 1/2*(1-s)*(1-r^2);
      end
   end
end