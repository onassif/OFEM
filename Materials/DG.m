classdef DG
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      numeq;
      ngp
      finiteDisp = 0;
      C0;
      
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
      function ob = DG(num, props, propsList, identity)
         ob.ndm   = num.ndm;
         ob.ndof  = num.ndof;
         ob.numeq = num.nen*num.ndm;
         if num.nen == 3 || num.nen == 4
            ob.ngp = 3;
            xiL = [-sqrt(0.6)  0   sqrt(0.6); -1 -1 -1]';
            xiR = [ sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1]';
            w = (1/18).*[5 8 5];
         elseif num.nen == 8
            ob.ngp = 4;
            xiL = [...
               -sqrt(1/3)  sqrt(1/3) -sqrt(1/3)  sqrt(1/3)
               -sqrt(1/3) -sqrt(1/3)  sqrt(1/3)  sqrt(1/3)
               -1         -1         -1         -1        ]';
            xiR = [...
               -sqrt(1/3) -sqrt(1/3)  sqrt(1/3)  sqrt(1/3)
               -sqrt(1/3)  sqrt(1/3) -sqrt(1/3)  sqrt(1/3)
               -1         -1         -1         -1        ]';
            w = [1 1 1 1]';
         end
         switch num.nen
            case 3
               xiL = (1+xiL)/2;
               xiR = (1+xiR)/2;
               ob.eGPL = T3(0,0,xiL,w);     ob.eGPR = T3(0,0,xiR,w);
               
               xi = [...
                  1/3 1/3
                  0.05971587179 0.47014206410
                  0.47014206410 0.05971587179
                  0.47014206410 0.47014206410
                  0.79742698540 0.10128650730
                  0.10128650730 0.79742698540
                  0.10128650730 0.10128650730];
               w = [0.1125 0.0662 0.0662 0.0662 0.0630 0.0630 0.0630];
               ob.bGP = T3(0, 0, xi, w);
            case 4
               ob.eGPL = Q4(0,0,xiL,w);     ob.eGPR = Q4(0,0,xiR,w);
               ob.sGP  = L3(0);
               
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1 0 +1 0 -1 0
                  -1 -1 +1 +1 -1 0 +1 0 0]';
               w = (1/81).*[25 25 25 25 40 40 40 40 64];
               ob.bGP = Q4(0, 0, xi, w);
            case 8
               ob.eGPL = Q8(0,0,xiL,w);     ob.eGPR = Q8(0,0,xiR,w);
               ob.sGP  = Q4(0);
               
               xi =  1/sqrt(3) .*[...
                  -1  1 -1  1 -1  1 -1  1
                  -1 -1  1  1 -1 -1  1  1
                  -1 -1 -1 -1  1  1  1  1]';
               w = [1 1 1 1 1 1 1 1]';
               ob.bGP = Q8(0, 0, xi, w);
            otherwise
               error("unimplemented shape");
         end
         
         for i = 1:2
            switch props{i,1}
               case 'L'
                  [ob.EL, ob.vL] = ob.getProps(propsList{props{i,2}});
               case 'R'
                  [ob.ER, ob.vR] = ob.getProps(propsList{props{i,2}});
            end
         end
         ob.lamdaL = ob.vL*ob.EL/((1+ob.vL)*(1-2*ob.vL));
         ob.lamdaR = ob.vR*ob.ER/((1+ob.vR)*(1-2*ob.vR));
         ob.muL    = ob.EL/(2*(1+ob.vL));
         ob.muR    = ob.ER/(2*(1+ob.vR));
         
         G =   0.5*ob.EL/(1+  ob.vL);
         K = (1/3)*ob.EL/(1-2*ob.vL);
         ob.C0 = K*identity.I4_bulk + 2*G*identity.I4_dev;
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, gp, ~, ~)
         numstr = (ob.ndm*ob.ndm + ob.ndm) /2;
         
         eps = zeros(numstr,1);
      end
      %% Sigma
      function [sigma_voigt, ob] = computeCauchy(ob, gp, ~, ~)
         sigma_voigt = gp.D*gp.eps;
      end
      %% Tangential stiffness
      function [D, ctan, ob] = computeTangentStiffness(ob, gp, el, ~)
         ndm  = ob.ndm;
         nen = size(el.conn(el.i,:),2)/2;
         elL = 1:nen;
         elR = nen+1:2*nen;
         I = eye(ndm);
         
         lamdaL = ob.lamdaL;    lamdaR = ob.lamdaR;
         muL = ob.muL;          muR = ob.muR;
         if ndm == 2
            DmatL = muL*diag([2 2 1]) + lamdaL*[1; 1; 0]*[1 1 0];
            DmatR = muR*diag([2 2 1]) + lamdaR*[1; 1; 0]*[1 1 0];
         elseif ndm == 3
            DmatL = muL*diag([2 2 2 1 1 1]) + lamdaL*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
            DmatR = muR*diag([2 2 2 1 1 1]) + lamdaR*[1; 1; 1; 0; 0; 0]*[1 1 1 0 0 0];
         end
         
         ulresL = reshape(el.w(elL,:)', numel(el.w(elL,:)),1);
         ulresR = reshape(el.w(elR,:)', numel(el.w(elR,:)),1);
         
         ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elL));
         tauL = ob.computeTau(ob.bGP, DmatL, ob.ndm, class(ob.bGP));
         
         ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elR));
         tauR = ob.computeTau(ob.bGP, DmatR, ob.ndm, class(ob.bGP));
         
         ob.eGPL.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elL)); ob.eGPL.iel = 1;
         ob.eGPR.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elR)); ob.eGPR.iel = 1;
         
         ob.eGPL.i = gp.i;   ob.eGPR.i = gp.i;
         
         TanL = ob.eGPL.dXdxi(:,1:end-1);
         TanR = ob.eGPR.dXdxi(:,1:end-1);
         
         [intedge, ob.c1, nvect] = ob.edgeInt(ob.sGP, TanL);
         
         eb = ob.edgeBubbleInt(ob.eGPL.xi, ob.c1, class(ob.eGPL));
         
         edgeK = (tauL*eb^2 + tauR*eb^2);
         gamL  = eb^2*(edgeK\tauL);
         gamR  = eb^2*(edgeK\tauR);
         ob.ep = ob.pencoeff*intedge*I/edgeK;
         
         ob.eGPL.i = gp.i;   ob.eGPR.i = gp.i;
         NL = ob.eGPL.N;  NR = ob.eGPR.N;
         pad = zeros(ndm, nen);
         
         ob.NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
         ob.NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
         
         BmatL = ob.eGPL.B;
         BmatR = ob.eGPR.B;
         
         ob.bnAdN1 = gamL*nvect*DmatL*BmatL;
         ob.bnAdN2 = gamR*nvect*DmatR*BmatR;
         
         ob.tvtr  = (ob.bnAdN1*ulresL + ob.bnAdN2*ulresR);
         ob.jumpu = ob.NmatR*ulresR - ob.NmatL*ulresL;
         
         D = ob.C0;
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if ob.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         c = ob.c1(i);
         
         if gp.i == 1
            ElemKLL = c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*ob.ep*NL));
            ElemKLR = c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*ob.ep*NR));
            ElemKRL = c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*ob.ep*NL));
            ElemKRR = c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*ob.ep*NR));
         else
            mid = size(el.K,1)/2;
            ElemKLL = el.K(    1:mid,     1:mid) + c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*ob.ep*NL));
            ElemKLR = el.K(    1:mid, mid+1:end) + c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*ob.ep*NR));
            ElemKRL = el.K(mid+1:end,     1:mid) + c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*ob.ep*NL));
            ElemKRR = el.K(mid+1:end, mid+1:end) + c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*ob.ep*NR));
         end
         Kel = [...
            ElemKLL ElemKLR
            ElemKRL ElemKRR];
         
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         c = ob.c1(i);
         
         tvtr  = ob.tvtr;
         jumpu = ob.jumpu;
         if gp.i == 1
            ElemFL = + c*( - NL'*(tvtr + ob.ep*jumpu) + bnAdN1'*jumpu);
            ElemFR = + c*( + NR'*(tvtr + ob.ep*jumpu) + bnAdN2'*jumpu);
         else
            mid = size(el.Fint,1)/2;
            ElemFL = el.Fint(    1:mid) + c*( - NL'*(tvtr + ob.ep*jumpu) + bnAdN1'*jumpu);
            ElemFR = el.Fint(mid+1:end) + c*( + NR'*(tvtr + ob.ep*jumpu) + bnAdN2'*jumpu);
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
               case 'nu'
                  v = props{j,2};
            end
         end
      end
      %% Compute tau
      function tau = computeTau(bGP, D, ndm, elType)
         ngp = size(bGP.xi,1);
         bGP.iel=1; tau = zeros(ndm, ndm);
         for i = 1:ngp
            bGP.i = i;
            B = DG.edgeBubbleB(bGP.xi(i,:), bGP.dXdxi, elType);

            tau  = tau  + bGP.J*bGP.w* (B'*D*B);
         end
         tau = inv(tau);
      end
      %% Integrate normal vectors
      function [intedge, c1, nvect] = edgeInt(sGP, surfTan)
         sGP.iel = 1;
         ndm = size(surfTan,2)+1;
         if ndm == 2
            J = norm(surfTan);
            n = cross([surfTan;0],[0;0;1] )/J;
            nvect = [...
               n(1) 0    n(2)
               0    n(2) n(1)];
         elseif ndm == 3
            J = abs(det(surfTan(2:3,:)));
            n = cross(surfTan(:,2),surfTan(:,1))/J;
            nvect = [...
               n(1) 0    0    n(2) 0    n(3)
               0    n(2) 0    n(1) n(3) 0
               0    0    n(3) 0    n(2) n(1)];
         end
         
         c1 = J * sGP.weights;
         intedge = sum(J * sGP.weights); % length of the edge
      end
      %% Compute Edge Bubble shape function' B matrix
      function B = edgeBubbleB(xi, dXdxi, elType)
         ndm = length(xi);
         switch elType
            case 'T3'
               r = xi(1); s = xi(2);
               dbdxi  = 4*[(1-2*r-s), -r];
               dbdX   = dbdxi / dXdxi;
            case 'Q4'
               r = xi(1); s = xi(2);
               dbdxi  = [r*(s-1), 1/2*(r^2-1)];
               dbdX   = dbdxi / dXdxi;
            case 'Q8'
               r = xi(1); s = xi(2); t = xi(3);
               dbdxi  =[-2*r*(1-s^2)*(1-t), -2*s*(1-r^2)*(1-t), -(1-r^2)*(1-s^2)];
               dbdX   = dbdxi / dXdxi;
         end
         if ndm == 2
            B = [...
               dbdX(1) 0       dbdX(2)
               0       dbdX(2) dbdX(1)]';
         elseif ndm == 3
            B = [...
               dbdX(1) 0       0       dbdX(2) 0       dbdX(3)
               0       dbdX(2) 0       dbdX(1) dbdX(3) 0      
               0       0       dbdX(3) 0       dbdX(2) dbdX(1)]';
         end
      end
      %% Compute integral of bubble at the edge
      function intb = edgeBubbleInt(xi, C1, elType)
         switch elType
            case 'T3'
               r = xi(:,1); s = xi(:,2);
               bubble = 4*(1-r-s).*r;
            case 'Q4'
               r = xi(:,1); s = xi(:,2);
               bubble = 1/2*(1-s).*(1-r.^2);
            case 'Q8'
               r = xi(:,1); s = xi(:,2); t = xi(:,3);
               bubble = (1-r.^2).*(1-s.^2).*(1-t);
         end
         intb = sum(C1.*bubble);
      end
   end
end