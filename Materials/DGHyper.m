classdef DGHyper
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      numeq;
      ngp
      numstr
      finiteDisp = 1;
      I
      I4_bulk
      
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
      c1L;
      c1R;
      C1;
      bnAdN1
      bnAdN2
      tvtr
      jumpu
      ep
      
      pencoeff = 4
      eGPL;
      eGPR;
      sGP;
      bGP;
      el;
      
      tauLHist;
      tauRHist;
      
      P1_2 = [...
         1 0 0 0
         0 1 0 0
         0 0 1 0];
      
      P2_2 = [...
         1 0 0    0
         0 0 1/2  1/2
         0 0 1/2 -1/2
         0 1 0    0];
      P3_2 = [...
         1 0 0    0
         0 0 1/2 -1/2
         0 0 1/2  1/2
         0 1 0    0];
      
      P1_3 = [eye(6,6) zeros(6,3)];
      P2_3 = [...
         1 0 0 0 0 0 0 0 0
         0 0 0 1/2 0 0 1/2 0 0
         0 0 0 0 0 1/2 0 0 -1/2
         0 0 0 1/2 0 0 -1/2 0 0
         0 1 0 0 0 0 0 0 0
         0 0 0 0 1/2 0 0 1/2 0
         0 0 0 0 0 1/2 0 0 1/2
         0 0 0 0 1/2 0 0 -1/2 0
         0 0 1 0 0 0 0 0 0];
      P3_3 = [...
         1 0 0 0 0 0 0 0 0
         0 0 0 1/2 0 0 -1/2 0 0
         0 0 0 0 0 1/2 0 0 1/2
         0 0 0 1/2 0 0 1/2 0 0
         0 1 0 0 0 0 0 0 0
         0 0 0 0 1/2 0 0 -1/2 0
         0 0 0 0 0 1/2 0 0 -1/2
         0 0 0 0 1/2 0 0 1/2 0
         0 0 1 0 0 0 0 0 0];
      
      term28L
      term28R
      term30L
      term30R
      
      jumpAddL
      jumpAddR
   end
   %%
   methods
      %% Construct
      function ob = DGHyper(num, props, propsList, identity)
         ob.ndm    = num.ndm;
         ob.ndof   = num.ndof;
         ob.numeq  = num.nen*num.ndm;
         ob.numstr = num.str;
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
               ob.bGP = T3(1, 0, xi, w);
            case 4
               ob.eGPL = Q4(1,0,xiL,w);     ob.eGPR = Q4(1,0,xiR,w);
               ob.sGP  = L3(0);
               
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1 0 +1 0 -1 0
                  -1 -1 +1 +1 -1 0 +1 0 0]';
               w = (1/81).*[25 25 25 25 40 40 40 40 64];
               ob.bGP = Q4(1, 0, xi, w);
            case 8
               ob.eGPL = Q8(1,0,xiL,w);     ob.eGPR = Q8(1,0,xiR,w);
               ob.sGP  = Q4(0);
               
               xi =  1/sqrt(3) .*[...
                  -1  1 -1  1 -1  1 -1  1
                  -1 -1  1  1 -1 -1  1  1
                  -1 -1 -1 -1  1  1  1  1]';
               w = [1 1 1 1 1 1 1 1]';
               ob.bGP = Q8(1, 0, xi, w);
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
         ob.I      = eye(ob.ndm);
         ob.I4_bulk= identity.I4_bulk;
         ob.tauLHist= zeros(ob.ndm, ob.ndm, size(ob.sGP.Nmat,1), num.el);
         ob.tauRHist= zeros(ob.ndm, ob.ndm, size(ob.sGP.Nmat,1), num.el);
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, ~, ~, ~)
         eps = zeros(ob.numstr,1);
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, el, ~)
         ndm  = ob.ndm;
         if ndm == 2
            P1 = ob.P1_2;  P2 = ob.P2_2;  P3 = ob.P3_2;
         elseif ndm == 3
            P1 = ob.P1_3;  P2 = ob.P2_3;  P3 = ob.P3_3;
         end
         nen = size(el.conn(el.i,:),2)/2;
         elL = 1:nen;
         elR = nen+1:2*nen;
         
         lamL = ob.lamdaL;    lamR = ob.lamdaR;
         muL  = ob.muL;        muR = ob.muR;
         
         ob.eGPL.U = el.Umt(elL,:);
         ob.eGPR.U = el.Umt(elR,:);
         ulresL = reshape(el.Umt(elL,:)', numel(el.Umt(elL,:)),1);
         ulresR = reshape(el.Umt(elR,:)', numel(el.Umt(elR,:)),1);
         
         
         iterset = 3;
         if el.iter < iterset
            ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elL));
            tauL = ob.computeTau(ob.bGP, muL, lamL, ob.ndm, class(ob.bGP), ob.eGPL.U);
            
            ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elR));
            tauR = ob.computeTau(ob.bGP, muR, lamR, ob.ndm, class(ob.bGP), ob.eGPR.U);
            
            ob.tauLHist(:, :, gp.i, el.i)= tauL;
            ob.tauRHist(:, :, gp.i, el.i)= tauR;
         else
            tauL = ob.tauLHist(:, :, gp.i, el.i);
            tauR = ob.tauRHist(:, :, gp.i, el.i);
         end
         
         ob.eGPL.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elL)); ob.eGPL.iel = 1;
         ob.eGPR.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elR)); ob.eGPR.iel = 1;
         
         ob.eGPL.i = gp.i;   ob.eGPR.i = gp.i;
         
         TanL = ob.eGPL.dXdxi(:,1:end-1);
         TanR = ob.eGPR.dXdxi(:,1:end-1);
         
         tanL = ob.eGPL.F*TanL;
         tanR = ob.eGPR.F*TanR;
         
         [intedge, ob.C1, ~] = edgeInt(ob.sGP, TanL);
         
         eb = ob.edgeBubbleInt(ob.eGPL.xi, ob.C1, class(ob.eGPL));
         
         [ ~, ob.c1L, nvectL] = edgeInt(ob.sGP, tanL);
         [ ~, ob.c1R, nvectR] = edgeInt(ob.sGP, tanR);
         
         edgeK  = (tauL*eb^2 + tauR*eb^2);
         gamL   = eb^2*(edgeK\tauL);
         gamR   = eb^2*(edgeK\tauR);
         ob.ep = ob.pencoeff*intedge*ob.I/edgeK;
         
         NL = ob.eGPL.N;  NR = ob.eGPR.N;
         pad = zeros(ndm, nen);
         
         ob.NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
         ob.NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
         
         BmatfL = ob.eGPL.Bf;    BmatfR = ob.eGPR.Bf;
         BmatL  = ob.eGPL.B;     BmatR  = ob.eGPR.B;
         
         [sigmaL,  cmatL]  = ob.SigmaCmat2(ob.eGPL.b, ob.eGPL.JxX, muL, lamL, ob.I, ob.I4_bulk);
         [sigmaR,  cmatR]  = ob.SigmaCmat2(ob.eGPR.b, ob.eGPR.JxX, muR, lamR, ob.I, ob.I4_bulk);
         % Make kirchhoff stress into cauchy:
         sigmaL = sigmaL/ob.eGPL.JxX;  cmatL = cmatL/ob.eGPL.JxX;
         sigmaR = sigmaR/ob.eGPR.JxX;  cmatR = cmatR/ob.eGPR.JxX;
         
         nvecL = diag(nvectL(1:ndm,1:ndm));
         nvecR = diag(nvectR(1:ndm,1:ndm));
         tL  = sigmaL*nvecL;   tR  = sigmaR*nvecR;
         DnL = nvectL*cmatL;   DnR = nvectR*cmatR;
         if ndm == 2
            SnL  = [tL zeros(ndm,1); zeros(ndm,1) tL];
            SnR  = [tR zeros(ndm,1); zeros(ndm,1) tR];
            
            cmatnBL=BmatfL'*P2'*[DnL zeros(2,3); zeros(2,3)  DnL];
            cmatnBR=BmatfR'*P2'*[DnR zeros(2,3); zeros(2,3)  DnR];
         elseif ndm == 3
            SnL  = [tL zeros(ndm,2); zeros(ndm,1) tL zeros(ndm,1); zeros(ndm,2) tL];
            SnR  = [tR zeros(ndm,2); zeros(ndm,1) tR zeros(ndm,1); zeros(ndm,2) tR];
            
            cmatnBL=BmatfL'*P2'*[DnL zeros(3,12); zeros(3,6) DnL zeros(3,6); zeros(3,12) DnL];
            cmatnBR=BmatfR'*P2'*[DnR zeros(3,12); zeros(3,6) DnR zeros(3,6); zeros(3,12) DnR];
         end
         
         term17L = P2'*SnL*gamL';
         term17R = P2'*SnR*gamR';
         
         term18L  = P1'*(gamL*nvectL*cmatL)';
         term18R  = P1'*(gamR*nvectR*cmatR)';
         
         ob.term28L = ob.NmatL'*(ob.c1L(gp.i)*gamL*sigmaL*nvecL-ob.c1R(gp.i)*gamR*sigmaR*nvecR);
         ob.term28R = ob.NmatR'*(ob.c1L(gp.i)*gamL*sigmaL*nvecL-ob.c1R(gp.i)*gamR*sigmaR*nvecR);
         
         ob.jumpu = ob.NmatL*ulresL - ob.NmatR*ulresR;
         
         ob.term30L = ob.NmatL'*ob.ep*ob.jumpu;
         ob.term30R = ob.NmatR'*ob.ep*ob.jumpu;
         
         gamajumpuL = gamL'*ob.jumpu;
         gamajumpuR = gamR'*ob.jumpu;
         
         sig8L2 = gamajumpuL'*DnL;
         sig8R2 = gamajumpuR'*DnR;
         if ndm == 2
            term5L=cmatnBL*[eye(3)*gamajumpuL(1) eye(3)*gamajumpuL(2)]'*BmatL;
            term5R=cmatnBR*[eye(3)*gamajumpuR(1) eye(3)*gamajumpuR(2)]'*BmatR;
            
            sig8L3 = [...
               sig8L2(1) 0         sig8L2(3) 0
               sig8L2(3) 0         sig8L2(2) 0
               0         sig8L2(1) 0         sig8L2(3)
               0         sig8L2(3) 0         sig8L2(2)];
            sig8R3 = [...
               sig8R2(1) 0         sig8R2(3) 0
               sig8R2(3) 0         sig8R2(2) 0
               0         sig8R2(1) 0         sig8R2(3)
               0         sig8R2(3) 0         sig8R2(2)];
         elseif ndm == 3
            term5L=cmatnBL*[eye(6)*gamajumpuL(1) eye(6)*gamajumpuL(2) eye(6)*gamajumpuL(3)]'*BmatL;
            term5R=cmatnBR*[eye(6)*gamajumpuR(1) eye(6)*gamajumpuR(2) eye(6)*gamajumpuR(3)]'*BmatR;
            
            sig8L3 = [...
               sig8L2(1) 0         0         sig8L2(4) 0         0         sig8L2(6) 0         0
               sig8L2(4) 0         0         sig8L2(2) 0         0         sig8L2(5) 0         0
               sig8L2(6) 0         0         sig8L2(5) 0         0         sig8L2(3) 0         0
               0         sig8L2(1) 0         0         sig8L2(4) 0         0         sig8L2(6) 0
               0         sig8L2(4) 0         0         sig8L2(2) 0         0         sig8L2(5) 0
               0         sig8L2(6) 0         0         sig8L2(5) 0         0         sig8L2(3) 0
               0         0         sig8L2(1) 0         0         sig8L2(4) 0         0         sig8L2(6)
               0         0         sig8L2(4) 0         0         sig8L2(2) 0         0         sig8L2(5)
               0         0         sig8L2(6) 0         0         sig8L2(5) 0         0         sig8L2(3)];
            sig8R3 = [...
               sig8R2(1) 0         0         sig8R2(4) 0         0         sig8R2(6) 0         0
               sig8R2(4) 0         0         sig8R2(2) 0         0         sig8R2(5) 0         0
               sig8R2(6) 0         0         sig8R2(5) 0         0         sig8R2(3) 0         0
               0         sig8R2(1) 0         0         sig8R2(4) 0         0         sig8R2(6) 0
               0         sig8R2(4) 0         0         sig8R2(2) 0         0         sig8R2(5) 0
               0         sig8R2(6) 0         0         sig8R2(5) 0         0         sig8R2(3) 0
               0         0         sig8R2(1) 0         0         sig8R2(4) 0         0         sig8R2(6)
               0         0         sig8R2(4) 0         0         sig8R2(2) 0         0         sig8R2(5)
               0         0         sig8R2(6) 0         0         sig8R2(5) 0         0         sig8R2(3)];
         end
         term8L = BmatfL'*P2'*sig8L3*(P3*BmatfL);
         term8R = BmatfR'*P2'*sig8R3*(P3*BmatfR);
         
         dmatL1 = ob.dmat2_no_p(ob.eGPL.JxX, muL, lamL, ndm);
         dmatR1 = ob.dmat2_no_p(ob.eGPR.JxX, muR, lamR, ndm);
         
         dmatL2 = reshape(dmatL1/ob.eGPL.JxX*nvectL'*gamL'*ob.jumpu, ob.numstr,ob.numstr);
         dmatR2 = reshape(dmatR1/ob.eGPR.JxX*nvectR'*gamR'*ob.jumpu, ob.numstr,ob.numstr);
         
         term7L = BmatL'*dmatL2*BmatL;
         term7R = BmatR'*dmatR2*BmatR;
         
         ob.bnAdN1 = (term17L+term18L)'*BmatfL;
         ob.bnAdN2 = (term17R+term18R)'*BmatfR;
         
         ob.jumpAddL = term5L+term5L'+term7L+term8L;
         ob.jumpAddR = term5R+term5R'+term7R+term8R;
         
         D       = zeros(6,6);
         sigma_v = zeros(ob.numstr,1);
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;
         bnAdN2 = ob.bnAdN2;
         cL = ob.c1L(i);
         cR = ob.c1R(i);
         C  = ob.C1(i);
         
         mid = size(el.K,1)/2;
         ElemKLL = el.K(    1:mid,     1:mid) + C*(NL'*ob.ep*NL) - cL*NL'*bnAdN1 - cL*bnAdN1'*NL - cL*ob.jumpAddL;
         ElemKLR = el.K(    1:mid, mid+1:end) - C*(NL'*ob.ep*NR) + cR*NL'*bnAdN2 + cL*bnAdN1'*NR;
         ElemKRL = el.K(mid+1:end,     1:mid) - C*(NR'*ob.ep*NL) + cL*NR'*bnAdN1 + cR*bnAdN2'*NL;
         ElemKRR = el.K(mid+1:end, mid+1:end) + C*(NR'*ob.ep*NR) - cR*NR'*bnAdN2 - cR*bnAdN2'*NR + cR*ob.jumpAddR;
         
         Kel = [...
            ElemKLL ElemKLR
            ElemKRL ElemKRR];
         
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         cL = ob.c1L(i);
         cR = ob.c1R(i);
         C = ob.C1(i);
         
         tvtr  = ob.tvtr;
         jumpu = ob.jumpu;
         
         mid = size(el.Fint,1)/2;
         ElemFL = el.Fint(    1:mid) + C*ob.term30L - cL*bnAdN1'*jumpu - ob.term28L;
         ElemFR = el.Fint(mid+1:end) - C*ob.term30R + cR*bnAdN2'*jumpu + ob.term28R;
         
         Fint = [ElemFL; ElemFR];
      end
   end
   methods (Static)
      function [E, nu] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E = props{j,2};
               case 'nu'
                  nu = props{j,2};
            end
         end
      end
      %% Compute tau
      function tau = computeTau(bGP, mu, lam, ndm, elType, U)
         ngp = size(bGP.xi,1);
         bGP.U = U';
         I4_bulk = zeros(6); I4_bulk(1:3,1:3) = ones(3);
         bGP.iel=1; tau = zeros(ndm, ndm); I = eye(ndm);
         for i = 1:ngp
            bGP.i = i;
            dxdxi = bGP.F*bGP.dXdxi;
            B = DGHyper.edgeBubbleB(bGP.xi(i,:), dxdxi, elType);
            
            [S, cmat] = DGHyper.SigmaCmat2(bGP.b, bGP.JxX, mu, lam, I, I4_bulk);
            if ndm == 2
               Dgeo = formGeo(S);
               Dmat = [cmat zeros(3,1); zeros(1,4)];
            elseif ndm == 3
               Dgeo = formGeo(S);
               Dmat = [cmat zeros(6,3); zeros(3,9)];
            end
            
            tau  = tau  + bGP.J *bGP.w* (B'*(Dgeo+Dmat)*B);
         end
         tau = inv(tau);
      end
      %% Compute Edge Bubble shape function' B matrix
      function B = edgeBubbleB(xi, dxdxi, elType)
         ndm = length(xi);
         switch elType
            case 'T3'
               r = xi(1); s = xi(2);
               dbdxi  = 4*[(1-2*r-s), -r];
               dbdx   = dbdxi / dxdxi;
            case 'Q4'
               r = xi(1); s = xi(2);
               dbdxi  = [r*(s-1), 1/2*(r^2-1)];
               dbdx   = dbdxi / dxdxi;
            case 'Q8'
               r = xi(1); s = xi(2); t = xi(3);
               dbdxi  =[-2*r*(1-s^2)*(1-t), -2*s*(1-r^2)*(1-t), -(1-r^2)*(1-s^2)];
               dbdx   = dbdxi / dxdxi;
         end
         if ndm == 2
            B = [...
               dbdx(1) 0       dbdx(2)  dbdx(2)
               0       dbdx(2) dbdx(1) -dbdx(1)]';
         elseif ndm == 3
            B = [...
               dbdx(1) 0       0       dbdx(2) 0       dbdx(3)  dbdx(2)  0       -dbdx(3)
               0       dbdx(2) 0       dbdx(1) dbdx(3) 0       -dbdx(1)  dbdx(3)  0
               0       0       dbdx(3) 0       dbdx(2) dbdx(1)  0       -dbdx(2)  dbdx(1)]';
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
      %% Compute sigma and Cmat
      function [sigma, cmat] = SigmaCmat2(b, JxX, mu, lam, I, I4_bulk)
         ndm  = size(I,1);
         matE = diag([2,2,2,1,1,1]);
         
         sigma = mu*(b-I) + lam*JxX*(JxX-1)*I;
         cmat  = mu*matE  + lam*JxX*( (2*JxX-1)*I4_bulk - (JxX-1)*matE );
         if ndm == 2
            sigma = sigma(1:2,1:2);
            cmat  = cmat([1,2,4],[1,2,4]);
         end
      end
      %% %d_ijklmn term
      function dmat = dmat2_no_p(JxX,mu,lam,ndm)
         if ndm == 2
            dpmat1 = [...
               1 1 0 1 1 0 0 0 0
               1 1 0 1 1 0 0 0 0
               0 0 0 0 0 0 0 0 0]';
            dpmat2 =-[...
               6 2 0 2 2 0 0 0 1
               2 2 0 2 6 0 0 0 1
               0 0 1 0 0 1 1 1 0]';
            dpmat3 = [...
               8 0 0 0 0 0 0 0 2
               0 0 0 0 8 0 0 0 2
               0 0 2 0 0 2 2 2 0]';
         elseif ndm == 3
            I1 = [1; 1; 1; 0; 0; 0];
            cpmat1 = I1*I1';
            dpmat1 = [cpmat1; cpmat1; cpmat1; zeros(18,6)];
            dpmat2 = [...
               -6	-2	-2	 0	 0	 0	-2	-2	 0	 0	 0	 0	-2	 0	-2	 0	 0	 0	 0	 0	 0	-1	0	0	0	0	0	0 -1 0  0  0  0 0	0 -1
               -2	-2	 0	 0	 0	 0	-2	-6	-2	 0	 0	 0	 0	-2	-2	 0	 0	 0	 0	 0	 0	-1	0	0  0	0	0	0 -1 0  0  0  0 0	0 -1
               -2	 0	-2	 0	 0	 0	 0	-2	-2	 0	 0	 0	-2	-2	-6	 0	 0	 0	 0	 0	 0	-1	0	0	0	0	0	0 -1 0  0  0  0 0	0 -1
               0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	-1	-1	-1	 0	0	0	0	0	0	0	0 0  0  0  0 0	0	0
               0	 0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	 0	 0	 0	0	0 -1 -1 -1	0	0 0  0  0  0 0	0	0
               0	 0	 0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	 0	 0	 0	-1	 0	 0	 0	 0	0	0	0  0	0	0	0 0 -1 -1 -1 0	0	0]';
            dpmat3 = [...
               8	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	0	0	2
               0	0	0	0	0	0	0	8	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	2	0	0	0	0	0	0	0
               0	0	0	0	0	0	0	0	0	0	0	0	0	0	8	0	0	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	0	2
               0	0	0	2	0	0	0	0	0	2	0	0	0	0	0	0	0	0	2	2	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0
               0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	0	2	0	0	0	0	0	0	1	0	2	2	0	0	0	0	0	0	1	0	0
               0	0	0	0	0	2	0	0	0	0	0	0	0	0	0	0	0	2	0	0	0	0	1	0	0	0	0	1	0	0	2	0	2	0	0	0]';
         end
         A = 4*JxX^2 - JxX;
         B = 2*JxX^2 - JxX;
         C = 1*JxX^2 - JxX;
         
         dmat = lam*(A*dpmat1 + B*dpmat2 + C*dpmat3) - mu*dpmat3;
      end
   end
end