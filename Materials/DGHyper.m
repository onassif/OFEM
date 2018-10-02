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
      dNdxi
      c1;
      C1;
      bnAdN1
      bnAdN2
      tvtr
      jumpu
      ep
      linear = true;
      
      pencoeff = 4
      eGPL;
      eGPR;
      sGP;
      bGP;
      el;
      
      P1 = [...
         1 0 0 0
         0 1 0 0
         0 0 1 0];
      P2 = [...
         1 0 0    0
         0 0 1/2  1/2
         0 0 1/2 -1/2
         0 1 0    0];
      P3 = [...
         1 0 0    0
         0 0 1/2 -1/2
         0 0 1/2  1/2
         0 1 0    0];
      
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
         if num.gp == 1 || num.gp == 4
            ob.ngp = 3;
            xiL = [-sqrt(0.6)  0  sqrt(0.6); -1 -1 -1]';
            xiR = [sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1]';
            w = (1/18).*[5 8 5];
         end
         switch num.gp
            case 1
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
               ob.dNdxi = [-1;1;0];
            case 4
               ob.eGPL = Q4(1,0,xiL,w);     ob.eGPR = Q4(1,0,xiR,w);
               ob.sGP = L3(0);
               
               xi = sqrt(0.6)*[...
                  -1 +1 +1 -1 0 +1 0 -1 0
                  -1 -1 +1 +1 -1 0 +1 0 0]';
               w = (1/81).*[25 25 25 25 40 40 40 40 64];
               ob.bGP = Q4(1, 0, xi, w);
               ob.dNdxi = [-0.394337567;0.394337567;0.105662433;-0.105662433];
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
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, ~, ~, ~)
         eps = zeros(ob.numstr,1);
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
         
         lamL = ob.lamdaL;    lamR = ob.lamdaR;
         muL  = ob.muL;        muR = ob.muR;
                  
         ob.eGPL.U = el.Umt(elL,:);
         ob.eGPR.U = el.Umt(elR,:);
         ulresL = reshape(el.Umt(elL,:)', numel(el.Umt(elL,:)),1);
         ulresR = reshape(el.Umt(elR,:)', numel(el.Umt(elR,:)),1);
         
         ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elL));
         [tauL, ~] = ob.computeTau(ob.bGP, muL, lamL, ob.ndm-1, class(ob.bGP), ob.eGPL.U);
         
         ob.bGP.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i, elR));
         [tauR, ~] = ob.computeTau(ob.bGP, muR, lamR, ob.ndm-1, class(ob.bGP), ob.eGPR.U);
         
         ob.eGPL.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elL)); ob.eGPL.iel = 1;
         ob.eGPR.mesh = struct('nodes', el.nodes, 'conn', el.conn(el.i,elR)); ob.eGPR.iel = 1;
         
         ob.eGPL.i = gp.i;   ob.eGPR.i = gp.i;
         
         TanL = el.nodes(el.conn(el.i, elL),:)'*ob.dNdxi;
         TanR = el.nodes(el.conn(el.i, elR),:)'*ob.dNdxi;
         
         tanL = ob.eGPL.F*TanL;
         tanR = ob.eGPR.F*TanR;
                 
         [eb, intedge, ob.C1, ~] = ob.edgeInt(ob.sGP, TanL);
         [ ~,       ~,     ~, ~] = ob.edgeInt(ob.sGP, TanR);
         
         [~, ~, ob.c1, nvectL] = ob.edgeInt(ob.sGP, tanL);
         [~, ~,     ~, nvectR] = ob.edgeInt(ob.sGP, tanR);
         
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
                  
         [sigmaL,  cmatL]  = ob.SigmaCmat(ob.eGPL.b, ob.eGPL.JxX, muL, lamL, ob.I, ob.I4_bulk);
         [sigmaR,  cmatR]  = ob.SigmaCmat(ob.eGPR.b, ob.eGPR.JxX, muR, lamR, ob.I, ob.I4_bulk);
         
         nvecL = nvectL(1,[1,3])';
         nvecR = nvectR(1,[1,3])';
         SnL  = [sigmaL*nvecL zeros(ndm,1); zeros(ndm,1) sigmaL*nvecL];
         SnR  = [sigmaR*nvecR zeros(ndm,1); zeros(ndm,1) sigmaR*nvecR];
         
         DnL = nvectL*cmatL;
         DnR = nvectR*cmatR;
         
         cmatnBL=BmatfL'*ob.P2'*[DnL zeros(2,3); zeros(2,3)  DnL];
         cmatnBR=BmatfR'*ob.P2'*[DnR zeros(2,3); zeros(2,3)  DnR];
         
         term17L = ob.P2'*SnL*gamL';
         term17R = ob.P2'*SnR*gamR';
         
         term18L  = ob.P1'*(gamL*nvectL*cmatL)';
         term18R  = ob.P1'*(gamR*nvectR*cmatR)';
         
         ob.term28L = ob.NmatL'*(gamL*sigmaL*nvecL-gamR*sigmaR*nvecR);
         ob.term28R = ob.NmatR'*(gamL*sigmaL*nvecL-gamR*sigmaR*nvecR);
         
         ob.jumpu = ob.NmatR*ulresR - ob.NmatL*ulresL;
         
         ob.term30L = ob.NmatL'*ob.ep*ob.jumpu;
         ob.term30R = ob.NmatR'*ob.ep*ob.jumpu;
         
         gamajumpuL = gamL'*ob.jumpu;
         gamajumpuR = gamR'*ob.jumpu;
         
         term5L=cmatnBL*[eye(3)*gamajumpuL(1) eye(3)*gamajumpuL(2)]'*BmatL;
         term5R=cmatnBR*[eye(3)*gamajumpuR(1) eye(3)*gamajumpuR(2)]'*BmatR;
         
         sig8L2 = [gamajumpuL(1) gamajumpuL(2)]*DnL;
         sig8R2 = [gamajumpuR(1) gamajumpuR(2)]*DnR;
         
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
         
         term8L = BmatfL'*ob.P2'*sig8L3*(ob.P3*BmatfL);
         term8R = BmatfR'*ob.P2'*sig8R3*(ob.P3*BmatfR);
         
         dmatL1 = ob.dmat2_no_p(ob.eGPL.JxX, muL, lamL);
         dmatR1 = ob.dmat2_no_p(ob.eGPR.JxX, muR, lamR);
         
         dmatL2 = reshape(dmatL1/ob.eGPL.JxX*nvectL'*gamL'*ob.jumpu, ob.numstr,ob.numstr);
         dmatR2 = reshape(dmatR1/ob.eGPR.JxX*nvectR'*gamR'*ob.jumpu, ob.numstr,ob.numstr);
         
         term7L = BmatL'*dmatL2*BmatL;
         term7R = BmatR'*dmatR2*BmatR;
         
         ob.bnAdN1 = (term17L+term18L)'*BmatfL;
         ob.bnAdN2 = (term17R+term18R)'*BmatfR;
         
         ob.jumpAddL = term5L+term5L'+term7L+term8L;
         ob.jumpAddR = term5R+term5R'+term7R+term8R;
         
         D = cmatL;
         D =[...
            D(1,1) D(1,2) D(1,2) 0      0      0
            D(2,1) D(2,2) D(1,2) 0      0      0
            D(2,1) D(2,1) D(2,2) 0      0      0
            0      0      0      D(3,3) 0      0
            0      0      0      0      D(3,3) 0
            0      0      0      0      0      D(3,3)];
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if ob.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    
         bnAdN2 = ob.bnAdN2;
         c = ob.c1(i);
         C = ob.C1(i);
         if gp.i == 1
            ElemKLL = + C*(NL'*ob.ep*NL) + c*( - NL'*bnAdN1 - bnAdN1'*NL - ob.jumpAddL);
            ElemKLR = - C*(NL'*ob.ep*NR) + c*( + NL'*bnAdN2 + bnAdN1'*NR);
            ElemKRL = - C*(NR'*ob.ep*NL) + c*( + NR'*bnAdN1 + bnAdN2'*NL);
            ElemKRR = + C*(NR'*ob.ep*NR) + c*( - NR'*bnAdN2 - bnAdN2'*NR + ob.jumpAddR);
         else
            mid = size(el.K,1)/2;
            ElemKLL = el.K(    1:mid,     1:mid) + C*(NL'*ob.ep*NL) + c*( - NL'*bnAdN1 - bnAdN1'*NL - ob.jumpAddL);
            ElemKLR = el.K(    1:mid, mid+1:end) - C*(NL'*ob.ep*NR) + c*( + NL'*bnAdN2 + bnAdN1'*NR);
            ElemKRL = el.K(mid+1:end,     1:mid) - C*(NR'*ob.ep*NL) + c*( + NR'*bnAdN1 + bnAdN2'*NL);
            ElemKRR = el.K(mid+1:end, mid+1:end) + C*(NR'*ob.ep*NR) + c*( - NR'*bnAdN2 - bnAdN2'*NR + ob.jumpAddR);
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
         C = ob.C1(i);
         
         tvtr  = ob.tvtr;
         jumpu = ob.jumpu;
         if gp.i == 1
            ElemFL = + C*ob.term30L + c*( - ob.term28L - bnAdN1'*jumpu);
            ElemFR = - C*ob.term30R + c*( + ob.term28R + bnAdN2'*jumpu);
         else
            mid = size(el.Fint,1)/2;
            ElemFL = el.Fint(    1:mid) + C*ob.term30L + c*( - ob.term28L - bnAdN1'*jumpu);
            ElemFR = el.Fint(mid+1:end) - C*ob.term30R + c*( + ob.term28R + bnAdN2'*jumpu);
         end
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
      function [tau, intb] = computeTau(bGP, mu, lam, ndm, elType, U)
         ngp = size(bGP.xi,1);
         bGP.U = U';
         I4_bulk = zeros(6); I4_bulk(1:3,1:3) = ones(3);
         bGP.iel=1; tau = zeros(ndm, ndm); intb = 0; I = eye(ndm+1);
         for i = 1:ngp
            bGP.i = i;
            [b,dbdX] = DGHyper.edgeBubble(bGP.xi(i,1), bGP.xi(i,2), bGP.dXdxi, elType);
            dbdx = dbdX/bGP.F;
            B = [...
               dbdx(1)  0
               0        dbdx(2)
               dbdx(2)  dbdx(1)
               dbdx(2) -dbdx(1)];
            [S, cmat] = DGHyper.SigmaCmat(bGP.b, bGP.JxX, mu, lam, I, I4_bulk);
            Dgeo = bGP.JxX.*[...
               S(1,1)    0         S(1,2)/2          S(1,2)/2
               0         S(2,2)    S(1,2)/2         -S(1,2)/2
               S(1,2)/2  S(1,2)/2 (S(2,2)+S(1,1))/4 (S(2,2)-S(1,1))/4
               S(1,2)/2 -S(1,2)/2 (S(2,2)-S(1,1))/4 (S(2,2)+S(1,1))/4];
            Dmat = bGP.JxX.*[cmat zeros(3,1); zeros(1,4)];
            
            tau  = tau  + bGP.J/bGP.JxX *bGP.w* (B'*(Dgeo+Dmat)*B);
            intb = intb + bGP.J/bGP.JxX *bGP.w* b;
         end
         tau = inv(tau);
      end
      %% Integrate normal vectors
      function [eb, intedge, c1, nvect] = edgeInt(sGP, surfTan)
         sGP.iel = 1;
         J = norm(surfTan);
         n = cross([surfTan;0],[0;0;1] )/J;
         nvect = [...
            n(1) 0    n(2)
            0    n(2) n(1)];
         c1 = J * sGP.weights;
         eb = sum(J * sGP.weights .* sGP.Nmat(:,2));
         intedge = sum(J * sGP.weights); % length of the edge
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
      %% Compute sigma and Cmat
      function [sigma, cmat] = SigmaCmat(b, JxX, mu, lam, I, I4_bulk)
         ndm  = size(I,1);
         matE = diag([2,2,2,1,1,1]);
         
         sigma = 1/JxX*mu*(b-I) + lam*(JxX-1)*I;
         cmat  = 1/JxX*mu*matE  + lam*( (2*JxX-1)*I4_bulk - (JxX-1)*matE );
         if ndm == 2
            sigma = sigma(1:2,1:2);
            cmat  = cmat([1,2,4],[1,2,4]);
         end
      end
      %% %d_ijklmn term
      function dmat = dmat2_no_p(JxX,mu,lam)
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
         
         A = 4*JxX^2 - JxX;
         B = 2*JxX^2 - JxX;
         C = 1*JxX^2 - JxX;
         
         dmat = lam*(A*dpmat1 + B*dpmat2 + C*dpmat3) - mu*dpmat3;
      end
   end
end