classdef DGCP
   properties (SetAccess = private)
      ndm;
      ndof;
      numeq;
      ngp
      numstr
      finiteDisp = 1;
      I
      I4_bulk
      listL = struct('D',[], 'Rp',[], 'R', [], 'S',[], 'tauT', []);
      listR = struct('D',[], 'Rp',[], 'R', [], 'S',[], 'tauT', []);
      
      matL
      matR
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
      function ob = DGCP(num, props, propsList, identity)
         ob.ndm    = num.ndm;
         ob.ndof   = num.ndof;
         ob.numeq  = num.nen*num.ndm;
         ob.numstr = num.str;
         [ob.eGPL,ob.eGPR,ob.bGP,ob.sGP,ob.ngp] = DGxi(num.nen,num.ndm,1);
         ob.listL.D = zeros(   6, 6, num.gp, num.el2, num.steps+1);
         ob.listL.S = zeros(      6, num.gp, num.el2, num.steps+1);
         ob.listL.tauT = zeros(      num.gp, num.el2, num.steps+1);
         ob.listL.Rp   = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.listL.R    = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.listL.Rp(1,1,:,:,1) = 1; ob.listL.Rp(2,2,:,:,1) = 1; ob.listL.Rp(3,3,:,:,1) = 1;
         ob.listL.R( 1,1,:,:,1) = 1; ob.listL.R( 2,2,:,:,1) = 1; ob.listL.R( 3,3,:,:,1) = 1;
         
         ob.listR.D = zeros(   6, 6, num.gp, num.el2, num.steps+1);
         ob.listR.S = zeros(      6, num.gp, num.el2, num.steps+1);
         ob.listR.tauT = zeros(      num.gp, num.el2, num.steps+1);
         ob.listR.Rp   = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.listR.R    = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.listR.Rp(1,1,:,:,1) = 1; ob.listR.Rp(2,2,:,:,1) = 1; ob.listR.Rp(3,3,:,:,1) = 1;
         ob.listR.R( 1,1,:,:,1) = 1; ob.listR.R( 2,2,:,:,1) = 1; ob.listR.R( 3,3,:,:,1) = 1;
         
         for i = 1:2
            switch props{i,1}
               case 'L'
                  ob.matL = props{i,2};
               case 'R'
                  ob.matR = props{i,2};
            end
         end
         
         ob.I       = eye(ob.ndm);
         ob.I4_bulk = identity.I4_bulk;
         ob.tauLHist= zeros(ob.ndm, ob.ndm, size(ob.sGP.Nmat,1), num.el);
         ob.tauRHist= zeros(ob.ndm, ob.ndm, size(ob.sGP.Nmat,1), num.el);
      end
      %% Epsilon
      function [eps, ob] = Strain(ob, ~, ~, ~)
         eps = zeros(ob.numstr,1);
      end
      %% Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, el, step)
         ndm  = ob.ndm;
         sigma_v = zeros(6,1);
         if ndm == 2
            P1 = ob.P1_2;  P2 = ob.P2_2;  P3 = ob.P3_2;
         elseif ndm == 3
            P1 = ob.P1_3;  P2 = ob.P2_3;  P3 = ob.P3_3;
         end
         nen = size(el.conn(el.i,:),2)/2;
         elL = 1:nen;
         elR = nen+1:2*nen;
         
         ob.eGPL.U = el.Umt(:,elL);
         ob.eGPR.U = el.Umt(:,elR);
         ulresL = reshape(el.Ures(:,elL), numel(el.Ures(:,elL)),1);
         ulresR = reshape(el.Ures(:,elR), numel(el.Ures(:,elR)),1);
         
         coorL = el.nodes(el.conn(el.i, elL),:)';
         coorR = el.nodes(el.conn(el.i, elR),:)';
         [xlintL,xlintR, drdrL,~, ob.eGPL,ob.eGPR] = intBounds2(coorL,coorR, ob.eGPL,ob.eGPR);
         
         iterset = 3;
         if el.iter < iterset
         [ob.bGP.det_dXdxi_list, ob.bGP.dNdX_list, ob.bGP.dXdxi_list] = shapeRef(...
            xlintL', 1:nen, ob.bGP.dNdxi_list);  
            ob.bGP.U=ob.eGPL.U;	ob.bGP.iel=1;
            ob.bGP.dU  = el.Ures(:,elL);
            tauL = ob.computeTau(ob.bGP, el.mat{ob.matL}, ob.ndm);
            
         [ob.bGP.det_dXdxi_list, ob.bGP.dNdX_list, ob.bGP.dXdxi_list] = shapeRef(...
            xlintR', 1:nen, ob.bGP.dNdxi_list);  
            ob.bGP.U=ob.eGPR.U;
            ob.bGP.dU  = el.Ures(:,elR);
            tauR = ob.computeTau(ob.bGP, el.mat{ob.matR}, ob.ndm);
            
            ob.tauLHist(:, :, gp.i, el.i)= tauL;
            ob.tauRHist(:, :, gp.i, el.i)= tauR;
         else
            tauL = ob.tauLHist(:, :, gp.i, el.i);
            tauR = ob.tauRHist(:, :, gp.i, el.i);
         end
         
         ob.eGPL.iel = 1;  ob.eGPL.i = gp.i;
         ob.eGPR.iel = 1;  ob.eGPR.i = gp.i;
         [ob.eGPL.det_dXdxi_list, ob.eGPL.dNdX_list, ob.eGPL.dXdxi_list] = shapeRef(...
            el.nodes, el.conn(el.i, elL), ob.eGPL.dNdxi_list);
         [ob.eGPR.det_dXdxi_list, ob.eGPR.dNdX_list, ob.eGPR.dXdxi_list] = shapeRef(...
            el.nodes, el.conn(el.i, elR), ob.eGPR.dNdxi_list);
         
         TanL = ob.eGPL.dXdxi(:,1:end-1);
         TanR = ob.eGPR.dXdxi(:,1:end-1);
         
         tanL = ob.eGPL.F*TanL;
         tanR = ob.eGPR.F*TanR;
         
         [intedge, ob.C1, ~] = edgeInt(ob.sGP, TanL, drdrL);
         
         eb = ob.eGPL.bubb*ob.C1;
         
         [ ~, ob.c1L, nvectL] = edgeInt(ob.sGP, tanL, drdrL);
         [ ~, ob.c1R, nvectR] = edgeInt(ob.sGP, tanR, drdrR);
         
         edgeK = (tauL*eb^2 + tauR*eb^2);
         gamL  = eb^2*(edgeK\tauL);
         gamR  = eb^2*(edgeK\tauR);
         ob.ep = ob.pencoeff*intedge*ob.I/edgeK;
         
         NL = ob.eGPL.N';  NR = ob.eGPR.N';
         pad = zeros(ndm, nen);
         
         ob.NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
         ob.NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
         
         BmatfL = ob.eGPL.Bf;    BmatfR = ob.eGPR.Bf;
         RmatL  = ob.eGPL.R;     RmatR  = ob.eGPR.R;
         
         ob.eGPL.U   = el.Umt_n(:,elL) + 1/2*el.w(:,elL);
         ob.eGPR.U   = el.Umt_n(:,elR) + 1/2*el.w(:,elR);
         
         Bmat2L  = ob.eGPL.B;    Bmat2R  = ob.eGPR.B;
         R2L     = ob.eGPL.R;    R2R     = ob.eGPR.R;
         
         cgn1L = ob.listL.S(:,gp.i,gp.iel,step);
         cgn1R = ob.listR.S(:,gp.i,gp.iel,step);
         if step == 1 && el.iter == 0
            QL = Qmat(el.mat{ob.matL}.gRot');
            QR = Qmat(el.mat{ob.matR}.gRot');
            qn1L   = Qmat(check2D(RmatL));	qn1R   = Qmat(check2D(RmatR));
            sigmaL = qn1L*cgn1L;             sigmaR = qn1R*cgn1R;
            
            cmatL = QL*el.mat{ob.matL}.C0*QL';
            cmatR = QR*el.mat{ob.matR}.C0*QR';
         else
            elm = struct('i',el.i,'iter',el.iter);
            if ndm == 2
               Q = Qmat([R2L zeros(2,1);zeros(1,2) 1]);
               Q = Q([1,2,4],[1,2,4]);
               de = Q'*Bmat2L*ulresL;
               el.mat{ob.matL}.de = [de(1) de(2) 0 de(3) 0 0]';
            elseif ndm == 3
               el.mat{ob.matL}.de = Qmat(R2L)'*Bmat2L*ulresL;
            end
            el.mat{ob.matL}.list = ob.listL;
            [sigmaL, cmatL, el.mat{ob.matL}] = el.mat{ob.matL}.SigmaCmat(ob.eGPL, elm, step);
            ob.listL = el.mat{ob.matL}.list;
            
            if ndm == 2
               Q = Qmat([R2R zeros(2,1);zeros(1,2) 1]);
               Q = Q([1,2,4],[1,2,4]);
               de = Q'*Bmat2R*ulresR;
               el.mat{ob.matR}.de = [de(1) de(2) 0 de(3) 0 0]';
            elseif ndm == 3
               el.mat{ob.matR}.de = Qmat(R2R)'*Bmat2R*ulresR;
            end
            el.mat{ob.matR}.list = ob.listR;
            [sigmaR, cmatR, el.mat{ob.matR}] = el.mat{ob.matR}.SigmaCmat(ob.eGPR, elm, step);
            ob.listR = el.mat{ob.matR}.list;
         end
         sigmaL = T1T2(sigmaL,1);
         sigmaR = T1T2(sigmaR,1);
         if ndm == 2
            sigmaL = sigmaL(1:2,1:2);        sigmaR = sigmaR(1:2,1:2);
            cmatL  = cmatL([1,2,4],[1,2,4]); cmatR  = cmatR([1,2,4],[1,2,4]);
         end
         
         nvecL = diag(nvectL(1:ndm,1:ndm));
         nvecR = diag(nvectR(1:ndm,1:ndm));
         tL  = sigmaL*nvecL;   tR  = sigmaR*nvecR;
         if ndm == 2
            SnL  = [tL zeros(ndm,1); zeros(ndm,1) tL];
            SnR  = [tR zeros(ndm,1); zeros(ndm,1) tR];
         elseif ndm == 3
            SnL  = [tL zeros(ndm,2); zeros(ndm,1) tL zeros(ndm,1); zeros(ndm,2) tL];
            SnR  = [tR zeros(ndm,2); zeros(ndm,1) tR zeros(ndm,1); zeros(ndm,2) tR];
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
         
         ob.bnAdN1 = (term17L+term18L)'*BmatfL;
         ob.bnAdN2 = (term17R+term18R)'*BmatfR;
         
         D = zeros(6);
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
         ElemKLL = el.K(    1:mid,     1:mid) + C*(NL'*ob.ep*NL) - cL*NL'*bnAdN1 - cL*bnAdN1'*NL;
         ElemKLR = el.K(    1:mid, mid+1:end) - C*(NL'*ob.ep*NR) + cR*NL'*bnAdN2 + cL*bnAdN1'*NR;
         ElemKRL = el.K(mid+1:end,     1:mid) - C*(NR'*ob.ep*NL) + cL*NR'*bnAdN1 + cR*bnAdN2'*NL;
         ElemKRR = el.K(mid+1:end, mid+1:end) + C*(NR'*ob.ep*NR) - cR*NR'*bnAdN2 - cR*bnAdN2'*NR;
         
         Kel = [...
            ElemKLL ElemKLR
            ElemKRL ElemKRR];
         
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, ~)
         i = gp.i;
         
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         cL = ob.c1L(i);
         cR = ob.c1R(i);
         C  = ob.C1(i);
         
         mid = size(el.Fint,1)/2;
         ElemFL = el.Fint(    1:mid) + C*ob.term30L - cL*bnAdN1'*ob.jumpu - ob.term28L;
         ElemFR = el.Fint(mid+1:end) - C*ob.term30R + cR*bnAdN2'*ob.jumpu + ob.term28R;
         
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
      function tau = computeTau(bGP, CP, ndm)
         Q    = Qmat(CP.gRot');
         cmat = Q*CP.C0*Q';
         
         ngp = size(bGP.xi,2);
         tau = zeros(ndm, ndm);
         for i = 1:ngp
            bGP.i = i;
            B = bGP.bubbB;
            D = formCombD(zeros(6,1), cmat, bGP.finiteDisp);
            if ndm == 2
               D = D([1,2,4,7],[1,2,4,7]);
            end
            tau  = tau  + bGP.J*bGP.w* (B'*D*B);
         end
         tau = inv(tau);
      end
   end
end