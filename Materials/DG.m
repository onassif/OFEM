classdef DG
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      numeq;
      ngp
      numstr
      finiteDisp = 0;
      C0L;
      C0R;
      
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
      C1;
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
         ob.numstr= num.str;
         [ob.eGPL,ob.eGPR,ob.bGP,ob.ngp] = DGxi(num.nen,num.ndm,0);
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
         
         GL =   0.5*ob.EL/(1+  ob.vL);
         GR =   0.5*ob.ER/(1+  ob.vR);
         KL = (1/3)*ob.EL/(1-2*ob.vL);
         KR = (1/3)*ob.ER/(1-2*ob.vR);
         ob.C0L = KL*identity.I4_bulk + 2*GL*identity.I4_dev;
         ob.C0R = KR*identity.I4_bulk + 2*GR*identity.I4_dev;
         if ob.ndm == 2
            ob.C0L = ob.C0L([1,2,4],[1,2,4]);
            ob.C0R = ob.C0R([1,2,4],[1,2,4]);
         end
      end
      %% Epsilon
      function [eps, ob] = Strain(ob, ~, ~, ~)
         eps = zeros(ob.numstr,1);
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, el, ~)
         ndm  = ob.ndm;
         nen = size(el.conn(el.i,:),2)/2;
         elL = 1:nen;
         elR = nen+1:2*nen;
         I = eye(ndm);

         DmatL = ob.C0L;
         DmatR = ob.C0R;
         
         ulresL = reshape(el.w(:,elL), numel(el.w(:,elL)),1);
         ulresR = reshape(el.w(:,elR), numel(el.w(:,elL)),1);
         
         coorL = el.nodes(el.conn(el.i, elL),:)';
         coorR = el.nodes(el.conn(el.i, elR),:)';
         [xlintL,xlintR, drdrL,~, ob.eGPL,ob.eGPR] = intBounds2(coorL,coorR, ob.eGPL,ob.eGPR);

         [ob.bGP.det_dXdxi_list, ob.bGP.dNdX_list, ob.bGP.dXdxi_list] = shapeRef(...
            xlintL', 1:nen, ob.bGP.dNdxi_list);
         tauL = ob.computeTau(ob.bGP, DmatL, ob.ndm);
         

         [ob.bGP.det_dXdxi_list, ob.bGP.dNdX_list, ob.bGP.dXdxi_list] = shapeRef(...
            xlintR', 1:nen, ob.bGP.dNdxi_list);
         tauR = ob.computeTau(ob.bGP, DmatR, ob.ndm);
         
         ob.eGPL.iel = 1;  ob.eGPL.i = gp.i;
         ob.eGPR.iel = 1;  ob.eGPR.i = gp.i;
         [ob.eGPL.det_dXdxi_list, ob.eGPL.dNdX_list, ob.eGPL.dXdxi_list] = shapeRef(...
            el.nodes, el.conn(el.i, elL), ob.eGPL.dNdxi_list);
         [ob.eGPR.det_dXdxi_list, ob.eGPR.dNdX_list, ob.eGPR.dXdxi_list] = shapeRef(...
            el.nodes, el.conn(el.i, elR), ob.eGPR.dNdxi_list);
         
         TanL = ob.eGPL.dXdxi(:,1:end-1);
         
         [intedge, ob.C1, nvect] = edgeInt(ob.eGPL, TanL, drdrL);
         
         eb = ob.eGPL.bubb*ob.C1;

         edgeK = (tauL*eb^2 + tauR*eb^2);
         gamL  = eb^2*(edgeK\tauL);
         gamR  = eb^2*(edgeK\tauR);
         ob.ep = ob.pencoeff*intedge*I/edgeK;
         
         ob.eGPL.i = gp.i;   ob.eGPR.i = gp.i;
         NL = ob.eGPL.N';  NR = ob.eGPR.N';
         pad = zeros(ndm, nen);
         
         ob.NmatL = reshape([ NL; repmat([pad; NL],ndm-1,1) ], ndm, ndm*nen);
         ob.NmatR = reshape([ NR; repmat([pad; NR],ndm-1,1) ], ndm, ndm*nen);
         
         BmatL = ob.eGPL.B;
         BmatR = ob.eGPR.B;
         
         ob.bnAdN1 = gamL*nvect*DmatL*BmatL;
         ob.bnAdN2 = gamR*nvect*DmatR*BmatR;
         
         ob.tvtr  = (ob.bnAdN1*ulresL + ob.bnAdN2*ulresR);
         ob.jumpu = ob.NmatR*ulresR - ob.NmatL*ulresL;
         
         D = zeros(6,6);
         sigma_v = zeros(ob.numstr,1);
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         c = ob.C1(i);
         
         mid = size(el.K,1)/2;
         ElemKLL = el.K(    1:mid,     1:mid) + c*( - NL'*bnAdN1 - bnAdN1'*NL + (NL'*ob.ep*NL));
         ElemKLR = el.K(    1:mid, mid+1:end) + c*( - NL'*bnAdN2 + bnAdN1'*NR - (NL'*ob.ep*NR));
         ElemKRL = el.K(mid+1:end,     1:mid) + c*( + NR'*bnAdN1 - bnAdN2'*NL - (NR'*ob.ep*NL));
         ElemKRR = el.K(mid+1:end, mid+1:end) + c*( + NR'*bnAdN2 + bnAdN2'*NR + (NR'*ob.ep*NR));
         
         Kel = [...
            ElemKLL ElemKLR
            ElemKRL ElemKRR];
         
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, ~)
         i = gp.i;
         
         NL     = ob.NmatL;         NR = ob.NmatR;
         bnAdN1 = ob.bnAdN1;    bnAdN2 = ob.bnAdN2;
         c = ob.C1(i);
         
         tvtr  = ob.tvtr;
         jumpu = ob.jumpu;
         
         mid = size(el.Fint,1)/2;
         ElemFL = el.Fint(    1:mid) + c*( - NL'*(tvtr + ob.ep*jumpu) + bnAdN1'*jumpu);
         ElemFR = el.Fint(mid+1:end) + c*( + NR'*(tvtr + ob.ep*jumpu) + bnAdN2'*jumpu);
         
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
      function tau = computeTau(bGP, D, ndm)
         ngp = size(bGP.xi,2);
         bGP.iel=1; tau = zeros(ndm, ndm);
         for i = 1:ngp
            bGP.i = i;
            B = bGP.bubbB;
            
            tau  = tau + bGP.J*bGP.w* (B'*D*B);
         end
         tau = inv(tau);
      end
   end
end