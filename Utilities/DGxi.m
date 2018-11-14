function [eGPL,eGPR,bGP,sGP,ngp] = DGxi(nen, finiteDisp)
%DGxi This computes the isoparametric coordinates related to DG

switch nen
   case 3
      ngp = 3;
      e_xiL = [-sqrt(0.6)  0   sqrt(0.6); -1 -1 -1];
      e_xiR = [ sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1];
      e_w = (1/18).*[5 8 5];
      
      b_xi = [...
         1/3 1/3
         0.05971587179 0.47014206410
         0.47014206410 0.05971587179
         0.47014206410 0.47014206410
         0.79742698540 0.10128650730
         0.10128650730 0.79742698540
         0.10128650730 0.10128650730]';
      b_w = [0.1125 0.0662 0.0662 0.0662 0.0630 0.0630 0.0630]';
      
      eGPL = T3(finiteDisp);
      eGPR = T3(finiteDisp);
      bGP  = T3(finiteDisp);
      sGP  = L3(0);
   case 4
      ngp = 3;
      e_xiL = [-sqrt(0.6)  0   sqrt(0.6); -1 -1 -1];
      e_xiR = [ sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1];
      e_w = (1/18).*[5 8 5];
      
      b_xi = sqrt(0.6)*[...
         -1 +1 +1 -1 0 +1 0 -1 0
         -1 -1 +1 +1 -1 0 +1 0 0];
      b_w = (1/81).*[25 25 25 25 40 40 40 40 64];
      
      eGPL = Q4(finiteDisp);
      eGPR = Q4(finiteDisp);
      bGP  = Q4(finiteDisp);
      sGP  = L3(0);
   case 8
      ngp = 4;
      e_xiL = [...
         -sqrt(1/3)  sqrt(1/3) -sqrt(1/3)  sqrt(1/3)
         -sqrt(1/3) -sqrt(1/3)  sqrt(1/3)  sqrt(1/3)
         -1         -1         -1         -1        ];
      e_xiR = [...
         -sqrt(1/3) -sqrt(1/3)  sqrt(1/3)  sqrt(1/3)
         -sqrt(1/3)  sqrt(1/3) -sqrt(1/3)  sqrt(1/3)
         -1         -1         -1         -1        ];
      e_w = [1 1 1 1]';
      
      b_xi =  1/sqrt(3) .*[...
         -1  1 -1  1 -1  1 -1  1
         -1 -1  1  1 -1 -1  1  1
         -1 -1 -1 -1  1  1  1  1];
      b_w = [1 1 1 1 1 1 1 1]';
      
      eGPL = Q8(finiteDisp);
      eGPR = Q8(finiteDisp);
      bGP  = Q8(finiteDisp);
      sGP  = Q4(0);
end
eGPL.xi = e_xiL; eGPL.weights = e_w;
eGPR.xi = e_xiR; eGPR.weights = e_w;
bGP.xi  = b_xi ; bGP.weights  = b_w;

eGPL = eGPL.shapeIso();
eGPR = eGPR.shapeIso();
bGP  = bGP.shapeIso( );
end