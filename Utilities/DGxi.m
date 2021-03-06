function [sGPL,sGPR,eGPL,eGPR,bGP,ngp] = DGxi(nen, ndm, finiteDisp)
%DGxi This computes the isoparametric coordinates related to DG

switch nen
   case 3
      ngp = 3;
      e_xiL = [-sqrt(0.6)  0   sqrt(0.6); -1 -1 -1];
      e_xiR = [ sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1];
      e_w = (1/9).*[5 8 5]';
      
      b_xi = [...
         1/3 1/3
         0.05971587179 0.47014206410
         0.47014206410 0.05971587179
         0.47014206410 0.47014206410
         0.79742698540 0.10128650730
         0.10128650730 0.79742698540
         0.10128650730 0.10128650730]';
      b_w = [0.1125 0.0662 0.0662 0.0662 0.0630 0.0630 0.0630]';
 
      sGPL = L2(0);
      sGPR = L2(0);
      eGPL = T3(finiteDisp);
      eGPR = T3(finiteDisp);
      bGP  = T3(finiteDisp);
   case 4
      if ndm == 2
         ngp = 3;
         e_xiL = [-sqrt(0.6)  0   sqrt(0.6); -1 -1 -1];
         e_xiR = [ sqrt(0.6)  0  -sqrt(0.6); -1 -1 -1];
         e_w = (1/9).*[5 8 5]';
      
         b_xi = sqrt(0.6)*[...
            -1 +1 +1 -1 0 +1 0 -1 0
            -1 -1 +1 +1 -1 0 +1 0 0];
         b_w = (1/81).*[25 25 25 25 40 40 40 40 64];
      
         sGPL = L3(0);
         sGPR = L3(0);
         eGPL = Q4(finiteDisp);
         eGPR = Q4(finiteDisp);
         bGP  = Q4(finiteDisp);
      elseif ndm == 3
         ngp = 3;
         e_xiL = 1/6 .* [...
            1	4	1
            1	1	4
            0	0	0];
         e_xiR = 1/6 .* [...
            1	1	4
            1	4	1
            0	0	0];
         e_w = 1/6 .* [1 1 1]';
         
         b_xi = [...
            0.1381966011 0.5854101966 0.1381966011 0.1381966011 
            0.1381966011 0.1381966011 0.5854101966 0.1381966011   
            0.1381966011 0.1381966011 0.1381966011 0.5854101966];
         b_w = 1/24 .*[1 1 1 1]; 
         
         sGPL = T3(0);
         sGPR = T3(0);
         eGPL = T4(finiteDisp);
         eGPR = T4(finiteDisp);
         bGP  = T4(finiteDisp);
      end
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
      
      sGPL = Q4(0);
      sGPR = Q4(0);
      eGPL = Q8(finiteDisp);
      eGPR = Q8(finiteDisp);
      bGP  = Q8(finiteDisp);
   case 10
      ngp = 13;
      e_xiL = [...
         0.06513010290 0.86973979420 0.06513010290 0.3128654960 0.6384441886 0.0486903154 0.6384441886 0.3128654960 0.0486903154 0.2603459661 0.4793080678 0.2603459661 1/3 
         0.06513010290 0.06513010290 0.86973979420 0.0486903154 0.3128654960 0.6384441886 0.0486903154 0.6384441886 0.3128654960 0.2603459661 0.2603459661 0.4793080678 1/3
         0             0             0             0            0            0            0            0            0            0            0            0            0];
      e_xiR = [...
         .06513010290 .86973979420 .06513010290 .3128654960 .6384441886 .0486903154 .6384441886 .3128654960 .0486903154 .2603459661 .4793080678 .2603459661 1/3 
         .06513010290 .06513010290 .86973979420 .0486903154 .3128654960 .6384441886 .0486903154 .6384441886 .3128654960 .2603459661 .2603459661 .4793080678 1/3
         0            0            0            0           0           0           0           0           0           0           0           0           0];
      e_xiR = e_xiR(:,[1,3,2,9,8,7,6,5,4,10,12,11,13]);
     e_w = [...
         0.02667361780 0.02667361780 0.02667361780 .03855688045 .03855688045 .03855688045 .03855688045 .03855688045 .03855688045 .08780762872 .08780762872 .08780762872 -0.07478502223]'; 
      
     b_xi = [...
.2146028713	.3561913862	.2146028713	.2146028713	.04067395853 .87797812440 .04067395853	.04067395853 .3223378901 .03298632957 .32233789010 .32233789010 .06366100188 .60300566480 .06366100188 .06366100188 .06366100188 .26967233150 .06366100188 .26967233150 .60300566480 .06366100188 .60300566480 .26967233150
.2146028713	.2146028713	.3561913862	.2146028713	.04067395853 .04067395853 .87797812440	.04067395853 .3223378901 .32233789010 .03298632957 .32233789010 .06366100188 .06366100188 .60300566480 .06366100188 .26967233150 .06366100188 .60300566480 .06366100188 .26967233150 .26967233150 .06366100188 .60300566480
.2146028713	.2146028713	.2146028713	.3561913862	.04067395853 .04067395853 .04067395853	.87797812440 .3223378901 .32233789010 .32233789010 .03298632957 .60300566480 .06366100188 .06366100188 .26967233150 .06366100188 .06366100188 .26967233150 .60300566480 .06366100188 .60300566480 .26967233150 .06366100188];
%      b_xi = b_xi(:,[3,2,1,4,6,5,7,10,9,8]);
     b_w = [...
        .00665379171 .00665379171 .00665379171 .00665379171 .001679535176 .001679535176 .001679535176 .001679535176 .009226196924 .009226196924 .009226196924 .009226196924 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286 .008035714286]';
              
     sGPL = T6(0);
     sGPR = T6(0); 
     eGPL = T10(finiteDisp);
     eGPR = T10(finiteDisp);
     bGP  = T10(finiteDisp);
end
eGPL.xi = e_xiL; eGPL.weights = e_w;
eGPR.xi = e_xiR; eGPR.weights = e_w;
bGP.xi  = b_xi ; bGP.weights  = b_w;

sGPL = sGPL.shapeIso();
sGPR = sGPR.shapeIso();
eGPL = eGPL.shapeIso();
eGPR = eGPR.shapeIso();
bGP  = bGP.shapeIso( );
if nen == 10
   for i = 1:24
      bGP.dNdxi_list(:,:,i) = bGP.dNdxi_list([3,2,1,4,6,5,7,10,9,8],:,i);
   end
end
end