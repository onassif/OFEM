function Dcomb = formCombD(S,D,finiteDisp)
%formGeo Forms the geometric part of D for finite displacements
if ~finiteDisp % Small Strains
   Dcomb = D;
elseif finiteDisp
   %% convert to voigt first
   if     size(S,2)==2 && min(size(S))~= 1
      S1=S(1,1); S2=S(2,2); S3=S(1,2);
   elseif length(S)==3 && min(size(S))== 1
      S1=S(1); S2=S(2); S3=S(3);
   elseif size(S,2)==3 && min(size(S))~= 1
      S1=S(1,1); S2=S(2,2); S3=S(3,3); S4=S(1,2); S5=S(2,3); S6=S(1,3);
   elseif length(S)==6 && min(size(S))== 1
      S1=S(1); S2=S(2); S3=S(3); S4=S(4); S5=S(5); S6=S(6);
   end
   %% Form the matrix
   if ~exist('S6','var')
      Dgeo = [...
         S1    0     S3/2      S3/2
         0     S2    S3/2     -S3/2
         S3/2  S3/2 (S2+S1)/4 (S2-S1)/4
         S3/2 -S3/2 (S2-S1)/4 (S2+S1)/4];
      if size(D,2) == 6
         D = D([1,2,4],[1,2,4]);
      end
      Dmat = [D zeros(3,1); zeros(1,4)];
   elseif exist('S6','var')
      Dgeo = [...
         +S1     0      0      S4/2       0          S6/2       S4/2       0         -S6/2
         +0      S2     0      S4/2       S5/2       0         -S4/2       S5/2       0
         +0      0      S3     0          S5/2       S6/2       0         -S5/2       S6/2
         +S4/2   S4/2   0     (S1+S2)/4   S6/4       S5/4      (S2-S1)/4   S6/4      -S5/4
         +0      S5/2   S5/2   S6/4      (S2+S3)/4   S4/4      -S6/4      (S3-S2)/4   S4/4
         +S6/2   0      S6/2   S5/4       S4/4      (S1+S3)/4   S5/4      -S4/4      (S1-S3)/4
         +S4/2  -S4/2   0     (S2-S1)/4  -S6/4       S5/4      (S1+S2)/4  -S6/4      -S5/4
         +0      S5/2  -S5/2   S6/4      (S3-S2)/4  -S4/4      -S6/4      (S2+S3)/4  -S4/4
         -S6/2   0      S6/2  -S5/4       S4/4      (S1-S3)/4  -S5/4      -S4/4      (S1+S3)/4];
      
      Dmat = [D zeros(6,3); zeros(3,9)];
   end
   Dcomb = Dmat + Dgeo;
end
end