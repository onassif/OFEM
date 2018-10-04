function [Dgeo] = formGeo(S)
%formGeo Forms the geometric part of D for finite displacements

%% convert to voigt 
if size(S,2)==2
   S = [S(1,1) S(2,2) S(1,2)]';
elseif size(S,2)==3
   S = [S(1,1) S(2,2) S(3,3) S(1,2) S(2,3) S(3,1)]';
end
%% Form the matrix
if (size(S,2)==1 && size(S,1)==3)
   Dgeo = [...
        S(1)       0        S(3)/2        S(3)/2
           0    S(2)        S(3)/2       -S(3)/2
      S(3)/2  S(3)/2 (S(2)+S(1))/4 (S(2)-S(1))/4
      S(3)/2 -S(3)/2 (S(2)-S(1))/4 (S(2)+S(1))/4];
elseif (size(S,2)==1 && size(S,1)==6)
   Dgeo = [...
      S(1)       0       0          S(4)/2               0          S(6)/2          S(4)/2               0         -S(6)/2
         0    S(2)       0          S(4)/2          S(5)/2               0         -S(4)/2          S(5)/2               0
         0       0    S(3)               0          S(5)/2          S(6)/2               0         -S(5)/2          S(6)/2
    S(4)/2  S(4)/2       0 S(1)/4 + S(2)/4          S(6)/4          S(5)/4 S(2)/4 - S(1)/4          S(6)/4         -S(5)/4
         0  S(5)/2  S(5)/2          S(6)/4 S(2)/4 + S(3)/4          S(4)/4         -S(6)/4 S(3)/4 - S(2)/4          S(4)/4
    S(6)/2       0  S(6)/2          S(5)/4          S(4)/4 S(1)/4 + S(3)/4          S(5)/4         -S(4)/4 S(1)/4 - S(3)/4
    S(4)/2 -S(4)/2       0 S(2)/4 - S(1)/4         -S(6)/4          S(5)/4 S(1)/4 + S(2)/4         -S(6)/4         -S(5)/4
         0  S(5)/2 -S(5)/2          S(6)/4 S(3)/4 - S(2)/4         -S(4)/4         -S(6)/4 S(2)/4 + S(3)/4         -S(4)/4
   -S(6)/2       0  S(6)/2         -S(5)/4          S(4)/4 S(1)/4 - S(3)/4         -S(5)/4         -S(4)/4 S(1)/4 + S(3)/4];
else
   error('unknown stress matrix dimensions');
end
end

