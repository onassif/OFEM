function [B] = T2T1(A,opt)
if length(A)==2 % 2D
   if opt==1 % stress-like voigt
      B = [A(1,1);A(2,2);A(1,2)];
   elseif opt==2 % strain-like voigt
      B = [A(1,1);A(2,2);2*A(1,2)];
   end
elseif length(A)==3 % 3D
   if opt==1 
      B = [A(1,1);A(2,2);A(3,3);A(1,2);A(2,3);A(3,1)];
   elseif opt == 2 
      B = [A(1,1);A(2,2);A(3,3);2*A(1,2);2*A(2,3);2*A(3,1)];
   end
end
end