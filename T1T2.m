function [B] = T1T2(A,opt)
if opt==1 %stress-like voigt
   B = [...
      A(1) A(4) A(6)
      A(4) A(2) A(5)
      A(6) A(5) A(3)];
elseif opt==2 %strain-like voigt
   B = [...
      A(1)     1/2*A(4) 1/2*A(6)
      1/2*A(4) A(2)     1/2*A(5)
      1/2*A(6) 1/2*A(5) A(3)];
end
end

