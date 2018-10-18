function [B] = T2T4(A)
if length(A) == 3 % 2D
   B = reshape(A([1,3,3,2],[1,3,3,2]),2,2,2,2);
elseif length(A) == 6 % 3D
   B = reshape(A([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
end
end

