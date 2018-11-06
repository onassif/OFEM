function B = check2D(A)
if size(A,2) == 2
   B = [A, zeros(2,1);zeros(1,2) 1];
elseif size(A,2) == 3
   B = A;
end
end

