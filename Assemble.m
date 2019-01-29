function [K, M, Fint] = Assemble(K, M, Fint, el)
ind = el.indices;
Fint(ind) = Fint(ind) + el.Fint;

K(ind, ind) = K(ind, ind) + el.K;

M(ind, ind) = M(ind, ind) + el.M;
end