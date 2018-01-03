function [K, Fint] = Assemble(K, Fint, el)
ind = el.indices;

Fint(ind) = Fint(ind) + el.Fint;
K(ind, ind) = K(ind, ind) + el.K;

end