function [K, Fint] = Assemble(K, Fint, el, iter)
ind = el.indices;
if iter ~=0
    Fint(ind) = Fint(ind) + el.Fint;
end
K(ind, ind) = K(ind, ind) + el.K;

end