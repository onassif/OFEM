function [globl, Fext, rmG, adG, indc]=ApplyConstraints_and_Loads(mult, globl, Fext, U,...
   inpt, ndof, step, iter, finiteDisp)

BC    = inpt.BC;
FORCE = inpt.FORCE;
indc = true(size(globl.K,1),1);

%%   Loads
for i = 1:size(FORCE,1)
   index=FORCE(i,1);
   direction=FORCE(i,2);
   
   %       Ignore any load that conflicts with BC
   mtch_indx = (BC(:,1) == index);
   if((find(BC(mtch_indx,2) == direction)))
      continue
   end
   
   f_index = ndof*(index-1) + direction; % A smart way to do it without if-statement
   Fext(f_index) = FORCE(i,3)*mult(FORCE(i,4),2);
end
G = Fext - globl.Fint;

%% BC
for i = 1:size(BC,1)
   index = BC(i,1);
   direction = BC(i,2);
   
   %       Locate the index in the global matrix
   k_index = ndof*(index-1) + direction;
   indc(k_index) = false;
   
   %       Change the Fs
   globl.Fint(k_index) = BC(i,3)*mult(BC(i,4),1) - U(k_index);
   G(k_index)    = BC(i,3)*mult(BC(i,4),1) - U(k_index);
end

%% Reduce K and G to only unknown indices
rmG = G(indc);
adG = G(~indc); 
if step == 1 && iter == 0 && finiteDisp
   rmG = rmG - globl.K(indc,~indc)*G(~indc);
end
globl.K = globl.K(indc,indc);
globl.M = globl.M(indc,indc);
end