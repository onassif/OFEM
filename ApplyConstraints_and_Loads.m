function [K, Fext, Fint]=ApplyConstraints_and_Loads(mult, K, Fext, Fint, inpt, ndof, iter)
BC = inpt.BC;
FORCE = inpt.FORCE;
%   BC
for i = 1:size(BC,1)
    index = BC(i,1);
    direction = BC(i,2);
    
    %       Locate the index in the global matrix
    k_index = ndof*(index-1) + direction;
    
    %       Change the global stiffness matrix
    K(k_index,:) = 0.0;
    K(k_index,k_index) = 1.0;
    
    %       Change the Fs
    Fext(k_index) = BC(i,3)*mult(BC(i,4),1);
    if (iter ~= 0 || BC(i,3)==0)
        Fint(k_index) =BC(i,3)*mult(BC(i,4),1); % Apply BC only if not the zeroth iteration
    end
end

%   Loads
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
end