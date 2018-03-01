function [K, Fext, Fint, G]=ApplyConstraints_and_Loads(mult, K, Fext, Fint, U, inpt, ndof)
BC = inpt.BC;
FORCE = inpt.FORCE;

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
G = Fext - Fint;

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
%     Fext(k_index) = BC(i,3)*mult(BC(i,4),1);
Fint(k_index) = BC(i,3)*mult(BC(i,4),1) - U(k_index);
G(k_index)    = BC(i,3)*mult(BC(i,4),1) - U(k_index);
% G(k_index) = BC(i,3)*mult(BC(i,4),1) - U(k_index) - w(k_index);
% G(k_index) = BC(i,3)*mult(BC(i,4),1) - w(k_index);
% G(k_index) = BC(i,3)*mult(BC(i,4),1);
end
end