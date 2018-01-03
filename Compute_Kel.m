function Kel = Compute_Kel(Kel, gp, ndof, ngp)
% Definitions
B=gp.B; 
D=gp.D;
sigma=[gp.sigma(1) gp.sigma(3)
       gp.sigma(3) gp.sigma(2)]; 
dNdx=gp.dNdx; 
K_geo = zeros(ndof*ngp,ndof*ngp);
I = eye(ndof);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for a =1:ngp
    for b=1:ngp
        for k=1:ndof
            for i=1:ndof
                for l=1:ndof
                    for j=1:ndof
                        aa=(a-1)*ndof; bb=(b-1)*ndof;
                        
                        K_geo(aa+i, bb+k) = K_geo(aa+i, bb+k) +...
                            I(l,k)*dNdx(a,l) *sigma(i,j)*dNdx(b,j);
                    end
                end
            end
        end
    end
end
Kel = Kel + ( (B'*D*B) + K_geo ) *gp.j*gp.w;
end

