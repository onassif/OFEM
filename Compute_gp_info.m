function gp = Compute_gp_info(gp, coor, U, iel)

%   Compute Shape function derivatives in reference coordinates
gp = ShpFun(gp, coor, iel);

%   Compute deformation gradient F = dx/dX 
gp.F = U'*gp.dNdX + eye(size(U,2));

%   Compute Shape function derivatives in spacial coordinates
gp.dNdx = gp.dNdX / gp.F ;

%   Compute left Cauchy-Green tensor
gp.b = gp.F*gp.F';

%   Compute J in both reference and spatial configurations
gp.det_F = det(gp.F);
gp.J = gp.det_dXdxi;
gp.j = gp.det_F * gp.J;
end 
    