function gp = Compute_gp_info(gp, coor, Umt, iel, num)
ndm = num.ndm;
%   Compute Shape function derivatives in reference coordinates
gp = ShpFun(gp, coor, iel);

%   Compute deformation gradient F = dx/dX 
gp.F = Umt(:, 1:num.ndm)'*gp.dNdX + eye(ndm); % dN/dx is calculated(if needed) inside the object

%   Compute J in both reference and spatial configurations
gp.det_F = det(gp.F);
gp.J = gp.det_dXdxi; %                            det(dX/dxi) = J
gp.j = gp.det_F * gp.J; %or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
end 
    