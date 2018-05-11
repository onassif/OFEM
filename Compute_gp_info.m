function gp = Compute_gp_info(gp, coor, Umt, Umt_n, iel, num)
%   Compute Shape function derivatives in reference coordinates
gp = ShpFun(gp, coor, iel);

%   Compute deformation gradient F = dx/dX ( calculated inside the object)
gp.U   = Umt  (:, 1:num.ndm)';
gp.U_n = Umt_n(:, 1:num.ndm)';

%   Compute J in both reference and spatial configurations (inside object)
end 
    