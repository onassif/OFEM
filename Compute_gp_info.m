function gp = Compute_gp_info(gp, coor, iel)
%   Compute Shape function derivatives in reference coordinates
gp = ShpFun(gp, coor, iel);
%   Compute J in both reference and spatial configurations (inside object)
end 
    