function gp = ShpFun(gp, coor, iel)
gp.iel = iel;

if (gp.det_dXdxi_list(gp.i,iel) == 0)
    dXdxi = gp.dNdxi*coor; %[dXdr dYdr; dXds dYds];
    
    gp.det_dXdxi_list(gp.i,iel) = det(dXdxi);
    
    gp.dNdX_list(:,:,gp.i,iel)  = gp.dNdxi' / dXdxi;
end
end

