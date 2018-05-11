function gp = ShpFun(gp, coor, iel)
gp.iel = iel;

if (sum(sum(gp.dNdX_list(:,:,gp.i,iel))) == 0)
    gp.dXdxi = (gp.dNdxi*coor)'; %[dXdr dYdr; dXds dYds];
    
    gp.det_dXdxi_list(iel) = det(gp.dXdxi);
    
    gp.dNdX_list(:,:,gp.i,iel)  = gp.dNdxi' / gp.dXdxi;
end
end

