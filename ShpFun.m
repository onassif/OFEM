function gp = ShpFun(gp, coor, iel)
gp.iel = iel;

% if size(gp.dNdX_list,4)>=iel
    gp.dXdxi = (gp.dNdxi*coor)'; %[dXdr dYdr; dXds dYds];
    
%     gp.det_dXdxi_list(iel) = det(gp.dXdxi);
    
%     gp.dNdX_list(:,:,gp.i,iel)  = gp.dNdxi' / gp.dXdxi;
% end
end

