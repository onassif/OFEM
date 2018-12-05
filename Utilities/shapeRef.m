function [det_dXdxi_list, dNdX_list, dXdxi_list] = shapeRef(nodes, conn, dNdxi_list)
numel = size(conn    , 1);
nen   = size(conn    , 2);
ndm   = size(dNdxi_list, 2);
ngp   = size(dNdxi_list, 3);


det_dXdxi_list = zeros(ngp, numel);
dNdX_list      = zeros(nen, ndm, ngp, numel);
dXdxi_list     = zeros(ndm, ndm, ngp, numel);

for i = 1:numel
   coor  = nodes(conn(i,:),:)';
   for j = 1:ngp
      dXdxi = coor*dNdxi_list(:,:,j);
      det_dXdxi_list(j,i) = det(dXdxi);
      if ndm == 1
         dNdX_list( :,1,j,i) = dNdxi_list(:,:,j) / dXdxi;
         dXdxi_list(1,1,j,i) = dXdxi;
      else
         dNdX_list( :,:,j,i) = dNdxi_list(:,:,j) / dXdxi;
         dXdxi_list(:,:,j,i) = dXdxi;
      end
   end
end
end

