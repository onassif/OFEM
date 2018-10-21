classdef Q8
   properties
      i
      iel
      
      b
      eps
      sigma
      D
      U
      U_n
      dU
      mesh
      xi
      weights = [1 1 1 1 1 1 1 1]';
   end
   
   properties (SetAccess = private)
      finiteDisp
      Nmat
      Ninv
      dNdX_list
      dNdX
      dNdx
      dNdxi_list
      d2Ndxi2_list
      dNdxi
      d2Ndxi2
      dXdxi
      dXdxi_list
      det_dXdxi_list
      B
      Bf
      N
      R
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
   end
   
   methods
      %% Construct
      function ob = Q8(varargin)
         ob.finiteDisp = varargin{1};
         
         if nargin >= 3
            ob.xi = varargin{3};
         else
            ob.xi = 1/sqrt(3) .*[...
               -1  1  1 -1 -1  1  1 -1
               -1 -1  1  1 -1 -1  1  1
               -1 -1 -1 -1  1  1  1  1]';
         end
         if nargin == 4
            ob.weights = varargin{4};
         end
         if nargin >= 2 && isstruct(varargin{2})
            ob.mesh = varargin{2};
         end
      end
      %% Get functions
      function value = get.N(ob)
         value = ob.Nmat(ob.i,:);
      end
      
      function value = get.dNdxi(ob)
         value = squeeze(ob.dNdxi_list(:,:,ob.i));
      end
      
      function value = get.d2Ndxi2(ob)
         value = squeeze(ob.d2Ndxi2_list(:,:,ob.i));
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(:,:,ob.i,ob.iel);
      end
      
      function value = get.w(ob)
         value = ob.weights(ob.i);
      end
      
      function value = get.J(ob)
         value = ob.det_dXdxi_list(ob.i,ob.iel);
      end
      
      function value = get.dNdX(ob)
         value = ob.dNdX_list(:,:,ob.i,ob.iel);
      end
      
      function value = get.F(ob)
         I     = eye(size(ob.U,1));
         value = ob.U*ob.dNdX + I;
      end
      
      function value = get.JxX(ob)
         value = det(ob.F);
      end
      
      function value = get.j(ob)
         value = ob.JxX * ob.J;
      end
      
      function value = get.dNdx(ob)
         value = ob.dNdX / ob.F;
      end
      
      function value = get.b(ob)
         value = ob.F*ob.F';
      end
      
      function value = get.B(ob)
         if (ob.finiteDisp)
            dx = ob.dNdx(:,1); dy = ob.dNdx(:,2); dz = ob.dNdx(:,3);
         else
            dx = ob.dNdX(:,1); dy = ob.dNdX(:,2); dz = ob.dNdX(:,3);
         end
         value =[...
            dx(1) 0     0     dx(2) 0     0     dx(3) 0     0     dx(4) 0     0     dx(5) 0     0     dx(6) 0     0     dx(7) 0     0     dx(8) 0     0
            0     dy(1) 0     0     dy(2) 0     0     dy(3) 0     0     dy(4) 0     0     dy(5) 0     0     dy(6) 0     0     dy(7) 0     0     dy(8) 0
            0     0     dz(1) 0     0     dz(2) 0     0     dz(3) 0     0     dz(4) 0     0     dz(5) 0     0     dz(6) 0     0     dz(7) 0     0     dz(8)
            dy(1) dx(1) 0     dy(2) dx(2) 0     dy(3) dx(3) 0     dy(4) dx(4) 0     dy(5) dx(5) 0     dy(6) dx(6) 0     dy(7) dx(7) 0     dy(8) dx(8) 0
            0     dz(1) dy(1) 0     dz(2) dy(2) 0     dz(3) dy(3) 0     dz(4) dy(4) 0     dz(5) dy(5) 0     dz(6) dy(6) 0     dz(7) dy(7) 0     dz(8) dy(8)
            dz(1) 0     dx(1) dz(2) 0     dx(2) dz(3) 0     dx(3) dz(4) 0     dx(4) dz(5) 0     dx(5) dz(6) 0     dx(6) dz(7) 0     dx(7) dz(8) 0     dx(8)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2); dz = ob.dNdx(:,3);
         value =[...
            dx(1)  0      0      dx(2)  0      0      dx(3)  0      0      dx(4)  0      0      dx(5)  0      0      dx(6)  0      0      dx(7)  0      0      dx(8)  0      0
            0      dy(1)  0      0      dy(2)  0      0      dy(3)  0      0      dy(4)  0      0      dy(5)  0      0      dy(6)  0      0      dy(7)  0      0      dy(8)  0
            0      0      dz(1)  0      0      dz(2)  0      0      dz(3)  0      0      dz(4)  0      0      dz(5)  0      0      dz(6)  0      0      dz(7)  0      0      dz(8)
            dy(1)  dx(1)  0      dy(2)  dx(2)  0      dy(3)  dx(3)  0      dy(4)  dx(4)  0      dy(5)  dx(5)  0      dy(6)  dx(6)  0      dy(7)  dx(7)  0      dy(8)  dx(8)  0
            0      dz(1)  dy(1)  0      dz(2)  dy(2)  0      dz(3)  dy(3)  0      dz(4)  dy(4)  0      dz(5)  dy(5)  0      dz(6)  dy(6)  0      dz(7)  dy(7)  0      dz(8)  dy(8)
            dz(1)  0      dx(1)  dz(2)  0      dx(2)  dz(3)  0      dx(3)  dz(4)  0      dx(4)  dz(5)  0      dx(5)  dz(6)  0      dx(6)  dz(7)  0      dx(7)  dz(8)  0      dx(8)
            dy(1) -dx(1)  0      dy(2) -dx(2)  0      dy(3) -dx(3)  0      dy(4) -dx(4)  0      dy(5) -dx(5)  0      dy(6) -dx(6)  0      dy(7) -dx(7)  0      dy(8) -dx(8)  0
            0      dz(1) -dy(1)  0      dz(2) -dy(2)  0      dz(3) -dy(3)  0      dz(4) -dy(4)  0      dz(5) -dy(5)  0      dz(6) -dy(6)  0      dz(7) -dy(7)  0      dz(8) -dy(8)
            -dz(1)  0      dx(1) -dz(2)  0      dx(2) -dz(3)  0      dx(3) -dz(4)  0      dx(4) -dz(5)  0      dx(5) -dz(6)  0      dx(6) -dz(7)  0      dx(7) -dz(8)  0      dx(8)];
      end
      
      function value = get.R(ob)
         [P, ~, Q] = svd(ob.F);
         value =  P*Q';
      end
      
      %% Set functions
      function ob = set.mesh(ob, val)
         ob.mesh = val;
         [ob.det_dXdxi_list, ob.dNdX_list, ob.dXdxi_list] = ...
            ob.computeJ_and_dNdX(val.nodes, val.conn, ob.dNdxi_list);
      end
      
      function ob = set.xi(ob, val)
         ob.xi = val;
         [ob.Nmat, ob.Ninv] = ob.compute_Nmat(   val);
         ob.dNdxi_list      = ob.compute_dNdxi(  val);
         ob.d2Ndxi2_list    = ob.compute_d2Ndxi2(val);
      end
   end
   
   methods (Static)
      function [Nmat, Ninv] = compute_Nmat(xi)
         Nmat = 1/8*[...
            (1-xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3)),...
            (1+xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3)),...
            (1+xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3)),...
            (1-xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3)),...
            (1-xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3)),...
            (1+xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3)),...
            (1+xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3)),...
            (1-xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3))];
         if size(Nmat,1) == size(Nmat,2)
            Ninv = inv(Nmat);
         else
            Ninv = 0;
         end
      end
      
      function dNdxi_list = compute_dNdxi(xi)
         ngp = size(xi,1);
         ndm = size(xi,2);
         nen = 8;
         dNdxi_list = zeros(ndm, nen, ngp);
         for i=1:size(xi,1)
            dNdxi_list(:,:,i) =1/8*[...
               -(xi(i,2)- 1)*(xi(i,3)- 1), -(xi(i,1)- 1)*(xi(i,3)- 1), -(xi(i,1)- 1)*(xi(i,2)- 1)
               +(xi(i,2)- 1)*(xi(i,3)- 1),  (xi(i,1)+ 1)*(xi(i,3)- 1),  (xi(i,1)+ 1)*(xi(i,2)- 1)
               -(xi(i,2)+ 1)*(xi(i,3)- 1), -(xi(i,1)+ 1)*(xi(i,3)- 1), -(xi(i,1)+ 1)*(xi(i,2)+ 1)
               +(xi(i,2)+ 1)*(xi(i,3)- 1),  (xi(i,1)- 1)*(xi(i,3)- 1),  (xi(i,1)- 1)*(xi(i,2)+ 1)
               +(xi(i,2)- 1)*(xi(i,3)+ 1),  (xi(i,1)- 1)*(xi(i,3)+ 1),  (xi(i,1)- 1)*(xi(i,2)- 1)
               -(xi(i,2)- 1)*(xi(i,3)+ 1), -(xi(i,1)+ 1)*(xi(i,3)+ 1), -(xi(i,1)+ 1)*(xi(i,2)- 1)
               +(xi(i,2)+ 1)*(xi(i,3)+ 1),  (xi(i,1)+ 1)*(xi(i,3)+ 1),  (xi(i,1)+ 1)*(xi(i,2)+ 1)
               -(xi(i,2)+ 1)*(xi(i,3)+ 1), -(xi(i,1)- 1)*(xi(i,3)+ 1), -(xi(i,1)- 1)*(xi(i,2)+ 1)]';
         end
         dNdxi_list = permute(dNdxi_list, [2 1 3]);
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2(xi)
         ngp = size(xi,1);
         nen = 8;
         d2Ndxi2_list = zeros(6, nen, ngp);
         for i=1:size(xi,1)
            d2Ndxi2_list(:,:,i) =1/8*[...
               0, 0, 0, -xi(i,3)+1, -xi(i,1)+1, -xi(i,2)+1
               0, 0, 0,  xi(i,3)-1,  xi(i,1)+1,  xi(i,2)-1
               0, 0, 0, -xi(i,3)+1, -xi(i,1)-1, -xi(i,2)-1
               0, 0, 0,  xi(i,3)-1,  xi(i,1)-1,  xi(i,2)+1
               0, 0, 0,  xi(i,3)+1,  xi(i,1)-1,  xi(i,2)-1
               0, 0, 0, -xi(i,3)-1, -xi(i,1)-1, -xi(i,2)+1
               0, 0, 0,  xi(i,3)+1,  xi(i,1)+1,  xi(i,2)+1
               0, 0, 0, -xi(i,3)-1, -xi(i,1)+1, -xi(i,2)-1]';
         end
         d2Ndxi2_list = permute(d2Ndxi2_list, [2 1 3]);
      end
      
      function [det_dXdxi_list, dNdX_list, dXdxi_list] = computeJ_and_dNdX(nodes, conn, dNdxi_list)
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
               dNdX_list( :,:,j,i) = dNdxi_list(:,:,j) / dXdxi;
               dXdxi_list(:,:,j,i) = dXdxi;
            end
         end
      end
      
   end
end