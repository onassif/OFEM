classdef Q9
   properties
      i
      iel

      b
      eps
      sigma
      ctan
      D
      U
      U_n
      mesh
   end
   
   properties (SetAccess = private)
      finiteDisp
      Nmat
      Ninv
      dNdX_list
      dNdX
      dNdx
      dNdxi_3D
      dNdxi
      det_dXdxi_list
      B
      N
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      xi = [...
         -sqrt(0.6)  -sqrt(0.6)
         +sqrt(0.6)  -sqrt(0.6)
         +sqrt(0.6)  +sqrt(0.6)
         -sqrt(0.6)  +sqrt(0.6)
         0           -sqrt(0.6)
         +sqrt(0.6)  0
         0           +sqrt(0.6)
         -sqrt(0.6)  0
         0           0];
      weights = (1/81).*[25 25 25 25 40 40 40 40 64];
   end
   
   methods
      %% Construct
      function obj = Q9(varargin)
         finiteDisp     = varargin{1};
         obj.finiteDisp = finiteDisp;
         
         if nargin == 3
            obj.xi = varargin{3};
         end
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat( obj);
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);

         if nargin >= 2 && isstruct(varargin{2})
            obj.mesh = varargin{2};
         end
      end
      %% Get functions
      function value = get.dNdxi(obj)
         value = squeeze(obj.dNdxi_3D(:,:,obj.i));
      end

      function value = get.N(obj)
         value = obj.Nmat(obj.i,:);
      end

      function value = get.w(obj)
         value = obj.weights(obj.i);
      end

      function value = get.J(obj)
         value = obj.det_dXdxi_list(obj.iel);
      end

      function value = get.dNdX(obj)
         value = obj.dNdX_list(:,:,obj.i,obj.iel);
      end

      function value = get.F(obj)
         I     = eye(size(obj.U,1));
         value = obj.U*obj.dNdX + I;
      end

      function value = get.j(obj)
         value = det(obj.F) * obj.J;
      end

      function value = get.dNdx(obj)
         value = obj.dNdX / obj.F;
      end

      function value = get.b(obj)
         if obj.finiteDisp
            value = obj.F*obj.F';
         end
      end

      function value = get.B(obj)
         if (obj.finiteDisp)
            dx = obj.dNdx(:,1);
            dy = obj.dNdx(:,2);
         else
            dx = obj.dNdX(:,1);
            dy = obj.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0, dx(5),   0.0, dx(6),   0.0, dx(7),   0.0, dx(8),   0.0, dx(9),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4),   0.0, dy(5),   0.0, dy(6),   0.0, dy(7),   0.0, dy(8),   0.0, dy(9)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4), dy(5), dx(5), dy(6), dx(6), dy(7), dx(7), dy(8), dx(8), dy(9), dx(9)];
      end
      %% Set functions
      function obj = set.U(obj, val)
         if size(val,3)==1 % Normal
            obj.U = val';
         elseif size(val,3) == 2 % DG
            obj.U = permute(val,[2 1 3]);
         end
      end

      function obj = set.U_n(obj, val)
         if size(val,3)==1 % Normal
            obj.U_n = val';
         elseif size(val,3) == 2 % DG
            obj.U_n = permute(val,[2 1 3]);
         end
      end

      function obj = set.mesh(obj, val)
         [obj.det_dXdxi_list, obj.dNdX_list] =...
            obj.computeJ_and_dNdX(val.nodes, val.conn, obj.dNdxi_3D);
      end
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         ngp = size(xi,1);
         ndm = size(xi,2);
         nen = 9;
         dNdxi_3D = zeros(ndm, nen, ngp);

         for i=1:ngp
            x1 = xi(i,1);  x2 = xi(i,2);
            
            dNdxi_3D(:,:,i) =1/4 *[...
               (-1 + x2)*x2*x1 + (-1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*x1 + (-1 + x2)*(-1 + x1)*x1
               (-1 + x2)*x2*x1 + (-1 + x2)*x2*( 1 + x1), x2*( 1 + x1)*x1 + (-1 + x2)*( 1 + x1)*x1
               ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*( 1 + x1), x2*( 1 + x1)*x1 + ( 1 + x2)*( 1 + x1)*x1
               ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*x1 + ( 1 + x2)*(-1 + x1)*x1
               (-1 + x2)*x2*(-2 - 2*x1) - 2*(-1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*(-2 - 2*x1) + (-1 + x2)*(-1 + x1)*(-2 - 2*x1)
               (-1 + x2)*(1 + x2)*(-2 - 2*x1) - 2*(-1 + x2)*(1 + x2)*x1, (1 + x2)*x1*(-2 - 2*x1) + (-1 + x2)*x1*(-2 - 2*x1)
               x2*( 1 + x2)*(-2 - 2*x1) - 2*x2*(1 + x2)*(-1 + x1), (1 + x2)*(-1 + x1)*(-2 - 2*x1) + x2*(-1 + x1)*(-2 - 2*x1)
               -2*(-1 + x2)*( 1 + x2)*x1 - 2*(-1 + x2)*(1 + x2)*(-1 + x1), -2*(1 + x2)*(-1 + x1)*x1 - 2*(-1 + x2)*(-1 + x1)*x1
               (-1 + x2)*(1 + x2)*(4 + 4*x1) + 4*(-1 + x2)*(1 + x2)*(-1 + x1), (1 + x2)*(-1 + x1)*(4 + 4*x1) + (-1 + x2)*(-1 + x1)*(4 + 4*x1)
               ]';
         end
         dNdxi_3D = permute(dNdxi_3D, [2 1 3]);
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         x1 = obj.xi(:,1);  x2 = obj.xi(:,2);
         
         Nmat = 1/4*[...
            x1 .* (x1-1) .*     x2 .* (x2-1),...
            x1 .* (x1+1) .*     x2 .* (x2-1),...
            x1 .* (x1+1) .*     x2 .* (x2+1),...
            x1 .* (x1-1) .*     x2 .* (x2+1),...
            -2 * (x1+1) .* (x1-1) .*     x2 .* (x2-1),...
            -2 * (x1+1) .*     x1 .* (x2+1) .* (x2-1),...
            -2 * (x1+1) .* (x1-1) .* (x2+1) .*     x2,...
            -2 *     x1 .* (x1-1) .* (x2+1) .* (x2-1),...
            +4 * (x1+1) .* (x1-1) .* (x2+1) .* (x2-1)];
         Ninv = inv(Nmat);
      end
      
      function [det_dXdxi_list, dNdX_list] = computeJ_and_dNdX(nodes, conn, dNdxi_3D)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ndm   = size(dNdxi_3D, 2);
         ngp   = size(dNdxi_3D, 3);
         
         det_dXdxi_list = zeros(numel,1);
         dNdX_list      = zeros(nen, ndm, ngp, numel);
         
         for i = 1:numel
            coor  = Q9.removePlane(nodes(conn(i,:),:)');
            dXdxi = coor*dNdxi_3D(:,:,1);
            det_dXdxi_list(i) = det(dXdxi);
            
            for j = 1:ngp
               dNdX_list(:,:,j,i) = dNdxi_3D(:,:,j) / dXdxi;
            end
         end
      end
      
      function fcoor = removePlane(coor)
         if size(coor, 1) == 2
            fcoor = coor;
         elseif size(coor, 1) == 3
            indc  = std(coor,0,2) > 1e-8;
            fcoor = coor(indc,:);
         end
      end
      
   end
end

