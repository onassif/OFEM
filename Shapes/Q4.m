classdef Q4
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
      xi = 1/sqrt(3) .*[...
         -1 +1 +1 -1
         -1 -1 +1 +1]';
      weights = [1 1 1 1]';
   end
   
   methods
      %% Construct
      function obj = Q4(varargin)
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
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4)];
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
         nen = 4;
         dNdxi_3D = zeros(ndm, nen, ngp);
         
         dNdxi_3D(1,:,:) = 1/4*[...
            -(1-xi(:,2)), ( 1-xi(:,2)), ( 1+xi(:,2)), -(1+xi(:,2))]';
         
         dNdxi_3D(2,:,:) = 1/4*[...
            -(1-xi(:,1)), -(1+xi(:,1)), ( 1+xi(:,1)), ( 1-xi(:,1))]';
         
         dNdxi_3D = permute(dNdxi_3D, [2 1 3]);
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         xi = obj.xi;
         Nmat = 1/4*[...
            (1-xi(:,1)).*(1-xi(:,2)),...
            (1+xi(:,1)).*(1-xi(:,2)),...
            (1+xi(:,1)).*(1+xi(:,2)),...
            (1-xi(:,1)).*(1+xi(:,2))];
         if size(Nmat,1) == size(Nmat,2)
            Ninv = inv(Nmat);
         else
            Ninv = 0;
         end
      end
      
      function [det_dXdxi_list, dNdX_list] = computeJ_and_dNdX(nodes, conn, dNdxi_3D)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ndm   = size(dNdxi_3D, 2);
         ngp   = size(dNdxi_3D, 3);
         
         det_dXdxi_list = zeros(numel,1);
         dNdX_list      = zeros(nen, ndm, ngp, numel);
         
         for i = 1:numel
            coor  = Q4.removePlane(nodes(conn(i,:),:)');
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

