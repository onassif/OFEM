classdef L2
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
      xi      = 1/sqrt(3) .*[-1; 1];
      weights = [1 1]';
   end
   
   methods
      %% Construct
      function obj = L2(varargin)
         finiteDisp     = varargin{1};
         obj.finiteDisp = finiteDisp;
         
         if nargin >= 3
            obj.xi = varargin{3};
         end
         if nargin == 4
            obj.weights = varargin{4};
         end
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat( obj);
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);

         if nargin >= 2 && isstruct(varargin{2})
            obj.mesh = varargin{2};
         end
      end
      %% Get functions
      function value = get.dNdxi(obj)
         value = obj.dNdxi_3D(:,obj.i);
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
         value = obj.dNdX_list(:,obj.i,obj.iel);
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
            dx = obj.dNdx;
         else
            dx = obj.dNdX;
         end
         value=[dx(1), dx(2)];
      end
      %% Set functions
      function obj = set.mesh(obj, val)
         [obj.det_dXdxi_list, obj.dNdX_list] =...
            obj.computeJ_and_dNdX(val.nodes, val.conn, obj.dNdxi_3D);
      end
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(~)
            dNdxi_3D = 1/2 *[...
                -1, -1
                +1, +1];
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         xi = obj.xi;
         
         Nmat = 1/2*[(1-xi(:,1)), (1+xi(:,1))];
         if size(Nmat,1) == size(Nmat,2)
            Ninv = inv(Nmat);
         else
            Ninv = 0;
         end
      end
      
      function [det_dXdxi_list, dNdX_list] = computeJ_and_dNdX(nodes, conn, dNdxi_3D)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ngp   = size(dNdxi_3D, 2);
         
         det_dXdxi_list = zeros(numel,1);
         dNdX_list      = zeros(nen, ngp, numel);
         
         for i = 1:numel
            coor  = nodes(conn(i,:),:)';
            
            dXdxi = coor*dNdxi_3D(:,1:size(coor,1));
            det_dXdxi_list(i) = det(dXdxi);
            
            for j = 1:ngp
               dNdX_list(:,j,i) = dNdxi_3D(:,j) / dXdxi;
            end
         end
      end
      
   end
end