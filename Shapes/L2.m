classdef L2
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
      N
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
      xi      = 1/sqrt(3) .*[-1; 1];
      weights = [1 1]';
   end
   
   methods
      %% Construct
      function obj = L2(varargin)
         obj.finiteDisp = varargin{1};
         
         if nargin >= 3
            obj.xi = varargin{3};
         end
         if nargin == 4
            obj.weights = varargin{4};
         end
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(   obj);
         obj.dNdxi_list       = obj.compute_dNdxi(  obj);
         obj.d2Ndxi2_list     = obj.compute_d2Ndxi2(obj);
         
         if nargin >= 2 && isstruct(varargin{2})
            obj.mesh = varargin{2};
         end
      end
      %% Get functions
      function value = get.N(obj)
         value = obj.Nmat(obj.i,:);
      end
      
      function value = get.dNdxi(obj)
         value = obj.dNdxi_list(:,obj.i);
      end
      
      function value = get.d2Ndxi2(obj)
         value = obj.d2Ndxi2_list(:,obj.i);
      end
      
      function value = get.dXdxi(obj)
         value = obj.dXdxi_list(obj.i,obj.iel);
      end
      
      function value = get.w(obj)
         value = obj.weights(obj.i);
      end
      
      function value = get.J(obj)
         value = obj.det_dXdxi_list(obj.i, obj.iel);
      end
      
      function value = get.dNdX(obj)
         value = obj.dNdX_list(:,obj.i,obj.iel);
      end
      
      function value = get.F(obj)
         I     = eye(size(obj.U,1));
         value = obj.U*obj.dNdX + I;
      end
      
      function value = get.JxX(obj)
         value = det(obj.F);
      end
      
      function value = get.j(obj)
         value = det(obj.F) * obj.J;
      end

      function value = get.dNdx(obj)
         value = obj.dNdX / obj.F;
      end

      function value = get.b(obj)
         value = obj.F*obj.F';
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
         [obj.det_dXdxi_list, obj.dNdX_list, obj.dXdxi_list] =...
            obj.computeJ_and_dNdX(val.nodes, val.conn, obj.dNdxi_list);
      end
   end
   
   methods (Static)
      function dNdxi_list = compute_dNdxi(~)
            dNdxi_list = 1/2 *[...
                -1, -1
                +1, +1];
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2(~)
         d2Ndxi2_list = zeros(2);
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
      
      function [det_dXdxi_list, dNdX_list, dXdxi_list] = computeJ_and_dNdX(nodes, conn, dNdxi_list)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ngp   = size(dNdxi_list, 2);
         
         det_dXdxi_list = zeros(ngp, numel);
         dNdX_list      = zeros(nen, ngp, numel);
         dXdxi_list     = zeros(ngp, numel);
         
         for i = 1:numel
            coor  = nodes(conn(i,:),:)';
            for j = 1:ngp
               dXdxi = coor*dNdxi_list(:,j);
               det_dXdxi_list(j,i) = det(dXdxi);
               dNdX_list( :,j,i) = dNdxi_list(:,j) / dXdxi;
               dXdxi_list(:,j,i) = dXdxi;
            end
         end
      end
      
   end
end