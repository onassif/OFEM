classdef T6
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
      Bf
      N
      R
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
      xi = 1/6*[...
         4 1 1
         1 1 4]';
      
      weights = 1/6* [1 1 1]';
   end
   
   methods
      %% Construct
      function obj = T6(varargin)
         finiteDisp     = varargin{1};
         obj.finiteDisp = finiteDisp;
         
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
         value = squeeze(obj.dNdxi_list(:,:,obj.i));
      end
      
      function value = get.dXdxi(obj)
         value = obj.dXdxi_list(:,:,obj.i,obj.iel);
      end
      
      function value = get.w(obj)
         value = obj.weights(obj.i);
      end
      
      function value = get.J(obj)
         value = obj.det_dXdxi_list(obj.i,obj.iel);
      end
      
      function value = get.dNdX(obj)
         value = obj.dNdX_list(:,:,obj.i,obj.iel);
      end
      
      function value = get.F(obj)
         I     = eye(size(obj.U,1));
         value = obj.U*obj.dNdX + I;
      end
      
      function value = get.JxX(obj)
         value = det(obj.F);
      end
      
      function value = get.j(obj)
         value = obj.JxX * obj.J;
      end
      
      function value = get.dNdx(obj)
         value = obj.dNdX / obj.F;
      end
      
      function value = get.b(obj)
         value = obj.F*obj.F';
      end
      
      function value = get.B(obj)
         if (obj.finiteDisp)
            dx = obj.dNdx(:,1);	dy = obj.dNdx(:,2);
         else
            dx = obj.dNdX(:,1);	dy = obj.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0, dx(5),   0.0, dx(6),   0.0
            0.000, dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4),   0.0, dy(5),   0.0, dy(6)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4), dy(5), dx(5), dy(6), dx(6)];
      end
      
      function value = get.Bf(obj)
         dx = obj.dNdx(:,1); dy = obj.dNdx(:,2);
         value =[...
            dx(1)  0      dx(2)  0      dx(3)  0      dx(4)  0      dx(5)  0      dx(6)  0    
            0      dy(1)  0      dy(2)  0      dy(3)  0      dy(4)  0      dy(5)  0      dy(6)
            dy(1)  dx(1)  dy(2)  dx(2)  dy(3)  dx(3)  dy(4)  dx(4)  dy(5)  dx(5)  dy(6)  dx(6)
            dy(1) -dx(1)  dy(2) -dx(2)  dy(3) -dx(3)  dy(4) -dx(4)  dy(5) -dx(5)  dy(6) -dx(6)];
      end
      
      function value = get.R(obj)
           [P, ~, Q] = svd(obj.F);
           value =  P*Q';
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
      
      function obj = set.dU(obj, val)
         if size(val,3)==1 % Normal
            obj.dU = val';
         elseif size(val,3) == 2 % DG
            obj.dU = permute(val,[2 1 3]);
         end
      end
      
      function obj = set.mesh(obj, val)
         [obj.det_dXdxi_list, obj.dNdX_list, obj.dXdxi_list] =...
            obj.computeJ_and_dNdX(val.nodes, val.conn, obj.dNdxi_list);
      end
   end
   
   methods (Static)
      function dNdxi_list = compute_dNdxi(obj)
         xi = obj.xi;
         ngp = size(xi,1);
         ndm = size(xi,2);
         nen = 6;
         dNdxi_list = zeros(nen, ndm, ngp);
         
         for i=1:ngp
            x1 = obj.xi(i,1);  x2 = obj.xi(i,2);
            
            dNdxi_list(:,:,i) =[...
               -3 + 4*(x1+x2), -3 + 4*(x1+x2)
               4*x1-1        ,  0
               0             ,  4*x2 - 1
               4-8*x1-4*x2   , -4*x1
               4*x2          ,  4*x1
               -4*x2         , -8*x2 + 4 - 4*x1];
         end
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2(obj)
         xi = obj.xi;
         ngp = size(xi,1);
         nen = 6;
         d2Ndxi2_list = zeros(nen, 3, ngp);
         
         for i=1:ngp
            d2Ndxi2_list(:,:,i) =[...
               +4,  4,  4
               +4,  0,  0
               +0,  4,  0
               -8,  0, -4
               +0,  0,  4
               +0, -8, -4];
         end
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         x1 = obj.xi(:,1);  x2 = obj.xi(:,2);
         
         Nmat = [...
            (1-x1-x2).*(2*(1-x1-x2)-1),...
            x1.*(2*x1-1),...
            x2.*(2*x2-1),...
            4*(1-x1-x2).*x1,...
            4*x1.*x2,...
            4*x2.*(1-x1-x2)];

         Ninv = ((Nmat*Nmat')\Nmat)';
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