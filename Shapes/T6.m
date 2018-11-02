classdef T6
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
      weights = 1/6* [1 1 1]';
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
      function ob = T6(varargin)
         finiteDisp     = varargin{1};
         ob.finiteDisp = finiteDisp;
         
         if nargin >= 3
            ob.xi = varargin{3};
         else
            ob.xi = 1/6*[...
               4 1 1
               1 1 4]';
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
            dx = ob.dNdx(:,1);	dy = ob.dNdx(:,2);
         else
            dx = ob.dNdX(:,1);	dy = ob.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0, dx(5),   0.0, dx(6),   0.0
            0.000, dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4),   0.0, dy(5),   0.0, dy(6)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4), dy(5), dx(5), dy(6), dx(6)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2);
         value =[...
            dx(1)  0      dx(2)  0      dx(3)  0      dx(4)  0      dx(5)  0      dx(6)  0
            0      dy(1)  0      dy(2)  0      dy(3)  0      dy(4)  0      dy(5)  0      dy(6)
            dy(1)  dx(1)  dy(2)  dx(2)  dy(3)  dx(3)  dy(4)  dx(4)  dy(5)  dx(5)  dy(6)  dx(6)
            dy(1) -dx(1)  dy(2) -dx(2)  dy(3) -dx(3)  dy(4) -dx(4)  dy(5) -dx(5)  dy(6) -dx(6)];
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
         x1 = xi(:,1);  x2 = xi(:,2);
         
         Nmat = [...
            (1-x1-x2).*(2*(1-x1-x2)-1),...
            x1.*(2*x1-1),...
            x2.*(2*x2-1),...
            4*(1-x1-x2).*x1,...
            4*x1.*x2,...
            4*x2.*(1-x1-x2)];
         
         Ninv = ((Nmat*Nmat')\Nmat)';
      end
      
      function dNdxi_list = compute_dNdxi(xi)
         ngp = size(xi,1);
         ndm = size(xi,2);
         nen = 6;
         dNdxi_list = zeros(nen, ndm, ngp);
         
         for i=1:ngp
            x1 = xi(i,1);  x2 = xi(i,2);
            
            dNdxi_list(:,:,i) =[...
               -3 + 4*(x1+x2), -3 + 4*(x1+x2)
               4*x1-1        ,  0
               0             ,  4*x2 - 1
               4-8*x1-4*x2   , -4*x1
               4*x2          ,  4*x1
               -4*x2         , -8*x2 + 4 - 4*x1];
         end
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2(xi)
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