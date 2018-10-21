classdef T3
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
      weights = 0.5;
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
      function ob = T3(varargin)
         finiteDisp     = varargin{1};
         ob.finiteDisp = finiteDisp;
         
         if nargin >= 3
            ob.xi = varargin{3};
         else
            ob.xi = 1/3 .*[1 1];
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
         if size(ob.xi, 1)>1
            value = ob.Nmat(ob.i,:);
         else
            value = ob.Nmat;
         end
      end
      
      function value = get.dNdxi(ob)
         value = ob.dNdxi_list;
      end

      function value = get.d2Ndxi2(ob)
         value = ob.d2Ndxi2_list;
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(:,:,ob.iel);
      end
      
      function value = get.w(ob)
         value = ob.weights(ob.i);
      end
      
      function value = get.J(ob)
         value = ob.det_dXdxi_list(ob.iel);
      end
      
      function value = get.dNdX(ob)
         value = ob.dNdX_list(:,:,ob.iel);
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
            dx = ob.dNdx(:,1);  dy = ob.dNdx(:,2);
         else
            dx = ob.dNdX(:,1);  dy = ob.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2);
         value =[...
            dx(1)  0      dx(2)  0      dx(3)  0    
            0      dy(1)  0      dy(2)  0      dy(3)
            dy(1)  dx(1)  dy(2)  dx(2)  dy(3)  dx(3)
            dy(1) -dx(1)  dy(2) -dx(2)  dy(3) -dx(3)];
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
         [ob.Nmat, ob.Ninv] = ob.compute_Nmat(val);
         ob.dNdxi_list      = ob.compute_dNdxi(  );
         ob.d2Ndxi2_list    = ob.compute_d2Ndxi2();
      end
   end
   
   methods (Static)
      function [Nmat, Ninv] = compute_Nmat(ob)
         xi = ob.xi;
         Nmat = [1-xi(:,1)-xi(:,2), xi(:,1), xi(:,2)];
         Ninv = ((Nmat*Nmat')\Nmat)';
      end
      
      function dNdxi_list = compute_dNdxi()
         dNdxi_list =[...
            -1 -1
            +1  0
            +0  1];
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2()
         d2Ndxi2_list =[...
            0 0 0
            0 0 0
            0 0 0];
      end
   
      function [det_dXdxi_list, dNdX_list, dXdxi_list] = computeJ_and_dNdX(nodes, conn, dNdxi_list)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ndm   = size(dNdxi_list, 2);
         ngp   = size(dNdxi_list, 3);
         
         det_dXdxi_list = zeros(numel,1);
         dNdX_list      = zeros(nen, ndm, numel);
         dXdxi_list     = zeros(ndm, ndm, numel);
         
         for i = 1:numel
            coor  = nodes(conn(i,:),:)';
            dXdxi = coor*dNdxi_list;
            det_dXdxi_list(i) = det(dXdxi);
            
            dNdX_list(:,:,i) = dNdxi_list / dXdxi;
            dXdxi_list(:,:,i) = dXdxi;
         end
      end
      
   end
end

