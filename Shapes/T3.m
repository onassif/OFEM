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
      Nmat
      dNdxi_list
      d2Ndxi2_list
      
      dNdX_list
      dXdxi_list
      det_dXdxi_list
   end
   
   properties (SetAccess = private)
      finiteDisp
      bubb
      bubbB
      Ninv
      dNdX
      dNdx
      dNdxi
      d2Ndxi2
      dXdxi
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
      
      function val = get.bubb(ob)
         r = ob.xi(:,1); s = ob.xi(:,2);
         val = 4*(1-r-s).*r;
      end
      
      function val = get.bubbB(ob)
         r = ob.xi(ob.i,1); s = ob.xi(ob.i,2);
         dbdxi = 4*[(1-2*r-s), -r];
         if (ob.finiteDisp)
            dxdxi = ob.F*ob.dXdxi;
            dbdx  = dbdxi / dxdxi;
            val = [...
               dbdx(1) 0       dbdx(2)  dbdx(2)
               0       dbdx(2) dbdx(1) -dbdx(1)]';
         else
            dbdX  = dbdxi / ob.dXdxi;
            val = [...
               dbdX(1) 0       dbdX(2)
               0       dbdX(2) dbdX(1)]';
         end
      end
      
      %% Set functions      
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
      
   end
end