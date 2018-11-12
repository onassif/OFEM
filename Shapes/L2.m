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
      xi      
      weights = [1 1]';
      Nmat
      dNdxi_list
      d2Ndxi2_list
      
      dNdX_list
      dXdxi_list
      det_dXdxi_list
   end
   
   properties (SetAccess = private)
      finiteDisp
      Ninv
      dNdX
      dNdx
      dNdxi
      d2Ndxi2
      dXdxi
      B
      N
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
   end
   
   methods
      %% Construct
      function ob = L2(varargin)
         ob.finiteDisp = varargin{1};
         
         if nargin >= 3
            ob.xi = varargin{3};
         else
            ob.xi = 1/sqrt(3) .*[-1; 1];
         end
         if nargin == 4
            ob.weights = varargin{4};
         end
      end
      %% Get functions
      function value = get.N(ob)
         value = ob.Nmat(ob.i,:);
      end
      
      function value = get.dNdxi(ob)
         value = ob.dNdxi_list(:,ob.i);
      end
      
      function value = get.d2Ndxi2(ob)
         value = ob.d2Ndxi2_list(:,ob.i);
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(ob.i,ob.iel);
      end
      
      function value = get.w(ob)
         value = ob.weights(ob.i);
      end
      
      function value = get.J(ob)
         value = ob.det_dXdxi_list(ob.i, ob.iel);
      end
      
      function value = get.dNdX(ob)
         value = ob.dNdX_list(:,ob.i,ob.iel);
      end
      
      function value = get.F(ob)
         I     = eye(size(ob.U,1));
         value = ob.U*ob.dNdX + I;
      end
      
      function value = get.JxX(ob)
         value = det(ob.F);
      end
      
      function value = get.j(ob)
         value = det(ob.F) * ob.J;
      end

      function value = get.dNdx(ob)
         value = ob.dNdX / ob.F;
      end

      function value = get.b(ob)
         value = ob.F*ob.F';
      end
      
      function value = get.B(ob)
         if (ob.finiteDisp)
            dx = ob.dNdx;
         else
            dx = ob.dNdX;
         end
         value=[dx(1), dx(2)];
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
      function [Nmat, Ninv] = compute_Nmat(xi)
         Nmat = 1/2*[(1-xi(:,1)), (1+xi(:,1))];
         if size(Nmat,1) == size(Nmat,2)
            Ninv = inv(Nmat);
         else
            Ninv = 0;
         end
      end
      
      function dNdxi_list = compute_dNdxi()
         dNdxi_list(:,1,:) = 1/2 *[-1, -1; +1, +1];
      end
      
      function d2Ndxi2_list = compute_d2Ndxi2()
         d2Ndxi2_list(:,1,:) = zeros(2);
      end      
   end
end