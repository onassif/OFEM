classdef L2
   
   properties (Constant)
      % Isoparametric gp
      xi = 1/sqrt(3) .*[-1; 1];
      weights = [1 1];
   end
   
   properties (SetAccess = public, GetAccess = public)
      i
      dXdxi;
      dxdX;
      dNdx;
      b;
      eps;
      sigma;
      ctan;
      D;
      U;
      U_n;
   end
   
   properties (Hidden)
      iel;
      det_dXdxi_list;
      dNdX_list;
   end
   
   properties (Hidden, SetAccess = private)
      dNdxi_3D;
      numel;
      finiteDisp;
      adof;
   end
   
   properties (SetAccess = private)
      Nmat;
      Ninv;
      dNdX;
      B;
      N;
      w;
      dNdxi;
      F;
      J; % det(dX/dxi) = J
      j; % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
   end
   
   methods
      function obj = L2(num, finiteDisp)
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
         obj.det_dXdxi_list   = zeros(num.el, 1);
         obj.dNdX_list        = zeros(2, 2, num.el);
         obj.finiteDisp       = finiteDisp;
      end
      
      function value = get.dNdxi(obj)
         value = obj.dNdxi_3D(obj.i,:);
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
      
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(~)
            dNdxi_3D = 1/2 *[...
                -1, +1
                -1, +1];
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         xi = obj.xi;
         Nmat = 1/2*[(1-xi), (1+xi)];
         Ninv = inv(Nmat);
      end
   end
end