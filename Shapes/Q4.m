classdef Q4
   
   properties (Constant)
      % Isoparametric gp
      xi = 1/sqrt(3) .*[...
         -1  1  1 -1
         -1 -1  1  1]';
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
   end
   
   properties (Hidden)
      iel;
      det_dXdxi_list;
      dNdX_list;
   end
   
   properties (Hidden, SetAccess = private)
      dNdxi_3D;
      numel;
      weights = [1 1 1 1];
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
      function obj = Q4(num, finiteDisp)
         obj.dNdxi_3D = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
         obj.det_dXdxi_list = zeros(4,num.el);
         obj.dNdX_list = zeros(4,2,4,num.el);
         obj.finiteDisp = finiteDisp;
         %             obj.adof = num.ndof - num.ndm;
         obj.adof = 0;
      end
      
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
         value = obj.det_dXdxi_list(obj.i,obj.iel);
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
         a = zeros(1,obj.adof);
         value=[...
            dx(1),   0.0, a, dx(2),   0.0, a, dx(3),   0.0, a, dx(4),   0.0 a
            0.0  , dy(1), a,   0.0, dy(2), a,   0.0, dy(3), a,   0.0, dy(4) a
            dy(1), dx(1), a, dy(2), dx(2), a, dy(3), dx(3), a, dy(4), dx(4) a];
      end
      
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         dNdxi_3D = zeros(2,4,4);
         for i=1:4
            dNdxi_3D(:,:,i) =...
               1/4 *[...
               -(1-xi(i,2)) -(1-xi(i,1))
               ( 1-xi(i,2)) -(1+xi(i,1))
               ( 1+xi(i,2))  (1+xi(i,1))
               -(1+xi(i,2))  (1-xi(i,1))]';
         end
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         xi = obj.xi;
         Nmat = 1/4*[...
            (1-xi(:,1)).*(1-xi(:,2)), (1+xi(:,1)).*(1-xi(:,2)),...
            (1+xi(:,1)).*(1+xi(:,2)), (1-xi(:,1)).*(1+xi(:,2))];
         Ninv = inv(Nmat);
      end
   end
end

