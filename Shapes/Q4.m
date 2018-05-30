classdef Q4
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
      xi = 1/sqrt(3) .*[...
         -1 -1
         +1 -1
         +1 +1
         -1 +1];
      weights = [1 1 1 1];
   end
   
   methods
      function obj = Q4(varargin)
         num        = varargin{1};
         finiteDisp = varargin{2};
         if nargin == 3
            obj.xi = varargin{3};
         end
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
         
         obj.det_dXdxi_list = zeros(num.el,1);
         obj.dNdX_list      = zeros(num.nen, num.ndm, num.gp, num.el);
         
         obj.finiteDisp = finiteDisp;
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
      
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         dNdxi_3D = zeros(2,4,size(xi,1));
         
         dNdxi_3D(1,:,:) = 1/4*[...
            -(1-xi(:,2)), ( 1-xi(:,2)), ( 1+xi(:,2)), -(1+xi(:,2))]';
         
         dNdxi_3D(2,:,:) = 1/4*[...
            -(1-xi(:,1)), -(1+xi(:,1)), ( 1+xi(:,1)), ( 1-xi(:,1))]';         
            

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
   end
end

