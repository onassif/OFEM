classdef Q9
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
      xi = [...
         -sqrt(0.6)  -sqrt(0.6)
         +sqrt(0.6)  -sqrt(0.6)
         +sqrt(0.6)  +sqrt(0.6)
         -sqrt(0.6)  +sqrt(0.6)
         0           -sqrt(0.6)
         +sqrt(0.6)  0
         0           +sqrt(0.6)
         -sqrt(0.6)  0
         0           0]; 
      weights = (1/81).*[25 25 25 25 40 40 40 40 64];
   end
   
   methods
      %% Construct
      function obj = Q9(varargin)
         num        = varargin{1};
         finiteDisp = varargin{2};
         if nargin == 3
            obj.xi = varargin{3};
         end
         obj.dNdxi_3D           = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv]   = obj.compute_Nmat(obj);
         obj.det_dXdxi_list     = zeros(num.el,1);
         obj.dNdX_list          = zeros(num.nen, num.ndm, num.gp, num.el);
         
         obj.finiteDisp         = finiteDisp;
      end
      %% Get functions
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
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0, dx(5),   0.0, dx(6),   0.0, dx(7),   0.0, dx(8),   0.0, dx(9),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4),   0.0, dy(5),   0.0, dy(6),   0.0, dy(7),   0.0, dy(8),   0.0, dy(9) 
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4), dy(5), dx(5), dy(6), dx(6), dy(7), dx(7), dy(8), dx(8), dy(9), dx(9)];
      end
      %% Set functions
      function obj = set.U(obj, value)
         if size(value,3)==1 % Normal
            obj.U = value';
         elseif size(value,3) == 2 % DG
            obj.U = permute(value,[2 1 3]);
         end
      end
      
      function obj = set.U_n(obj, value)
         if size(value,3)==1 % Normal
            obj.U_n = value';
         elseif size(value,3) == 2 % DG
            obj.U_n = value';
         end
      end
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         dNdxi_3D = zeros(2,9,9);
         for i=1:9
             x1 = xi(i,1);  x2 = xi(i,2);
             
             dNdxi_3D(:,:,i) =1/4 *[...
             (-1 + x2)*x2*x1 + (-1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*x1 + (-1 + x2)*(-1 + x1)*x1
             (-1 + x2)*x2*x1 + (-1 + x2)*x2*( 1 + x1), x2*( 1 + x1)*x1 + (-1 + x2)*( 1 + x1)*x1
             ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*( 1 + x1), x2*( 1 + x1)*x1 + ( 1 + x2)*( 1 + x1)*x1
             ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*x1 + ( 1 + x2)*(-1 + x1)*x1
             (-1 + x2)*x2*(-2 - 2*x1) - 2*(-1 + x2)*x2*(-1 + x1), x2*(-1 + x1)*(-2 - 2*x1) + (-1 + x2)*(-1 + x1)*(-2 - 2*x1)
             (-1 + x2)*(1 + x2)*(-2 - 2*x1) - 2*(-1 + x2)*(1 + x2)*x1, (1 + x2)*x1*(-2 - 2*x1) + (-1 + x2)*x1*(-2 - 2*x1)
             x2*( 1 + x2)*(-2 - 2*x1) - 2*x2*(1 + x2)*(-1 + x1), (1 + x2)*(-1 + x1)*(-2 - 2*x1) + x2*(-1 + x1)*(-2 - 2*x1)
             -2*(-1 + x2)*( 1 + x2)*x1 - 2*(-1 + x2)*(1 + x2)*(-1 + x1), -2*(1 + x2)*(-1 + x1)*x1 - 2*(-1 + x2)*(-1 + x1)*x1
             (-1 + x2)*(1 + x2)*(4 + 4*x1) + 4*(-1 + x2)*(1 + x2)*(-1 + x1), (1 + x2)*(-1 + x1)*(4 + 4*x1) + (-1 + x2)*(-1 + x1)*(4 + 4*x1)           
             ]';
             
%                          dNdxi_3D2(:,:,i) =...
%                 1/4 *[...
%                +x2 - x2^2 + 2*x1*(1-x2) - 2*x1*(1-x2^2), +x1 - x1^2  + 2*x2*(1-x1) - 2*x2*(1-x1^2)
%                -x2 + x2^2 + 2*x1*(1-x2) - 2*x1*(1-x2^2), -x1 - x1^2  + 2*x2*(1+x1) - 2*x2*(1-x1^2)
%                +x2 + x2^2 + 2*x1*(1+x2) - 2*x1*(1-x2^2), +x1 + x1^2  + 2*x2*(1+x1) - 2*x2*(1-x1^2)
%                -x2 - x2^2 + 2*x1*(1+x2) - 2*x1*(1-x2^2), -x1 + x1^2  + 2*x2*(1-x1) - 2*x2*(1-x1^2)
%                -4*x1*(1-x2) + 4*x1*(1-x2^2), -2 + 2*x1^2       + 4*x2*(1-x1^2)
%                +2 - 2*x2^2       + 4*x1*(1-x2^2), -4*x2*(1+x1) + 4*x2*(1-x1^2)
%                -4*x1*(1+x2) + 4*x1*(1-x2^2), +2 - 2*x1^2       + 4*x2*(1-x1^2)
%                -2 + 2*x2^2       + 4*x1*(1-x2^2), -4*x2*(1-x1) + 4*x2*(1-x1^2)
%                -8*x1*(1-x2^2), -8*x2*(1-x1^2)
%                ]';
         end
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         x1 = obj.xi(:,1);  x2 = obj.xi(:,2);
         
         Nmat = 1/4*[...
                      x1 .* (x1-1) .*     x2 .* (x2-1),...
                      x1 .* (x1+1) .*     x2 .* (x2-1),...
                      x1 .* (x1+1) .*     x2 .* (x2+1),...
                      x1 .* (x1-1) .*     x2 .* (x2+1),...
             -2 * (x1+1) .* (x1-1) .*     x2 .* (x2-1),...
             -2 * (x1+1) .*     x1 .* (x2+1) .* (x2-1),...
             -2 * (x1+1) .* (x1-1) .* (x2+1) .*     x2,...
             -2 *     x1 .* (x1-1) .* (x2+1) .* (x2-1),...
             +4 * (x1+1) .* (x1-1) .* (x2+1) .* (x2-1)];
         Ninv = inv(Nmat);
      end
   end
end

