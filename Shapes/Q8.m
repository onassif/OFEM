classdef Q8
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
      mesh;
   end
   
   properties
      iel;
      det_dXdxi_list;
      dNdX_list;
   end
   
   properties (SetAccess = private)
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
         -1  1  1 -1 -1  1  1 -1
         -1 -1  1  1 -1 -1  1  1
         -1 -1 -1 -1  1  1  1  1]';
      weights = [1 1 1 1 1 1 1 1];
   end
   
   methods
      %% Consturct
      function obj = Q8(varargin)
         finiteDisp     = varargin{1};
         obj.finiteDisp = finiteDisp;
         
         if nargin == 3
            obj.xi = varargin{3};
         end
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat( obj);
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);

         if nargin >= 2
            obj.mesh = varargin{2};
         end
      end
      %% Get functions
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
            dx = obj.dNdx(:,1); dy = obj.dNdx(:,2); dz = obj.dNdx(:,3);
         else
            dx = obj.dNdX(:,1); dy = obj.dNdX(:,2); dz = obj.dNdX(:,3);
         end
         value =[...
            dx(1) 0     0     dy(1) 0     dz(1)
            0     dy(1) 0     dx(1) dz(1) 0
            0     0     dz(1) 0     dy(1) dx(1)
            dx(2) 0     0     dy(2) 0     dz(2)
            0     dy(2) 0     dx(2) dz(2) 0
            0     0     dz(2) 0     dy(2) dx(2)
            dx(3) 0     0     dy(3) 0     dz(3)
            0     dy(3) 0     dx(3) dz(3) 0
            0     0     dz(3) 0     dy(3) dx(3)
            dx(4) 0     0     dy(4) 0     dz(4)
            0     dy(4) 0     dx(4) dz(4) 0
            0     0     dz(4) 0     dy(4) dx(4)
            dx(5) 0     0     dy(5) 0     dz(5)
            0     dy(5) 0     dx(5) dz(5) 0
            0     0     dz(5) 0     dy(5) dx(5)
            dx(6) 0     0     dy(6) 0     dz(6)
            0     dy(6) 0     dx(6) dz(6) 0
            0     0     dz(6) 0     dy(6) dx(6)
            dx(7) 0     0     dy(7) 0     dz(7)
            0     dy(7) 0     dx(7) dz(7) 0
            0     0     dz(7) 0     dy(7) dx(7)
            dx(8) 0     0     dy(8) 0     dz(8)
            0     dy(8) 0     dx(8) dz(8) 0
            0     0     dz(8) 0     dy(8) dx(8)]';
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
            obj.U_n = permute(value,[2 1 3]);
         end
      end
      
      function obj = set.mesh(obj, val)
         [obj.det_dXdxi_list, obj.dNdX_list] =...
            obj.computeJ_and_dNdX(val.nodes, val.conn, obj.dNdxi_3D);
      end
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         ngp = size(xi,1);
         ndm = size(xi,2);
         nen = 8;
         dNdxi_3D = zeros(ndm, nen, ngp);
         for i=1:size(xi,1)
            dNdxi_3D(:,:,i) =1/8*...
               [ -(xi(i,2)- 1)*(xi(i,3)- 1), -(xi(i,1)- 1)*(xi(i,3)- 1), -(xi(i,1)- 1)*(xi(i,2)- 1)
               (xi(i,2)- 1)*(xi(i,3)- 1),  (xi(i,1)+ 1)*(xi(i,3)- 1),  (xi(i,1)+ 1)*(xi(i,2)- 1)
               -(xi(i,2)+ 1)*(xi(i,3)- 1), -(xi(i,1)+ 1)*(xi(i,3)- 1), -(xi(i,1)+ 1)*(xi(i,2)+ 1)
               (xi(i,2)+ 1)*(xi(i,3)- 1),  (xi(i,1)- 1)*(xi(i,3)- 1),  (xi(i,1)- 1)*(xi(i,2)+ 1)
               (xi(i,2)- 1)*(xi(i,3)+ 1),  (xi(i,1)- 1)*(xi(i,3)+ 1),  (xi(i,1)- 1)*(xi(i,2)- 1)
               -(xi(i,2)- 1)*(xi(i,3)+ 1), -(xi(i,1)+ 1)*(xi(i,3)+ 1), -(xi(i,1)+ 1)*(xi(i,2)- 1)
               (xi(i,2)+ 1)*(xi(i,3)+ 1),  (xi(i,1)+ 1)*(xi(i,3)+ 1),  (xi(i,1)+ 1)*(xi(i,2)+ 1)
               -(xi(i,2)+ 1)*(xi(i,3)+ 1), -(xi(i,1)- 1)*(xi(i,3)+ 1), -(xi(i,1)- 1)*(xi(i,2)+ 1)]';
         end
         dNdxi_3D = permute(dNdxi_3D, [2 1 3]);
      end
      
      function [Nmat, Ninv] = compute_Nmat(obj)
         xi = obj.xi;
         Nmat = 1/8*[...
            (1-xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3)),...
            (1+xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3)),...
            (1+xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3)),...
            (1-xi(:,1)).*(1+xi(:,2)).*(1-xi(:,3)),...
            (1-xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3)),...
            (1+xi(:,1)).*(1-xi(:,2)).*(1+xi(:,3)),...
            (1+xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3)),...
            (1-xi(:,1)).*(1+xi(:,2)).*(1+xi(:,3))];
         Ninv = inv(Nmat);
      end
      
      function [det_dXdxi_list, dNdX_list] = computeJ_and_dNdX(nodes, conn, dNdxi_3D)
         numel = size(conn    , 1);
         nen   = size(conn    , 2);
         ndm   = size(dNdxi_3D, 2);
         ngp   = size(dNdxi_3D, 3);
         
         det_dXdxi_list = zeros(numel,1);
         dNdX_list      = zeros(nen, ndm, ngp, numel);
         for i = 1:numel
            coor  = nodes(conn(i,:),:)';
            dXdxi = coor*dNdxi_3D(:,:,1);
            det_dXdxi_list(i) = det(dXdxi);
            
            for j = 1:ngp
               dNdX_list(:,:,j,i) = dNdxi_3D(:,:,j) / dXdxi';
            end
         end
      end
      
   end
end

