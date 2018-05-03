classdef Q8Crys
   
   properties (Constant)
      % Isoparametric gp
      xi = 1/sqrt(3) .*[...
         -1  1  1 -1 -1  1  1 -1
         -1 -1  1  1 -1 -1  1  1
         -1 -1 -1 -1  1  1  1  1]';
   end
   
   properties (SetAccess = public, GetAccess = public)
      i
      dXdxi;
      dxdX;
      dNdx;
      F;
      det_F;
      b;
      J;
      j;
      eps;
      sigma;
      ctan;
      D;
   end
   
   properties (Hidden)
      iel;
      det_dXdxi_list;
      dNdX_list;
   end
   
   properties (Hidden, SetAccess = private)
      dNdxi_3D;
      numel;
      weights = [1 1 1 1 1 1 1 1];
      finiteDisp;
      adof;
   end
   
   properties (SetAccess = private)
      Nmat;
      Ninv;
      det_dXdxi;
      dNdX;
      B;
      N;
      w;
      dNdxi;
      
      gRot
      R
      np1_q
      mat;
   end
   
   methods
      function obj = Q8Crys(num, mat)
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
         
         obj.det_dXdxi_list = zeros(8,num.el);
         obj.dNdX_list      = zeros(8,3,8,num.el);
         
         obj.mat = mat;
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
      
      function value = get.det_dXdxi(obj)
         value = obj.det_dXdxi_list(obj.i,obj.iel);
      end
      
      function value = get.dNdX(obj)
         value = obj.dNdX_list(:,:,obj.i,obj.iel);
      end
      
      function value = get.dNdx(obj)
         value = obj.dNdX / obj.F;
      end
      
      function value = get.b(obj)
         if obj.finiteDisp
            value = obj.F*obj.F';
         end
      end
 
      function value = get.R(obj)
         [P,~,Q] = svd(obj.F);
         value   = P*Q';
      end
      
      function value = get.np1_q(obj)
         r = obj.R;
         value= [...
            r(1,1)^2, r(1,2)^2, r(1,3)^2, ...
            2*r(1,1)*r(1,2), ...
            2*r(1,3)*r(1,2), ...
            2*r(1,1)*r(1,3);
            ...
            r(2,1)^2, r(2,2)^2, r(2,3)^2, ...
            two*r(2,1)*r(2,2), ...
            two*r(2,3)*r(2,2), ...
            two*r(2,1)*r(2,3);
            ...
            r(3,1)^2, r(3,2)^2, r(3,3)^2, ...
            two*r(3,1)*r(3,2), ...
            two*r(3,3)*r(3,2), ...
            two*r(3,1)*r(3,3);
            ...
            r(1,1)*r(2,1), r(1,2)*r(2,2), r(1,3)*r(2,3), ...
            r(1,1)*r(2,2)+r(2,1)*r(1,2), ...
            r(1,2)*r(2,3)+r(1,3)*r(2,2), ...
            r(1,1)*r(2,3)+r(1,3)*r(2,1);
            ...
            r(2,1)*r(3,1), r(3,2)*r(2,2), r(2,3)*r(3,3), ...
            r(2,1)*r(3,2)+r(2,2)*r(3,1), ...
            r(2,2)*r(3,3)+r(3,2)*r(2,3), ...
            r(2,1)*r(3,3)+r(2,3)*r(3,1);
            ...
            r(1,1)*r(3,1), r(1,2)*r(3,2), r(1,3)*r(3,3), ...
            r(1,1)*r(3,2)+r(1,2)*r(3,1), ...
            r(1,2)*r(3,3)+r(1,3)*r(3,2), ...
            r(1,1)*r(3,3)+r(3,1)*r(1,3)];
      end
            
      function value = get.B(obj)
         if (obj.finiteDisp)
            dx = obj.dNdx(:,1); dy = obj.dNdx(:,2); dz = obj.dNdx(:,3);
         else
            dx = obj.dNdX(:,1); dy = obj.dNdX(:,2); dz = obj.dNdX(:,3);
         end
         a = zeros(obj.adof,1);
         value =[...
            dx(1) 0     0     dy(1) 0     dz(1)
            0     dy(1) 0     dx(1) dz(1) 0
            0     0     dz(1) 0     dy(1) dx(1)
            a     a     a     a     a     a
            dx(2) 0     0     dy(2) 0     dz(2)
            0     dy(2) 0     dx(2) dz(2) 0
            0     0     dz(2) 0     dy(2) dx(2)
            a     a     a     a     a     a
            dx(3) 0     0     dy(3) 0     dz(3)
            0     dy(3) 0     dx(3) dz(3) 0
            0     0     dz(3) 0     dy(3) dx(3)
            a     a     a     a     a     a
            dx(4) 0     0     dy(4) 0     dz(4)
            0     dy(4) 0     dx(4) dz(4) 0
            0     0     dz(4) 0     dy(4) dx(4)
            a     a     a     a     a     a
            dx(5) 0     0     dy(5) 0     dz(5)
            0     dy(5) 0     dx(5) dz(5) 0
            0     0     dz(5) 0     dy(5) dx(5)
            a     a     a     a     a     a
            dx(6) 0     0     dy(6) 0     dz(6)
            0     dy(6) 0     dx(6) dz(6) 0
            0     0     dz(6) 0     dy(6) dx(6)
            a     a     a     a     a     a
            dx(7) 0     0     dy(7) 0     dz(7)
            0     dy(7) 0     dx(7) dz(7) 0
            0     0     dz(7) 0     dy(7) dx(7)
            a     a     a     a     a     a
            dx(8) 0     0     dy(8) 0     dz(8)
            0     dy(8) 0     dx(8) dz(8) 0
            0     0     dz(8) 0     dy(8) dx(8)
            a     a     a     a     a     a    ]';
      end
      
      function value = get.gRot(obj)
         if strcmp(obj.mat.hard.angleConvention, 'kocks')
            
            if (strcmp(obj.mat.hard.angleType, 'degrees'))
               psi   = obj.mat.hard.angles(1) * (pi/180);
               theta = obj.mat.hard.angles(2) * (pi/180);
               phi   = obj.mat.hard.angles(3) * (pi/180);
            elseif (strcmp(atype, 'radians'))
               psi   = obj.mat.hard.angles(1);
               theta = obj.mat.hard.angles(2);
               phi   = obj.mat.hard.angles(3);
            end
            
            value = [...
               -sin(psi)*sin(phi)-cos(psi)*cos(phi)*cos(theta),...
               cos(psi)*sin(phi)-sin(psi)*cos(phi)*cos(theta) ,...
               cos(phi)*sin(theta);
               ...
               sin(psi)*cos(phi)-cos(psi)*sin(phi)*cos(theta) ,...
               -cos(psi)*cos(phi)-sin(psi)*sin(phi)*cos(theta),...
               sin(phi)*sin(theta);
               ...
               cos(psi)*sin(theta),...
               sin(psi)*sin(theta),...
               cos(theta)];
         else
            error("Wrong angle convention");
         end
      end
   end
   
   methods (Static)
      function dNdxi_3D = compute_dNdxi(obj)
         xi = obj.xi;
         dNdxi_3D = zeros(3,8,8);
         for i=1:8
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
   end
end

