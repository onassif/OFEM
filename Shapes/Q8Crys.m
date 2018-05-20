classdef Q8Crys
   properties (SetAccess = public, GetAccess = public)
      i
      dXdxi;
      dxdX;
      dNdx;
      dNdx_hf;
      b;
      eps;
      sigma;
      ctan;
      D;
      U;
      U_n;
      angles = [0, 0, 0];
      Rpn_list
      Rn_list;
      dXdxi_list;
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
      I;
      bi;
      ni;
      nslip;
   end
   
   properties (SetAccess = private)
      Nmat;
      Ninv;
      dNdX;
      B;
      B_hf;
      N;
      w;
      dNdxi;
      F;
      R;
      q;
      F_hf;
      R_hf;
      q_hf;
      sigma_un;
      D_un;
      J; % det(dX/dxi) = J
      j; % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      gRot;
      Rpn;
      
      ms;
      qs;
      q_cr;
      xi = 1/sqrt(3) .*[...
         -1  1  1 -1 -1  1  1 -1
         -1 -1  1  1 -1 -1  1  1
         -1 -1 -1 -1  1  1  1  1]';
     weights = [1 1 1 1 1 1 1 1]; 
   end
   
   methods
      %% Construct
      function obj = Q8Crys(varargin)
         num        = varargin{1};
         finiteDisp = varargin{2};
         if nargin == 3
            obj.xi = varargin{3};
         end
         obj.dNdxi_3D         = obj.compute_dNdxi(obj);
         [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
         
         obj.det_dXdxi_list = zeros(8,num.el);
         obj.dXdxi_list     = zeros(3,3,num.el);
         obj.dNdX_list      = zeros(8,3,8,num.el);
         
         obj.Rpn_list = zeros(3,3,8,num.el);
         obj.Rpn_list(1,1,:,:) = 1;
         obj.Rpn_list(2,2,:,:) = 1;
         obj.Rpn_list(3,3,:,:) = 1;
         
         obj.Rn_list = obj.Rpn_list;
         
         obj.bi = slip.b;
         obj.ni = slip.n;
         obj.nslip = slip.nslip;
         
         obj.I = eye(num.ndm);
         
         if (strcmp(angles.conv, 'kocks') && strcmp(angles.type, 'degrees'))
            obj.angles = angles.val;
         else
            error("Wrong angle convention/type");
         end
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
      
      function value = get.Rpn(obj)
         value = obj.Rpn_list(:,:,obj.i,obj.iel);
      end
      
      function value = get.F(obj)
        value = obj.U*obj.dNdX + obj.I;
      end
      
      function value = get.R(obj)
           [P, ~, Q] = svd(obj.F);
           value =  P*Q';
      end
      
      function value = get.Rn_list(obj)
           value(:,:,obj.i,obj.iel) = obj.R;
      end
      
      function value = get.dXdxi_list(obj)
           value(:,:,obj.iel) = obj.dXdxi;
      end
      
      function value = get.q(obj)
         R11 = obj.R(1,1); R12 = obj.R(1,2); R13 = obj.R(1,3);
         R21 = obj.R(2,1); R22 = obj.R(2,2); R23 = obj.R(2,3);
         R31 = obj.R(3,1); R32 = obj.R(3,2); R33 = obj.R(3,3);
         value = [...
            R11^2   	R12^2   	R13^2   	2*R11*R12        	2*R12*R13        	2*R11*R13
            R21^2   	R22^2   	R23^2   	2*R21*R22        	2*R22*R23        	2*R21*R23 
            R31^2    R32^2   	R33^2   	2*R31*R32        	2*R32*R33        	2*R31*R33
            R11*R21  R12*R22  R13*R23  R11*R22+R12*R21	R12*R23+R13*R22   R11*R23+R13*R21
            R21*R31  R32*R22  R23*R33  R21*R32+R22*R31	R22*R33+R23*R32   R21*R33+R23*R31
            R11*R31  R12*R32  R13*R33  R11*R32+R12*R31	R12*R33+R13*R32   R11*R33+R13*R31            
            ];
      end

      function value = get.F_hf(obj)
        value = (1/2)*(obj.U + obj.U_n)*obj.dNdX + obj.I;
      end
      
      function value = get.R_hf(obj)
           [P, ~, Q] = svd(obj.F_hf);
           value =  P*Q';
      end
      
      function value = get.q_hf(obj)
         R11 = obj.R_hf(1,1); R12 = obj.R_hf(1,2); R13 = obj.R_hf(1,3);
         R21 = obj.R_hf(2,1); R22 = obj.R_hf(2,2); R23 = obj.R_hf(2,3);
         R31 = obj.R_hf(3,1); R32 = obj.R_hf(3,2); R33 = obj.R_hf(3,3);
         value = [...
            R11^2   	R12^2   	R13^2   	2*R11*R12        	2*R12*R13        	2*R11*R13
            R21^2   	R22^2   	R23^2   	2*R21*R22        	2*R22*R23        	2*R21*R23 
            R31^2    R32^2   	R33^2   	2*R31*R32        	2*R32*R33        	2*R31*R33
            R11*R21  R12*R22  R13*R23  R11*R22+R12*R21	R12*R23+R13*R22   R11*R23+R13*R21
            R21*R31  R32*R22  R23*R33  R21*R32+R22*R31	R22*R33+R23*R32   R21*R33+R23*R31
            R11*R31  R12*R32  R13*R33  R11*R32+R12*R31	R12*R33+R13*R32   R11*R33+R13*R31            
            ];
      end
      
      function value = get.sigma_un(obj)
         value = [...
            obj.q*obj.sigma
            0
            0
            0];
      end
      
      function value = get.j(obj)     
         value = det(obj.F) * obj.J;
      end
      
      function value = get.dNdx(obj)
         value = obj.dNdX / obj.F;
      end

      function value = get.dNdx_hf(obj)
         value = obj.dNdX / obj.F_hf;
      end
      
      function value = get.b(obj)
         if obj.finiteDisp
            value = obj.F*obj.F';
         end
      end
      
      function value = get.gRot(obj)
            r = obj.angles(1);
            s = obj.angles(2);
            t = obj.angles(3);
         
            value = [...
               -sin(r)*sin(t)-cos(r)*cos(t)*cos(s), cos(r)*sin(t)-sin(r)*cos(t)*cos(s), cos(t)*sin(s)
                sin(r)*cos(t)-cos(r)*sin(t)*cos(s),-cos(r)*cos(t)-sin(r)*sin(t)*cos(s), sin(t)*sin(s)
                cos(r)*sin(s)                     , sin(r)*sin(s)                     , cos(s)       ];   
      end
      
      function value = get.ms(obj)       
         value = zeros(3,3,obj.nslip);
         bs = obj.bi * obj.gRot;
         ns = obj.ni * obj.gRot;
         
         for s = 1:obj.nslip
            A = bs(s,:)'*ns(s,:);
            value(:,:,s) = obj.Rpn'*(1/2)*(A+A')*obj.Rpn;
         end
      end

      function value = get.qs(obj)
         value = zeros(3,3,obj.nslip);
         bs = obj.bi * obj.gRot;
         ns = obj.ni * obj.gRot;
         
         for s = 1:obj.nslip
            A = bs(s,:)'*ns(s,:);
            value(:,:,s) = obj.Rpn'*(1/2)*(A-A')*obj.Rpn;
         end
      end

      function value = get.q_cr(obj)
         value = obj.R * obj.qs * obj.R';
      end
      
      function value = get.B(obj)
            dx = obj.dNdx(:,1); dy = obj.dNdx(:,2); dz = obj.dNdx(:,3);
            
         value =[...
            dx(1)  0      0      dx(2)  0      0      dx(3)  0      0      dx(4)  0      0      dx(5)  0      0      dx(6)  0      0      dx(7)  0      0      dx(8)  0      0     
            0      dy(1)  0      0      dy(2)  0      0      dy(3)  0      0      dy(4)  0      0      dy(5)  0      0      dy(6)  0      0      dy(7)  0      0      dy(8)  0    
            0      0      dz(1)  0      0      dz(2)  0      0      dz(3)  0      0      dz(4)  0      0      dz(5)  0      0      dz(6)  0      0      dz(7)  0      0      dz(8)
            dy(1)  dx(1)  0      dy(2)  dx(2)  0      dy(3)  dx(3)  0      dy(4)  dx(4)  0      dy(5)  dx(5)  0      dy(6)  dx(6)  0      dy(7)  dx(7)  0      dy(8)  dx(8)  0     
            0      dz(1)  dy(1)  0      dz(2)  dy(2)  0      dz(3)  dy(3)  0      dz(4)  dy(4)  0      dz(5)  dy(5)  0      dz(6)  dy(6)  0      dz(7)  dy(7)  0      dz(8)  dy(8)  
            dz(1)  0      dx(1)  dz(2)  0      dx(2)  dz(3)  0      dx(3)  dz(4)  0      dx(4)  dz(5)  0      dx(5)  dz(6)  0      dx(6)  dz(7)  0      dx(7)  dz(8)  0      dx(8) 
            dy(1) -dx(1)  0      dy(2) -dx(2)  0      dy(3) -dx(3)  0      dy(4) -dx(4)  0      dy(5) -dx(5)  0      dy(6) -dx(6)  0      dy(7) -dx(7)  0      dy(8) -dx(8)  0     
            0      dz(1) -dx(1)  0      dz(2) -dx(2)  0      dz(3) -dx(3)  0      dz(4) -dx(4)  0      dz(5) -dx(5)  0      dz(6) -dx(6)  0      dz(7) -dx(7)  0      dz(8) -dx(8)  
           -dz(1)  0      dx(1) -dz(2)  0      dx(2) -dz(3)  0      dx(3) -dz(4)  0      dx(4) -dz(5)  0      dx(5) -dz(6)  0      dx(6) -dz(7)  0      dx(7) -dz(8)  0      dx(8) 
           ];
      end
      
      function value = get.B_hf(obj)
            dx = obj.dNdx_hf(:,1); dy = obj.dNdx_hf(:,2); dz = obj.dNdx_hf(:,3);
            
         value =[...
            dx(1)  0      0      dx(2)  0      0      dx(3)  0      0      dx(4)  0      0      dx(5)  0      0      dx(6)  0      0      dx(7)  0      0      dx(8)  0      0     
            0      dy(1)  0      0      dy(2)  0      0      dy(3)  0      0      dy(4)  0      0      dy(5)  0      0      dy(6)  0      0      dy(7)  0      0      dy(8)  0    
            0      0      dz(1)  0      0      dz(2)  0      0      dz(3)  0      0      dz(4)  0      0      dz(5)  0      0      dz(6)  0      0      dz(7)  0      0      dz(8)
            dy(1)  dx(1)  0      dy(2)  dx(2)  0      dy(3)  dx(3)  0      dy(4)  dx(4)  0      dy(5)  dx(5)  0      dy(6)  dx(6)  0      dy(7)  dx(7)  0      dy(8)  dx(8)  0     
            0      dz(1)  dy(1)  0      dz(2)  dy(2)  0      dz(3)  dy(3)  0      dz(4)  dy(4)  0      dz(5)  dy(5)  0      dz(6)  dy(6)  0      dz(7)  dy(7)  0      dz(8)  dy(8)  
            dz(1)  0      dx(1)  dz(2)  0      dx(2)  dz(3)  0      dx(3)  dz(4)  0      dx(4)  dz(5)  0      dx(5)  dz(6)  0      dx(6)  dz(7)  0      dx(7)  dz(8)  0      dx(8) 
            dy(1) -dx(1)  0      dy(2) -dx(2)  0      dy(3) -dx(3)  0      dy(4) -dx(4)  0      dy(5) -dx(5)  0      dy(6) -dx(6)  0      dy(7) -dx(7)  0      dy(8) -dx(8)  0     
            0      dz(1) -dx(1)  0      dz(2) -dx(2)  0      dz(3) -dx(3)  0      dz(4) -dx(4)  0      dz(5) -dx(5)  0      dz(6) -dx(6)  0      dz(7) -dx(7)  0      dz(8) -dx(8)  
           -dz(1)  0      dx(1) -dz(2)  0      dx(2) -dz(3)  0      dx(3) -dz(4)  0      dx(4) -dz(5)  0      dx(5) -dz(6)  0      dx(6) -dz(7)  0      dx(7) -dz(8)  0      dx(8) 
           ];
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

