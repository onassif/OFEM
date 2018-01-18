classdef Q8
    
    properties (Constant)
        % Isoparametric gp
        xi = 1/sqrt(3) .*...
            [-1  1  1 -1 -1  1  1 -1
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
    end
    
    methods
        function obj = Q8(num, finiteDisp)
            obj.dNdxi_3D = obj.compute_dNdxi(obj);
            [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
            obj.det_dXdxi_list = zeros(8,num.el);
            obj.dNdX_list = zeros(8,3,8,num.el);
            obj.finiteDisp = finiteDisp;
            obj.adof = num.ndof - num.ndm;
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
            value = obj.F*obj.F';
        end
        
        function value = get.B(obj)
            if (obj.finiteDisp)
                dx = obj.dNdx(:,1); dy = obj.dNdx(:,2); dz = obj.dNdx(:,3);
            else
                dx = obj.dNdX(:,1); dy = obj.dNdX(:,2); dz = obj.dNdX(:,3);
            end
            a = zeros(obj.adof,1);
            value =[dx(1) 0     0     dy(1) 0     dz(1)
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
            Nmat = 1/8*[(1-xi(:,1)).*(1-xi(:,2)).*(1-xi(:,3)),...
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

