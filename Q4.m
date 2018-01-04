classdef Q4
    
    properties (Constant)
        % Isoparametric gp
        xi = 1/sqrt(3) .*...
            [-1  1  1 -1
            -1 -1  1  1]';
    end
    
    properties (SetAccess = public, GetAccess = public)
        i
        dXdxi;
        dxdX;
        dNdx;
        F;
        b;
        J;
        j;
        j0;
        eps;
        sigma;
        ctan;
        D;
        B;
        det_F;
    end
    
    properties (Dependent)
        N;
        dNdxi;
        w;
        det_dXdxi;
        dNdX;
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
    end
    
    properties (SetAccess = private)
        Nmat;
        Ninv;
    end
    
    methods
        function obj = Q4(numel)
            obj.dNdxi_3D = obj.compute_dNdxi(obj);
            [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
            obj.det_dXdxi_list = zeros(4,numel);
            obj.dNdX_list = zeros(4,2,4,numel);
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
    end
    
    methods (Static)
        function dNdxi_3D = compute_dNdxi(obj)
            xi = obj.xi;
            dNdxi_3D = zeros(2,4,4);
            for i=1:4
                dNdxi_3D(:,:,i) =...
                    1/4 *[...
                    -(1-xi(i,2)) -(1-xi(i,1))
                    (1-xi(i,2)) -(1+xi(i,1))
                    (1+xi(i,2))  (1+xi(i,1))
                    -(1+xi(i,2))  (1-xi(i,1))]';
            end
        end
        
        function [Nmat, Ninv] = compute_Nmat(obj)
            xi = obj.xi;
            Nmat = 1/4*[(1-xi(:,1)).*(1-xi(:,2)), (1+xi(:,1)).*(1-xi(:,2)),...
                (1+xi(:,1)).*(1+xi(:,2)), (1-xi(:,1)).*(1+xi(:,2))];
            Ninv = inv(Nmat);
        end
    end
end

