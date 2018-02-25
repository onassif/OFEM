classdef T3
    
    properties (Constant)
        % Isoparametric gp
        xi = 1/3 .*[1 1];
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
        deps;
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
        finiteDisp;
    end
    
    properties (SetAccess = private)
        Nmat;
        Ninv;
        det_dXdxi;
        dNdX;
        B;
        N;
        w=0.5;
        dNdxi;
    end
    
    methods
        function obj = T3(numel, finiteDisp)
            obj.dNdxi_3D = obj.compute_dNdxi(obj);
            [obj.Nmat, obj.Ninv] = obj.compute_Nmat(obj);
            obj.det_dXdxi_list = zeros(1,numel);
            obj.dNdX_list = zeros(3,2,1,numel);
            obj.finiteDisp = finiteDisp;
        end
        
        function value = get.dNdxi(obj)
            value = squeeze(obj.dNdxi_3D(:,:,obj.i));
        end
        
        function value = get.N(obj)
            value = obj.Nmat;
        end
        
        
        function value = get.det_dXdxi(obj)
            value = obj.det_dXdxi_list(obj.iel);
        end
        
        function value = get.dNdX(obj)
            value = obj.dNdX_list(:,:,obj.iel);
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
                dx(1),   0.0, dx(2),   0.0, dx(3),   0.0
                0.0,   dy(1),   0.0, dy(2),   0.0, dy(3)
                dy(1), dx(1), dy(2), dx(2), dy(3), dx(3)];
        end
        
    end
    
    methods (Static)
        function dNdxi_3D = compute_dNdxi(~)
            dNdxi_3D =[...
                -1 1 0
                -1 0 1];
        end
        
        function [Nmat, Ninv] = compute_Nmat(obj)
            xi = obj.xi;
            Nmat = [1-xi(1)-xi(2), xi(1), xi(2)];
            Ninv = Nmat';
        end
    end
end

