classdef Elements < handle
    
    properties (SetAccess = public, GetAccess = public)
        i
        K
        Fint
        list
    end    
    properties (SetAccess = private, GetAccess = public)
        Uvc
        indices
    end
    properties (SetAccess = public, Hidden)
        U_global
    end
    properties (SetAccess = private, Hidden)
        elements
        nodes
        numel
        numnp
        ndm
        ndof
        nen  
        property_num
        Props
    end
    properties (Dependent)
        Umt
        coor
        props
    end
    
    methods
        function obj = Elements(elements, nodes, numel, numnp, ndm, ndof,...
                nen, Props, U_global)
            
            obj.elements        = elements(:,1:nen);
            obj.nodes           = nodes;
            obj.numel           = numel;
            obj.numnp           = numnp;
            obj.ndm             = ndm;
            obj.ndof            = ndof;
            obj.nen             = nen;
            obj.property_num    = elements(:,end);
            obj.Props           = Props;
            obj.U_global        = U_global;
            obj.list            = cell(numel,1);
        end
        
        function value = get.coor(obj)
            value = obj.nodes(obj.elements(obj.i,:),:);
        end
        
        function value = get.Umt(obj)
            reshaped_U = reshape(obj.U_global,obj.ndof,obj.numnp);
            indc = obj.elements(obj.i,:);
            value = reshaped_U(:,indc)';
            
            obj.Uvc = reshape(value',obj.ndof*obj.nen,1);
            
            obj.indices = sub2ind(size(reshaped_U),...
                repmat((1:obj.ndof)',length(indc),1),...
                reshape(repmat(indc,obj.ndof,1),obj.ndof*obj.nen,1));
        end
        
        function clear_K_Fint(obj)
            obj.K    = zeros(obj.ndof*obj.nen, obj.ndof*obj.nen);
            obj.Fint = zeros(obj.ndof*obj.nen,1);
        end
        
        function value = get.props(obj)
            value = obj.Props( obj.property_num(obj.i),:);
        end
        
    end
end

