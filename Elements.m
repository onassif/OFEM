classdef Elements
    
    properties (SetAccess = public, GetAccess = public)
        i
        K
        Fint
    end    
    properties (SetAccess = private, GetAccess = public)
        Umt
        Uvc
        coor
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
        props
    end
    
    methods
        %% Construct
        function obj = Elements(elements, nodes, num, Props, U_global)
            obj.elements     = elements(:,1:num.nen);
            obj.nodes        = nodes;
            obj.numel        = num.el;
            obj.numnp        = num.np;
            obj.ndm          = num.ndm;
            obj.ndof         = num.ndof;
            obj.nen          = num.nen;
            obj.property_num = elements(:,end);
            obj.props        = Props;
            obj.U_global     = U_global;
            if num.ndof == num.ndm
                obj.K            = zeros(num.ndof*num.nen, num.ndof*num.nen);
                obj.Fint         = zeros(num.ndof*num.nen, 1);
            elseif (num.ndof - num.ndm) == 1
                sz = (num.ndof-1)*num.nen + 1;
                obj.K            = zeros(sz, sz);
                obj.Fint         = zeros(sz,  1);
            else
                error("ndof-ndm > 1")
            end
        end
        %% get functions
        function value = get.coor(obj)
            value = obj.nodes(obj.elements(obj.i,:),:);
        end
        function value = get.Umt(obj)
            reshaped_U = reshape(obj.U_global(1:obj.ndm*obj.numnp), obj.ndm, obj.numnp);
            value = reshaped_U(:, obj.elements(obj.i,:))';
        end      
        function value = get.Uvc(obj)
            value = reshape(obj.Umt', obj.ndm*obj.nen, 1);
        end
        function value = get.indices(obj)
            indc = obj.elements(obj.i,:);
            
%             dofs_col = repmat((1:obj.ndof)', obj.nen,1);
            dofs_col = repmat((1:obj.ndm)', obj.nen,1);
            
%             repeated_el_conn = repmat(indc, obj.ndof, 1);
            repeated_el_conn = repmat(indc, obj.ndm, 1);
%             conns_col= reshape(repeated_el_conn, obj.ndof*obj.nen, 1);
            conns_col= reshape(repeated_el_conn, obj.ndm*obj.nen, 1);
            
%             value = sub2ind([obj.ndof obj.numnp], dofs_col, conns_col);
            value = sub2ind([obj.ndm  obj.numnp], dofs_col, conns_col);
            if (obj.ndof - obj.ndm) == 1
                value = [value; obj.numnp*obj.ndm+obj.i];
            end
        end
        
        %% set functions
        function obj = set.K(obj, value)
            if isempty(obj.K) || sum(size(value)==size(obj.K))==2 % 1st time or same size
                obj.K = value;
            elseif value == 0
                obj.K = zeros(size(obj.K));
            else
                error("Size of element K has changed")
            end
        end
        function obj = set.Fint(obj, value)
            if isempty(obj.Fint) || sum(size(value)==size(obj.Fint))==2
                obj.Fint = value;
            elseif value == 0
                obj.Fint = zeros(size(obj.Fint));
            else
                error("Size of element Fint has changed") 
            end
        end
    end
end

