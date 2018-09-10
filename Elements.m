classdef Elements
   properties
      i
      K
      Fint
      iter = 0;
      U_global
      w_global
      U_glb_n
      ulres
   end
   properties (SetAccess = private, GetAccess = public)
      w
      Umt
      Umt_n
      Uvc
      Uvc_n
      coor
      indices
      conn
      nconn
      im
      elements
      nodes
   end
   properties (SetAccess = private, Hidden)
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
      function obj = Elements(elements, nodes, num, Props, U_global, hist)
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
         obj.conn         = hist.conn;
         obj.nconn        = hist.nconn;
         if num.ndof == num.ndm
            obj.K            = zeros(num.ndof*num.nen, num.ndof*num.nen);
            obj.Fint         = zeros(num.ndof*num.nen, 1);
         elseif (num.ndof - num.ndm) == 1
            switch num.nen
               case {3,4}
                  sz = (num.ndof-1)*num.nen + 1;
               case 6
                  sz = (num.ndof-1)*num.nen + 3;
               case 9
                  sz = (num.ndof-1)*num.nen + 4;
            end
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
         if (obj.i <= obj.numel) || ~rem(obj.i - obj.numel,2)
            value = reshaped_U(:, obj.elements(obj.i,:))';
         else % DG
               value(:,:,1) = reshaped_U(:, obj.elements(obj.i  ,:))';
               value(:,:,2) = reshaped_U(:, obj.elements(obj.i+1,:))';
         end
      end
      function value = get.w(obj)
         reshaped_w = reshape(obj.w_global(1:obj.ndm*obj.numnp), obj.ndm, obj.numnp);
         if (obj.i <= obj.numel) || ~rem(obj.i - obj.numel,2)
            value = reshaped_w(:, obj.elements(obj.i,:))';
         else % DG
            value(:,:,1) = reshaped_w(:, obj.elements(obj.i  ,:))';
            value(:,:,2) = reshaped_w(:, obj.elements(obj.i+1,:))';
         end
      end
      function value = get.Umt_n(obj)
         reshaped_U = reshape(obj.U_glb_n(1:obj.ndm*obj.numnp), obj.ndm, obj.numnp);
         value = reshaped_U(:, obj.elements(obj.i,:))';
      end
      function value = get.Uvc(obj)
         if (obj.i <= obj.numel) || ~rem(obj.i - obj.numel,2)
            value = reshape(obj.Umt', obj.ndm*obj.nen, 1);
         else % DG
            value(:,1) = reshape(obj.Umt(:,:,1)', obj.ndm*obj.nen, 1);
            value(:,2) = reshape(obj.Umt(:,:,2)', obj.ndm*obj.nen, 1);
         end
      end
      function value = get.Uvc_n(obj)
         if (obj.i <= obj.numel) || ~rem(obj.i - obj.numel,2)
            value = reshape(obj.Umt_n', obj.ndm*obj.nen, 1);
         else % DG
            value(:,1) = reshape(obj.Umt_n(:,:,1)', obj.ndm*obj.nen, 1);
            value(:,2) = reshape(obj.Umt_n(:,:,2)', obj.ndm*obj.nen, 1);
         end
      end
      function value = get.im(obj)
         value = obj.property_num(obj.i);
      end
      function value = get.indices(obj)
         if (obj.i <= obj.numel) || ~rem(obj.i - obj.numel,2)
            indc = obj.elements(obj.i,:);
            dofs_col = repmat((1:obj.ndm)', obj.nen,1);
         else
            indc = [obj.elements(obj.i,:), obj.elements(obj.i+1,:)];
            dofs_col = repmat((1:obj.ndm)', 2*obj.nen,1);
         end
         
         repeated_el_conn = repmat(indc, obj.ndm, 1);

         conns_col= reshape(repeated_el_conn, numel(repeated_el_conn), 1);
         
         value = sub2ind([obj.ndm  obj.numnp], dofs_col, conns_col);
         if (obj.ndof - obj.ndm) == 1
            switch obj.nen
               case {3,4}
                  extra = obj.i;
               case 6
                  extra = (1 + 3*(obj.i-1):3*obj.i)';
               case 9
                  extra = (1 + 4*(obj.i-1):4*obj.i)';
            end
            value = [value; obj.numnp*obj.ndm+extra];
         end
      end
      
      %% set functions
      function obj = set.U_global(obj, value)
         if obj.iter == 0
            if ~isempty(obj.U_global)
               obj.U_glb_n  = obj.U_global;
            end
         end
         obj.U_global = value;
      end
   end
end

