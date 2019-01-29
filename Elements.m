classdef Elements
   properties
      i
      K
      Fint
      M
      iter = 0;
      U_global
      w_global
      U_glb_n
      Ures_glb
      mat
   end
   properties (SetAccess = private)
      w
      Ures
      UresVc
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
      nen
      property_num
   end
   properties (SetAccess = private, Hidden)
      numnp
      ndm
      ndof
   end
   
   methods
      %% Construct
      function ob = Elements(elements, nodes, num, U_global, hist)
         ob.elements     = elements(:,1:end-1);
         ob.nodes        = nodes;
         ob.numnp        = num.np;
         ob.ndm          = num.ndm;
         ob.ndof         = num.ndof;
         ob.nen          = num.nen;
         ob.property_num = elements(:,end);
         ob.U_global     = U_global;
         ob.conn         = hist.conn;
         ob.nconn        = hist.nconn;

         ob.Ures_glb = zeros(num.ndm*num.np,1);
         ob.w_global = zeros(num.ndm*num.np,1);
         ob.U_global = zeros(num.ndm*num.np,1);
         ob.U_glb_n  = zeros(num.ndm*num.np,1);
      end
      %% get functions
      function value = get.coor(ob)
         value = ob.nodes(nonzeros(ob.elements(ob.i,:)),:);
      end
      function value = get.Umt(ob)
         reshaped_U = reshape(ob.U_global(1:ob.ndm*ob.numnp), ob.ndm, ob.numnp);
         value = reshaped_U(:, nonzeros(ob.elements(ob.i,:)));
      end
      function value = get.w(ob)
         reshaped_w = reshape(ob.w_global(1:ob.ndm*ob.numnp), ob.ndm, ob.numnp);
         value = reshaped_w(:, nonzeros(ob.elements(ob.i,:)));
      end
      function value = get.Ures(ob)
         reshaped = reshape(ob.Ures_glb(1:ob.ndm*ob.numnp), ob.ndm, ob.numnp);
         value = reshaped(:, nonzeros(ob.elements(ob.i,:)));
      end
      function value = get.Umt_n(ob)
         reshaped_U = reshape(ob.U_glb_n(1:ob.ndm*ob.numnp), ob.ndm, ob.numnp);
         value = reshaped_U(:, nonzeros(ob.elements(ob.i,:)));
      end
      function value = get.Uvc(ob)
         value = reshape(ob.Umt, numel(ob.Umt), 1);
      end
      function value = get.Uvc_n(ob)
         value = reshape(ob.Umt_n, numel(ob.Umt_n), 1);
      end
      function value = get.UresVc(ob)
         value = reshape(ob.Ures, numel(ob.Ures), 1);
      end
      function value = get.im(ob)
         value = ob.property_num(ob.i);
      end
      function value = get.indices(ob)
         indc = nonzeros(ob.elements(ob.i,:))';
         dofs_col = repmat((1:ob.ndm)', length(indc),1);
         
         repeated_el_conn = repmat(indc, ob.ndm, 1);

         conns_col= reshape(repeated_el_conn, numel(repeated_el_conn), 1);
         
         value = sub2ind([ob.ndm  ob.numnp], dofs_col, conns_col);
         if (ob.ndof - ob.ndm) == 1
            switch ob.nen
               case {3,4}
                  extra = ob.i;
               case 6
                  extra = (1 + 3*(ob.i-1):3*ob.i)';
               case 9
                  extra = (1 + 4*(ob.i-1):4*ob.i)';
            end
            value = [value; ob.numnp*ob.ndm+extra];
         end
      end
      
   end
end

