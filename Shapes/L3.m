classdef L3
   properties
      i
      iel
      
      b
      eps
      sigma
      D
      U
      U_n
      dU
      mesh
      xi = [-sqrt(0.6) 0 sqrt(0.6)];
      weights = (1/9).*[5 8 5]';
      Nmat
      dNdxi_list
      d2Ndxi2_list
      
      dNdX_list
      dXdxi_list
      det_dXdxi_list
   end
   
   properties (SetAccess = private)
      finiteDisp
      Ninv
      dNdX
      dNdx
      dNdxi
      d2Ndxi2
      dXdxi
      B
      N
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
   end
   
   methods
      %% Construct
      function ob = L3(varargin)
         ob.finiteDisp = varargin{1};
      end
      %% Get functions
      function value = get.N(ob)
         value = ob.Nmat(:,ob.i);
      end
      
      function value = get.dNdxi(ob)
         value = ob.dNdxi_list(:,:,ob.i);
      end
      
      function value = get.d2Ndxi2(ob)
         value = ob.d2Ndxi2_list(:,:,ob.i);
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(ob.i,ob.iel);
      end
      
      function value = get.w(ob)
         value = ob.weights(ob.i);
      end
      
      function value = get.J(ob)
         value = ob.det_dXdxi_list(ob.i, ob.iel);
      end
      
      function value = get.dNdX(ob)
         value = ob.dNdX_list(:,:,ob.i,ob.iel);
      end
      
      function value = get.F(ob)
         I     = eye(size(ob.U,1));
         value = ob.U*ob.dNdX + I;
      end
      
      function value = get.JxX(ob)
         value = det(ob.F);
      end
      
      function value = get.j(ob)
         value = det(ob.F) * ob.J;
      end
      
      function value = get.dNdx(ob)
         value = ob.dNdX / ob.F;
      end
      
      function value = get.b(ob)
         value = ob.F*ob.F';
      end
      
      function value = get.B(ob)
         if (ob.finiteDisp)
            dx = ob.dNdx;
         else
            dx = ob.dNdX;
         end
         value=[dx(1), dx(2), dx(3)];
      end
      %% xi-dependant functions
      function ob = shapeIso(ob)
         ndm = size(ob.xi,1); 
         ngp = size(ob.xi,2);
         nen = 3;
         x1  = ob.xi(1,:);
         
         ob.Nmat = 1/2*[...
            x1.*(x1-1)
            -2*(x1+1).*(x1-1)
            x1.*(x1+1)];
         if size(ob.Nmat,1) == size(ob.Nmat,2)
            ob.Ninv = inv(ob.Nmat);
         else
            ob.Ninv = 0;
         end
         ob.dNdxi_list = zeros(nen,ndm,ngp);
         ob.dNdxi_list(:,1,:) = 1/2 *[2*x1-1; -4*x1; 2*x1+1];

         vec = ones(1,ngp);
         ob.d2Ndxi2_list = [...
            1*vec
            -2*vec
            1*vec];
         ob.d2Ndxi2_list = permute(ob.d2Ndxi2_list,[1,3,2]);
      end
   end
end