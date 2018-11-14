classdef Q9
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
      xi = sqrt(0.6)*[...
         -1 +1 +1 -1  0 +1  0 -1  0
         -1 -1 +1 +1 -1  0 +1  0  0];
      weights = (1/81).*[25 25 25 25 40 40 40 40 64]';
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
      Bf
      N
      R
      w
      F
      J % det(dX/dxi) = J
      j % or det( dx/dX*dX/dxi ) = det(dx/dxi) = j
      JxX
   end
   
   methods
      %% Construct
      function ob = Q9(varargin)
         ob.finiteDisp = varargin{1};
      end
      %% Get functions
      function value = get.N(ob)
         value = ob.Nmat(:,ob.i);
      end
      
      function value = get.dNdxi(ob)
         value = squeeze(ob.dNdxi_list(:,:,ob.i));
      end
      
      function value = get.d2Ndxi2(ob)
         value = squeeze(ob.d2Ndxi2_list(:,:,ob.i));
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(:,:,ob.i,ob.iel);
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
         value = ob.JxX * ob.J;
      end
      
      function value = get.dNdx(ob)
         value = ob.dNdX / ob.F;
      end
      
      function value = get.b(ob)
         value = ob.F*ob.F';
      end
      
      function value = get.B(ob)
         if (ob.finiteDisp)
            dx = ob.dNdx(:,1);	dy = ob.dNdx(:,2);
         else
            dx = ob.dNdX(:,1);	dy = ob.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0, dx(4),   0.0, dx(5),   0.0, dx(6),   0.0, dx(7),   0.0, dx(8),   0.0, dx(9),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3),   0.0, dy(4),   0.0, dy(5),   0.0, dy(6),   0.0, dy(7),   0.0, dy(8),   0.0, dy(9)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3), dy(4), dx(4), dy(5), dx(5), dy(6), dx(6), dy(7), dx(7), dy(8), dx(8), dy(9), dx(9)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2);
         value =[...
            dx(1)  0      dx(2)  0      dx(3)  0      dx(4)  0      dx(5)  0      dx(6)  0      dx(7)  0      dx(8)  0        dx(9)  0
            0      dy(1)  0      dy(2)  0      dy(3)  0      dy(4)  0      dy(5)  0      dy(6)  0      dy(7)  0      dy(8)    0      dy(9)
            dy(1)  dx(1)  dy(2)  dx(2)  dy(3)  dx(3)  dy(4)  dx(4)  dy(5)  dx(5)  dy(6)  dx(6)  dy(7)  dx(7)  dy(8)  dx(8)    dy(9)  dx(9)
            dy(1) -dx(1)  dy(2) -dx(2)  dy(3) -dx(3)  dy(4) -dx(4)  dy(5) -dx(5)  dy(6) -dx(6)  dy(7) -dx(7)  dy(8) -dx(8)    dy(9) -dx(9)];
      end
      
      function value = get.R(ob)
         [P, ~, Q] = svd(ob.F);
         value =  P*Q';
      end
      
      %% xi-dependant functions
      function ob = shapeIso(ob)
         x1 = ob.xi(1,:);  x2 = ob.xi(2,:);
         ndm = size(ob.xi,1);
         ngp = size(ob.xi,2);
         nen = 9;
         
         ob.Nmat = 1/4*[...
            x1 .* (x1-1) .*     x2 .* (x2-1)
            x1 .* (x1+1) .*     x2 .* (x2-1)
            x1 .* (x1+1) .*     x2 .* (x2+1)
            x1 .* (x1-1) .*     x2 .* (x2+1)
            -2 * (x1+1) .* (x1-1) .*     x2 .* (x2-1)
            -2 * (x1+1) .*     x1 .* (x2+1) .* (x2-1)
            -2 * (x1+1) .* (x1-1) .* (x2+1) .*     x2
            -2 *     x1 .* (x1-1) .* (x2+1) .* (x2-1)
            +4 * (x1+1) .* (x1-1) .* (x2+1) .* (x2-1)];
         if size(ob.Nmat, 1) == size(ob.Nmat, 2)
            ob.Ninv = inv(ob.Nmat);
         else
            ob.Ninv = 0;
         end
         %
         ob.dNdxi_list = zeros(nen, ngp, ndm);
         ob.dNdxi_list(:,:,1) = 1/4 *[...
            x2.*(2*x1-1).*(x2-1)
            x2.*(2*x1+1).*(x2-1)
            x2.*(2*x1+1).*(x2+1)
            x2.*(2*x1-1).*(x2+1)
            -4*x1.*x2.*(x2-1)
            -2*(2*x1+1).*(x2.^2-1)
            -4*x1.*x2.*(x2+1)
            -2*(2*x1-1).*(x2.^2-1)
            8*x1.*(x2.^2-1)];
         ob.dNdxi_list(:,:,2) = 1/4 *[...
            x1.*(2*x2-1).*(x1-1)
            x1.*(2*x2-1).*(x1+1)
            x1.*(2*x2+1).*(x1+1)
            x1.*(2*x2+1).*(x1-1)
            -2*(2*x2-1).*(x1.^2-1)
            -4*x1.*x2.*(x1+1)
            -2*(2*x2+1).*(x1.^2-1)
            -4*x1.*x2.*(x1-1)
            8*x2.*(x1.^2-1)];
         ob.dNdxi_list = permute(ob.dNdxi_list,[1 3 2]);
         %
         ob.d2Ndxi2_list = zeros(nen, ngp, 3);
         ob.d2Ndxi2_list(:,:,1) =1/4 *[2*x2.*(x2-1); 2*x2.*(x2-1); 2*x2.*(x2+1); 2*x2.*(x2+1); -4*x2.*(x2-1); 4-4*x2.^2; -4*x2.*(x2+1); 4-4*x2.^2; 8*x2.^2-8];
         ob.d2Ndxi2_list(:,:,2) =1/4 *[2*x1.*(x1-1); 2*x1.*(x1-1); 2*x1.*(x1+1); 2*x1.*(x1+1); -4*x1.*(x1-1); 4-4*x1.^2; -4*x1.*(x1+1); 4-4*x1.^2; 8*x1.^2-8];
         ob.d2Ndxi2_list(:,:,3) =1/4 *[...
            (2*x1-1).*(2*x2-1)
            (2*x1+1).*(2*x2-1)
            (2*x1+1).*(2*x2+1)
            (2*x1-1).*(2*x2+1)
            -4*x1.*(2*x2-1)
            -4*x2.*(2*x1+1)
            -4*x1.*(2*x2+1)
            -4*x2.*(2*x1-1)
            16*x1.*x2];
         ob.d2Ndxi2_list = permute(ob.d2Ndxi2_list,[1 3 2]);
      end
   end
end