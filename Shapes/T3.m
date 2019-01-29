classdef T3
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
      xi = 1/3 .*[1; 1];
      weights = 0.5;
      Nmat
      dNdxi_list
      d2Ndxi2_list
      
      dNdX_list
      dXdxi_list
      det_dXdxi_list
   end
   
   properties (SetAccess = private)
      finiteDisp
      bubb
      bubbB
      Ninv
      dNdX
      dNdx
      dNdxi
      d2Ndxi2
      dXdxi
      B
      Bf
      NN
      NN_2
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
      function ob = T3(varargin)
         ob.finiteDisp = varargin{1};
      end
      %% Get functions
      function value = get.N(ob)
         value = ob.Nmat(:,ob.i);
      end
      
      function value = get.dNdxi(ob)
         value = ob.dNdxi_list;
      end
      
      function value = get.d2Ndxi2(ob)
         value = ob.d2Ndxi2_list;
      end
      
      function value = get.dXdxi(ob)
         value = ob.dXdxi_list(:,:,ob.iel);
      end
      
      function value = get.w(ob)
         value = ob.weights(ob.i);
      end
      
      function value = get.J(ob)
         value = ob.det_dXdxi_list(ob.iel);
      end
      
      function value = get.dNdX(ob)
         value = ob.dNdX_list(:,:,ob.iel);
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
            dx = ob.dNdx(:,1);  dy = ob.dNdx(:,2);
         else
            dx = ob.dNdX(:,1);  dy = ob.dNdX(:,2);
         end
         value=[...
            dx(1),   0.0, dx(2),   0.0, dx(3),   0.0
            0.0  , dy(1),   0.0, dy(2),   0.0, dy(3)
            dy(1), dx(1), dy(2), dx(2), dy(3), dx(3)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2);
         value =[...
            dx(1)  0      dx(2)  0      dx(3)  0
            0      dy(1)  0      dy(2)  0      dy(3)
            dy(1)  dx(1)  dy(2)  dx(2)  dy(3)  dx(3)
            dy(1) -dx(1)  dy(2) -dx(2)  dy(3) -dx(3)];
      end
      
      function value = get.R(ob)
         [P, ~, Q] = svd(ob.F);
         value =  P*Q';
      end
      
      function val = get.bubb(ob)
         r = ob.xi(1,:); s = ob.xi(2,:);
         val = 4*(1-r-s).*r;
      end
      
      function val = get.bubbB(ob)
         r = ob.xi(1,ob.i); s = ob.xi(2,ob.i);
         dbdxi = 4*[(1-2*r-s), -r];
         if (ob.finiteDisp)
            dxdxi = ob.F*ob.dXdxi;
            dbdx  = dbdxi / dxdxi;
            val = [...
               dbdx(1) 0       dbdx(2)  dbdx(2)
               0       dbdx(2) dbdx(1) -dbdx(1)]';
         else
            dbdX  = dbdxi / ob.dXdxi;
            val = [...
               dbdX(1) 0       dbdX(2)
               0       dbdX(2) dbdX(1)]';
         end
      end
      
      function val = get.NN(ob)
         N1 = ob.N(1);  N2 = ob.N(2);  N3 = ob.N(3);
         
         val = [...
            N1*N1	0     N1*N2	0     N1*N3	0     
            0   	N1*N1	0     N1*N2	0     N1*N3
            N2*N1	0     N2*N2	0     N2*N3	0    
            0   	N2*N1	0    	N2*N2	0    	N2*N3
            N3*N1	0     N3*N2	0     N3*N3	0    
            0   	N3*N1	0    	N3*N2	0    	N3*N3];
      end
      
      function val = get.NN_2(ob)
         N1 = ob.N(1);  N2 = ob.N(2);  N3 = ob.N(3);
         
         val = diag([...
            N1*N1 + N1*N2 + N1*N3     
            N1*N1 + N1*N2 + N1*N3
            N2*N1 + N2*N2 + N2*N3
            N2*N1 + N2*N2 + N2*N3
            N3*N1 + N3*N2 + N3*N3  
            N3*N1 + N3*N2 + N3*N3]);
      end      
      %% xi-dependant functions
      function ob = shapeIso(ob,varargin)
         ndm = size(ob.xi,1); 
         ngp = size(ob.xi,2);
         nen = 3;
         x1  = ob.xi(1,:); x2 = ob.xi(2,:);

         ob.Nmat = [1-x1-x2; x1; x2];
         ob.Ninv = ((ob.Nmat'*ob.Nmat)\ob.Nmat')';
         %
         vec = ones(1,ngp);
         ob.dNdxi_list = zeros(nen,ngp,ndm);
         ob.dNdxi_list(:,:,1) =[-1*vec; +1*vec; 0*vec];
         ob.dNdxi_list(:,:,2) =[-1*vec;  0*vec; 1*vec];
         ob.dNdxi_list = permute(ob.dNdxi_list,[1,3,2]);

         %
         ob.d2Ndxi2_list = zeros(nen,ngp,3);
         ob.d2Ndxi2_list = permute(ob.d2Ndxi2_list,[1,3,2]);
      end
   end
end