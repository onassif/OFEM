classdef T10
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
      xi = [...
         0.1381966011 0.5854101966 0.1381966011 0.1381966011 
         0.1381966011 0.1381966011 0.5854101966 0.1381966011
         0.1381966011 0.1381966011 0.1381966011 0.5854101966];
      weights = 1/24*[1 1 1 1]';
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
      function ob = T10(varargin)
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
         value = ob.det_dXdxi_list(ob.i,ob.iel);
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
            dx = ob.dNdx(:,1); dy = ob.dNdx(:,2); dz = ob.dNdx(:,3);
         else
            dx = ob.dNdX(:,1); dy = ob.dNdX(:,2); dz = ob.dNdX(:,3);
         end
         value =[...
            dx(1) 0     0     dx(2) 0     0     dx(3) 0     0     dx(4) 0     0     dx(5) 0     0     dx(6) 0     0     dx(7) 0     0     dx(8) 0     0     dx(9) 0     0     dx(10) 0      0    
            0     dy(1) 0     0     dy(2) 0     0     dy(3) 0     0     dy(4) 0     0     dy(5) 0     0     dy(6) 0     0     dy(7) 0     0     dy(8) 0     0     dy(9) 0     0      dy(10) 0 
            0     0     dz(1) 0     0     dz(2) 0     0     dz(3) 0     0     dz(4) 0     0     dz(5) 0     0     dz(6) 0     0     dz(7) 0     0     dz(8) 0     0     dz(9) 0      0      dz(10)
            dy(1) dx(1) 0     dy(2) dx(2) 0     dy(3) dx(3) 0     dy(4) dx(4) 0     dy(5) dx(5) 0     dy(6) dx(6) 0     dy(7) dx(7) 0     dy(8) dx(8) 0     dy(9) dx(9) 0     dy(10) dx(10) 0
            0     dz(1) dy(1) 0     dz(2) dy(2) 0     dz(3) dy(3) 0     dz(4) dy(4) 0     dz(5) dy(5) 0     dz(6) dy(6) 0     dz(7) dy(7) 0     dz(8) dy(8) 0     dz(9) dy(9) 0      dz(10) dy(10)
            dz(1) 0     dx(1) dz(2) 0     dx(2) dz(3) 0     dx(3) dz(4) 0     dx(4) dz(5) 0     dx(5) dz(6) 0     dx(6) dz(7) 0     dx(7) dz(8) 0     dx(8) dz(9) 0     dx(9) dz(10) 0      dx(10)];
      end
      
      function value = get.Bf(ob)
         dx = ob.dNdx(:,1); dy = ob.dNdx(:,2); dz = ob.dNdx(:,3);
         value =[...
            dx(1)  0      0      dx(2)  0      0      dx(3)  0      0      dx(4)  0      0     dx(5)  0      0      dx(6)  0      0      dx(7)  0      0      dx(8)  0      0      dx(9)  0      0      dx(10)  0       0      
            0      dy(1)  0      0      dy(2)  0      0      dy(3)  0      0      dy(4)  0     0      dy(5)  0      0      dy(6)  0      0      dy(7)  0      0      dy(8)  0      0      dy(9)  0      0       dy(10)  0    
            0      0      dz(1)  0      0      dz(2)  0      0      dz(3)  0      0      dz(4) 0      0      dz(5)  0      0      dz(6)  0      0      dz(7)  0      0      dz(8)  0      0      dz(9)  0       0       dz(10)
            dy(1)  dx(1)  0      dy(2)  dx(2)  0      dy(3)  dx(3)  0      dy(4)  dx(4)  0     dy(5)  dx(5)  0      dy(6)  dx(6)  0      dy(7)  dx(7)  0      dy(8)  dx(8)  0      dy(9)  dx(9)  0      dy(10)  dx(10)  0    
            0      dz(1)  dy(1)  0      dz(2)  dy(2)  0      dz(3)  dy(3)  0      dz(4)  dy(4) 0      dz(5)  dy(5)  0      dz(6)  dy(6)  0      dz(7)  dy(7)  0      dz(8)  dy(8)  0      dz(9)  dy(9)  0       dz(10)  dy(10)
            dz(1)  0      dx(1)  dz(2)  0      dx(2)  dz(3)  0      dx(3)  dz(4)  0      dx(4) dz(5)  0      dx(5)  dz(6)  0      dx(6)  dz(7)  0      dx(7)  dz(8)  0      dx(8)  dz(9)  0      dx(9)  dz(10)  0       dx(10)
            dy(1) -dx(1)  0      dy(2) -dx(2)  0      dy(3) -dx(3)  0      dy(4) -dx(4)  0     dy(5) -dx(5)  0      dy(6) -dx(6)  0      dy(7) -dx(7)  0      dy(8) -dx(8)  0      dy(9) -dx(9)  0      dy(10) -dx(10)  0 
            0      dz(1) -dy(1)  0      dz(2) -dy(2)  0      dz(3) -dy(3)  0      dz(4) -dy(4) 0      dz(5) -dy(5)  0      dz(6) -dy(6)  0      dz(7) -dy(7)  0      dz(8) -dy(8)  0      dz(9) -dy(9)  0       dz(10) -dy(10)
            -dz(1) 0      dx(1) -dz(2)  0      dx(2) -dz(3)  0      dx(3) -dz(4)  0      dx(4) -dz(5) 0      dx(5) -dz(6)  0      dx(6) -dz(7)  0      dx(7) -dz(8)  0      dx(8) -dz(9)  0      dx(9) -dz(10)  0       dx(10)];
      end
      
      function value = get.R(ob)
         [P, ~, Q] = svd(ob.F);
         value =  P*Q';
      end
      
      function val = get.bubb(ob)
         r = ob.xi(1,:); s = ob.xi(2,:); t = ob.xi(3,:);
         u = 1-r-s-t;
         val = (1-(2*r-1).^4).*(1-(2*s-1).^4).*(1-(2*u-1).^4);
      end
      
      function val = get.bubbB(ob)
         r = ob.xi(1,ob.i); s = ob.xi(2,ob.i); t = ob.xi(3,ob.i);
         dbdxi = [...
            - 8*(2*r-1)^3*((2*(r+s+t)-1)^4-1)*((2*s-1)^4-1) - 8*((2*r-1)^4-1)*((2*s-1)^4-1)*(2*(r+s+t)-1)^3,...
            - 8*(2*s-1)^3*((2*(r+s+t)-1)^4-1)*((2*r-1)^4-1) - 8*((2*r-1)^4-1)*((2*s-1)^4-1)*(2*(r+s+t)-1)^3,...
            - 8*((2*r-1)^4-1)*((2*s-1)^4-1)*(2*(r+s+t)-1)^3];
         if (ob.finiteDisp)
            dxdxi = ob.F*ob.dXdxi;
            dbdx  = dbdxi / dxdxi;
            val   = [...
               dbdx(1) 0       0       dbdx(2) 0       dbdx(3)  dbdx(2)  0       -dbdx(3)
               0       dbdx(2) 0       dbdx(1) dbdx(3) 0       -dbdx(1)  dbdx(3)  0
               0       0       dbdx(3) 0       dbdx(2) dbdx(1)  0       -dbdx(2)  dbdx(1)]';
         else
            dbdX = dbdxi / ob.dXdxi;
            val  = [...
               dbdX(1) 0       0       dbdX(2) 0       dbdX(3)
               0       dbdX(2) 0       dbdX(1) dbdX(3) 0
               0       0       dbdX(3) 0       dbdX(2) dbdX(1)]';
         end
      end

      function val = get.NN(ob)
         N1 = ob.N(1);  N2 = ob.N(2);  N3 = ob.N(3);  N4 = ob.N(4);	N5  = ob.N(5);
         N6 = ob.N(6);  N7 = ob.N(7);  N8 = ob.N(8);  N9 = ob.N(9);  N10 = ob.N(10);
         
         val = [...
            N1*N1	 0      0      N1*N2	 0      0      N1*N3	 0      0      N1*N4	 0      0      N1*N5	 0      0      N1*N6	 0      0      N1*N7	 0      0      N1*N8 	0    0      N1*N9	 0      0      N1*N10	0        0
            0   	 N1*N1  0    	0      N1*N2  0    	0      N1*N3  0    	0      N1*N4  0      0      N1*N5  0      0      N1*N6  0      0      N1*N7  0      0      N1*N8  0      0      N1*N9  0      0        N1*N10	0
            0      0      N1*N1	0    	 0      N1*N2	0    	 0      N1*N3	0    	 0      N1*N4	0      0      N1*N5	0      0      N1*N6	0      0      N1*N7	0      0      N1*N8	0      0      N1*N9	0        0        N1*N10
            N2*N1	 0      0      N2*N2	 0      0      N2*N3	 0      0      N2*N4	 0      0      N2*N5	 0      0      N2*N6	 0      0      N2*N7	 0      0      N2*N8	 0      0      N2*N9	 0      0      N2*N10	0        0
            0   	 N2*N1  0    	0      N2*N2  0    	0      N2*N3  0    	0      N2*N4  0      0      N2*N5  0      0      N2*N6  0      0      N2*N7  0      0      N2*N8  0      0      N2*N9  0      0        N2*N10	0
            0      0      N2*N1	0    	 0      N2*N2	0    	 0      N2*N3	0    	 0      N2*N4	0      0      N2*N5	0      0      N2*N6	0      0      N2*N7	0      0      N2*N8	0      0      N2*N9	0        0        N2*N10
            N3*N1	 0      0      N3*N2	 0      0      N3*N3	 0      0      N3*N4	 0      0      N3*N5	 0      0      N3*N6	 0      0      N3*N7	 0      0      N3*N8	 0      0      N3*N9	 0      0      N3*N10	0        0
            0   	 N3*N1  0    	0      N3*N2  0    	0      N3*N3  0    	0      N3*N4  0      0      N3*N5  0      0      N3*N6  0      0      N3*N7  0      0      N3*N8  0      0      N3*N9  0      0        N3*N10	0
            0      0      N3*N1	0    	 0      N3*N2	0    	 0      N3*N3	0    	 0      N3*N4	0      0      N3*N5	0      0      N3*N6	0      0      N3*N7	0      0      N3*N8	0      0      N3*N9	0        0        N3*N10
            N4*N1	 0      0      N4*N2	 0      0      N4*N3	 0      0      N4*N4	 0      0      N4*N5	 0      0      N4*N6	 0      0      N4*N7	 0      0      N4*N8	 0      0      N4*N9	 0      0      N4*N10	0        0
            0   	 N4*N1  0    	0      N4*N2  0    	0      N4*N3  0    	0      N4*N4  0      0      N4*N5  0      0      N4*N6  0      0      N4*N7  0      0      N4*N8  0      0      N4*N9  0      0        N4*N10	0
            0      0      N4*N1	0    	 0      N4*N2	0    	 0      N4*N3	0    	 0      N4*N4	0      0      N4*N5	0      0      N4*N6	0      0      N4*N7	0      0      N4*N8	0      0      N4*N9	0        0        N4*N10
            N5*N1	 0      0      N5*N2	 0      0      N5*N3	 0      0      N5*N4	 0      0      N5*N5	 0      0      N5*N6	 0      0      N5*N7	 0      0      N5*N8	 0      0      N5*N9	 0      0      N5*N10	0        0
            0   	 N5*N1  0    	0      N5*N2  0    	0      N5*N3  0    	0      N5*N4  0      0      N5*N5  0      0      N5*N6  0      0      N5*N7  0      0      N5*N8  0      0      N5*N9  0      0        N5*N10	0
            0      0      N5*N1	0    	 0      N5*N2	0    	 0      N5*N3	0    	 0      N5*N4	0      0      N5*N5	0      0      N5*N6	0      0      N5*N7	0      0      N5*N8	0      0      N5*N9	0        0        N5*N10
            N6*N1	 0      0      N6*N2	 0      0      N6*N3	 0      0      N6*N4	 0      0      N6*N5	 0      0      N6*N6	 0      0      N6*N7	 0      0      N6*N8	 0      0      N6*N9	 0      0      N6*N10	0        0
            0   	 N6*N1  0    	0      N6*N2  0    	0      N6*N3  0    	0      N6*N4  0      0      N6*N5  0      0      N6*N6  0      0      N6*N7  0      0      N6*N8  0      0      N6*N9  0      0        N6*N10	0
            0      0      N6*N1	0    	 0      N6*N2	0    	 0      N6*N3	0    	 0      N6*N4	0      0      N6*N5	0      0      N6*N6	0      0      N6*N7	0      0      N6*N8	0      0      N6*N9	0        0        N6*N10
            N7*N1	 0      0      N7*N2	 0      0      N7*N3	 0      0      N7*N4	 0      0      N7*N5	 0      0      N7*N6	 0      0      N7*N7	 0      0      N7*N8	 0      0      N7*N9	 0      0      N7*N10	0        0
            0   	 N7*N1  0    	0      N7*N2  0    	0      N7*N3  0    	0      N7*N4  0      0      N7*N5  0      0      N7*N6  0      0      N7*N7  0      0      N7*N8  0      0      N7*N9  0      0        N7*N10	0
            0      0      N7*N1	0    	 0      N7*N2	0    	 0      N7*N3	0    	 0      N7*N4	0      0      N7*N5	0      0      N7*N6	0      0      N7*N7	0      0      N7*N8	0      0      N7*N9	0        0        N7*N10
            N8*N1	 0      0      N8*N2	 0      0      N8*N3	 0      0      N8*N4	 0      0      N8*N5	 0      0      N8*N6	 0      0      N8*N7	 0      0      N8*N8	 0      0      N8*N9	 0      0      N8*N10	0        0
            0   	 N8*N1  0    	0      N8*N2  0    	0      N8*N3  0    	0      N8*N4  0      0      N8*N5  0      0      N8*N6  0      0      N8*N7  0      0      N8*N8  0      0      N8*N9  0      0        N8*N10	0
            0      0      N8*N1	0    	 0      N8*N2	0    	 0      N8*N3	0    	 0      N8*N4	0      0      N8*N5	0      0      N8*N6	0      0      N8*N7	0      0      N8*N8	0      0      N8*N9	0        0        N8*N10
            N9*N1  0      0      N9*N2	 0      0      N9*N3	 0      0      N9*N4	 0      0      N9*N5	 0      0      N9*N6	 0      0      N9*N7	 0      0      N9*N8	 0      0      N9*N9	 0      0      N9*N10	0        0
            0      N9*N1  0    	0      N9*N2  0    	0      N9*N3  0    	0      N9*N4  0      0      N9*N5  0      0      N9*N6  0      0      N9*N7  0      0      N9*N8  0      0      N9*N9  0      0        N9*N10	0
            0      0      N9*N1	0    	 0      N9*N2	0    	 0      N9*N3	0    	 0      N9*N4	0      0      N9*N5	0      0      N9*N6	0      0      N9*N7	0      0      N9*N8	0      0      N9*N9	0        0        N9*N10
            N10*N1 0      0      N10*N2 0      0      N10*N3 0      0      N10*N4 0      0      N10*N5 0      0      N10*N6 0      0      N10*N7 0      0      N10*N8 0      0      N10*N9 0      0      N10*N10	0        0
            0      N10*N1 0    	0      N10*N2 0    	0      N10*N3 0    	0      N10*N4 0      0      N10*N5 0      0      N10*N6 0      0      N10*N7 0      0      N10*N8 0      0      N10*N9 0      0        N10*N10	0
            0      0      N10*N1	0    	 0      N10*N2	0    	 0      N10*N3	0    	 0      N10*N4	0      0      N10*N5	0      0      N10*N6	0      0      N10*N7	0      0      N10*N8	0      0      N10*N9	0        0        N10*N10];
      end

      function val = get.NN_2(ob)
         N1 = ob.N(1);  N2 = ob.N(2);  N3 = ob.N(3);  N4 = ob.N(4);	N5  = ob.N(5);
         N6 = ob.N(6);  N7 = ob.N(7);  N8 = ob.N(8);  N9 = ob.N(9);  N10 = ob.N(10);
         
         val = diag([...
            N1*N1 + N1*N2 + N1*N3 + N1*N4 + N1*N5 + N1*N6 + N1*N7 + N1*N8 + N1*N9 + N1*N10
            N1*N1 + N1*N2 + N1*N3 + N1*N4 + N1*N5 + N1*N6 + N1*N7 + N1*N8 + N1*N9 + N1*N10
            N1*N1 + N1*N2 + N1*N3 + N1*N4 + N1*N5 + N1*N6 + N1*N7 + N1*N8 + N1*N9 + N1*N10
            N2*N1 + N2*N2 + N2*N3 + N2*N4 + N2*N5 + N2*N6 + N2*N7 + N2*N8 + N2*N9 + N2*N10
            N2*N1 + N2*N2 + N2*N3 + N2*N4 + N2*N5 + N2*N6 + N2*N7 + N2*N8 + N2*N9 + N2*N10
            N2*N1 + N2*N2 + N2*N3 + N2*N4 + N2*N5 + N2*N6 + N2*N7 + N2*N8 + N2*N9 + N2*N10
            N3*N1 + N3*N2 + N3*N3 + N3*N4 + N3*N5 + N3*N6 + N3*N7 + N3*N8 + N3*N9 + N3*N10
            N3*N1 + N3*N2 + N3*N3 + N3*N4 + N3*N5 + N3*N6 + N3*N7 + N3*N8 + N3*N9 + N3*N10
            N3*N1 + N3*N2 + N3*N3 + N3*N4 + N3*N5 + N3*N6 + N3*N7 + N3*N8 + N3*N9 + N3*N10
            N4*N1 + N4*N2 + N4*N3 + N4*N4 + N4*N5 + N4*N6 + N4*N7 + N4*N8 + N4*N9 + N4*N10
            N4*N1 + N4*N2 + N4*N3 + N4*N4 + N4*N5 + N4*N6 + N4*N7 + N4*N8 + N4*N9 + N4*N10
            N4*N1 + N4*N2 + N4*N3 + N4*N4 + N4*N5 + N4*N6 + N4*N7 + N4*N8 + N4*N9 + N4*N10
            N5*N1 + N5*N2 + N5*N3 + N5*N4 + N5*N5 + N5*N6 + N5*N7 + N5*N8 + N5*N9 + N5*N10
            N5*N1 + N5*N2 + N5*N3 + N5*N4 + N5*N5 + N5*N6 + N5*N7 + N5*N8 + N5*N9 + N5*N10
            N5*N1 + N5*N2 + N5*N3 + N5*N4 + N5*N5 + N5*N6 + N5*N7 + N5*N8 + N5*N9 + N5*N10
            N6*N1 + N6*N2 + N6*N3 + N6*N4 + N6*N5 + N6*N6 + N6*N7 + N6*N8 + N6*N9 + N6*N10
            N6*N1 + N6*N2 + N6*N3 + N6*N4 + N6*N5 + N6*N6 + N6*N7 + N6*N8 + N6*N9 + N6*N10
            N6*N1 + N6*N2 + N6*N3 + N6*N4 + N6*N5 + N6*N6 + N6*N7 + N6*N8 + N6*N9 + N6*N10
            N7*N1 + N7*N2 + N7*N3 + N7*N4 + N7*N5 + N7*N6 + N7*N7 + N7*N8 + N7*N9 + N7*N10
            N7*N1 + N7*N2 + N7*N3 + N7*N4 + N7*N5 + N7*N6 + N7*N7 + N7*N8 + N7*N9 + N7*N10
            N7*N1 + N7*N2 + N7*N3 + N7*N4 + N7*N5 + N7*N6 + N7*N7 + N7*N8 + N7*N9 + N7*N10
            N8*N1 + N8*N2 + N8*N3 + N8*N4 + N8*N5 + N8*N6 + N8*N7 + N8*N8 + N8*N9 + N8*N10
            N8*N1 + N8*N2 + N8*N3 + N8*N4 + N8*N5 + N8*N6 + N8*N7 + N8*N8 + N8*N9 + N8*N10
            N8*N1 + N8*N2 + N8*N3 + N8*N4 + N8*N5 + N8*N6 + N8*N7 + N8*N8 + N8*N9 + N8*N10
            N9*N1 + N9*N2 + N9*N3 + N9*N4 + N9*N5 + N9*N6 + N9*N7 + N9*N8 + N9*N9 + N9*N10
            N9*N1 + N9*N2 + N9*N3 + N9*N4 + N9*N5 + N9*N6 + N9*N7 + N9*N8 + N9*N9 + N9*N10
            N9*N1 + N9*N2 + N9*N3 + N9*N4 + N9*N5 + N9*N6 + N9*N7 + N9*N8 + N9*N9 + N9*N10
            N10*N1 + N10*N2 + N10*N3 + N10*N4 + N10*N5 + N10*N6 + N10*N7 + N10*N8 + N10*N9 + N10*N10
            N10*N1 + N10*N2 + N10*N3 + N10*N4 + N10*N5 + N10*N6 + N10*N7 + N10*N8 + N10*N9 + N10*N10
            N10*N1 + N10*N2 + N10*N3 + N10*N4 + N10*N5 + N10*N6 + N10*N7 + N10*N8 + N10*N9 + N10*N10]);
      end      
      %% xi-dependant functions
      function ob = shapeIso(ob)
         x1  = ob.xi(1,:); x2 = ob.xi(2,:); x3 = ob.xi(3,:);
         x0  = 1 - x1 - x2 - x3; 
         ndm = size(ob.xi,1);
         ngp = size(ob.xi,2);
         nen = 10;
   
         ob.Nmat = [...
            2*x0.*(x0-0.5)
            2*x1.*(x1-0.5)
            2*x2.*(x2-0.5)
            2*x3.*(x3-0.5)
            4*x1.*x0
            4*x1.*x2
            4*x2.*x0
            4*x3.*x0
            4*x1.*x3
            4*x2.*x3];
         if size(ob.Nmat,2) ~= 13 && size(ob.Nmat,2) ~= 24
            ob.Ninv = ((ob.Nmat'*ob.Nmat)\ob.Nmat')';
         else
            ob.Ninv = 0;
         end
         %
         vec = ones(1,ngp);
         ob.dNdxi_list = zeros(nen, ngp, ndm);
         ob.dNdxi_list(:,:,1) = 4*[1/4-x0; x1-1/4; 0*vec; 0*vec; x0-x1; x2; -x2; -x3; x3; 0*vec];
         ob.dNdxi_list(:,:,2) = 4*[1/4-x0; 0*vec; x2-1/4; 0*vec; -x1; x1; x0-x2; -x3; 0*vec; x3];
         ob.dNdxi_list(:,:,3) = 4*[1/4-x0; 0*vec; 0*vec; x3-1/4; -x1; 0*vec; -x2; x0-x3; x1; x2];
         ob.dNdxi_list = permute(ob.dNdxi_list,[1,3,2]);
         
         %
         ob.d2Ndxi2_list = zeros(nen,ngp,3);
         ob.d2Ndxi2_list = permute(ob.d2Ndxi2_list,[1,3,2]);
      end
   end
end