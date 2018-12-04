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
         0.1381966011250105 0.5854101966249685 0.1381966011250105 0.1381966011250105 
         0.1381966011250105 0.1381966011250105 0.5854101966249685 0.1381966011250105
         0.1381966011250105 0.1381966011250105 0.1381966011250105 0.5854101966249685];
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
         val = 27*r*s*u;
      end
      
      function val = get.bubbB(ob)
         r = ob.xi(1,ob.i); s = ob.xi(2,ob.i); t = ob.xi(3,ob.i);
         u = 1-r-s-t;
         dbdxi = 27*[(s*u - r*s), (r*u - r*s), -r*s];
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
      
      %% xi-dependant functions
      function ob = shapeIso(ob)
         x1  = ob.xi(1,:); x2 = ob.xi(2,:); x3 = ob.xi(3,:);
         ndm = size(ob.xi,1);
         ngp = size(ob.xi,2);
         nen = 10;
   
         ob.Nmat = [...
            2*(1-x1-x2-x3).*((1-x1-x2-x3)-0.5)
            2*x1.*(x1-0.5)
            2*x2.*(x2-0.5)
            2*x3.*(x3-0.5)
            4*x1.*(1-x1-x2-x3)
            4*x1.*x2
            4*x2.*(1-x1-x2-x3)
            4*x3.*(1-x1-x2-x3)
            4*x1.*x3
            4*x2.*x3];
         
         ob.Ninv = ((ob.Nmat'*ob.Nmat)\ob.Nmat')';
         %
         vec = ones(1,ngp);
         ob.dNdxi_list = zeros(nen, ngp, ndm);
         ob.dNdxi_list(:,:,1) = [...
            4*(x1+x2+x3)-3; 4*x1-1; 0*vec; 0*vec; 4*(1-x2-x3-2*x1); 4*x2; -4*x2; -4*x3; 4*x3; 0*vec];
         ob.dNdxi_list(:,:,2) = [...
            4*(x1+x2+x3)-3; 0*vec; 4*x2-1; 0*vec; -4*x1; 4*x1; 4*(1-2*x2-x3-x1); -4*x3; 0*vec; 4*x3];
         ob.dNdxi_list(:,:,3) = [...
            4*(x1+x2+x3)-3; 0*vec; 0*vec; 4*x3-1; -4*x1; 0*vec; -4*x2; 4*(1-x2-2*x3-x1); 4*x1; 4*x2];
         ob.dNdxi_list = permute(ob.dNdxi_list,[1,3,2]);
         
         %
         ob.d2Ndxi2_list = zeros(nen,ngp,3);
         ob.d2Ndxi2_list = permute(ob.d2Ndxi2_list,[1,3,2]);
      end
   end
end