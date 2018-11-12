function [N,dNdxi] = shapeF(xi,nen)
x1 = xi(1); x2 = xi(2);
switch nen
   case 3
      N = [...
         1-x1-x2
         x1
         x2];
      dNdxi =[...
         -1 -1
         +1  0
         +0  1];
   case 6
      N = [...
         (1-x1-x2)*(2*(1-x1-x2)-1)
         x1*(2*x1-1)
         x2*(2*x2-1)
         4*(1-x1-x2)*x1
         4*x1*x2
         4*x2*(1-x1-x2)];
      dNdxi = [...
         -3 + 4*(x1+x2)	, -3 + 4*(x1+x2)
         +4*x1-1        ,  0
         +0             ,  4*x2 - 1
         +4-8*x1-4*x2	, -4*x1
         +4*x2          ,  4*x1
         -4*x2          , -8*x2 + 4 - 4*x1];
   case 4
      N = 1/4*[...
         (1-x1).*(1-x2)
         (1+x1).*(1-x2)
         (1+x1).*(1+x2)
         (1-x1).*(1+x2)];
      dNdxi = 1/4*[...
         -(1-x2), -(1-x1)
         +(1-x2), -(1+x1)
         +(1+x2), +(1+x1)
         -(1+x2), +(1-x1)];
   case 9
         N = 1/4*[...
            x1 * (x1-1) *     x2 * (x2-1)
            x1 * (x1+1) *     x2 * (x2-1)
            x1 * (x1+1) *     x2 * (x2+1)
            x1 * (x1-1) *     x2 * (x2+1)
            -2 * (x1+1) * (x1-1) *     x2 * (x2-1)
            -2 * (x1+1) *     x1 * (x2+1) * (x2-1)
            -2 * (x1+1) * (x1-1) * (x2+1) *     x2
            -2 *     x1 * (x1-1) * (x2+1) * (x2-1)
            +4 * (x1+1) * (x1-1) * (x2+1) * (x2-1)];
         dNdxi = 1/4*[...
               (-1 + x2)*x2*x1 + (-1 + x2)*x2*(-1 + x1)                       , x2*(-1 + x1)*x1 + (-1 + x2)*(-1 + x1)*x1
               (-1 + x2)*x2*x1 + (-1 + x2)*x2*( 1 + x1)                       , x2*( 1 + x1)*x1 + (-1 + x2)*( 1 + x1)*x1
               ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*( 1 + x1)                       , x2*( 1 + x1)*x1 + ( 1 + x2)*( 1 + x1)*x1
               ( 1 + x2)*x2*x1 + ( 1 + x2)*x2*(-1 + x1)                       , x2*(-1 + x1)*x1 + ( 1 + x2)*(-1 + x1)*x1
               x2*(-1 + x2)*(-2 - 2*x1) - 2*x2*(-1 + x2)*(-1 + x1)            , x2*(-1 + x1)*(-2 - 2*x1) + (-1 + x2)*(-1 + x1)*(-2 - 2*x1)
               (-1 + x2)*(1 + x2)*(-2 - 2*x1) - 2*(-1 + x2)*(1 + x2)*x1       , (1 + x2)*x1*(-2 - 2*x1) + (-1 + x2)*x1*(-2 - 2*x1)
               x2*( 1 + x2)*(-2 - 2*x1) - 2*x2*( 1 + x2)*(-1 + x1)            , (1 + x2)*(-1 + x1)*(-2 - 2*x1) + x2*(-1 + x1)*(-2 - 2*x1)
               -2*(-1 + x2)*( 1 + x2)*x1 - 2*(-1 + x2)*(1 + x2)*(-1 + x1)     , -2*(1 + x2)*(-1 + x1)*x1 - 2*(-1 + x2)*(-1 + x1)*x1
               (-1 + x2)*(1 + x2)*(4 + 4*x1) + 4*(-1 + x2)*(1 + x2)*(-1 + x1) , (1 + x2)*(-1 + x1)*(4 + 4*x1) + (-1 + x2)*(-1 + x1)*(4 + 4*x1)];
   case 8
      x3 = xi(3);
      N = 1/8*[...
         (1-x1).*(1-x2).*(1-x3)
         (1+x1).*(1-x2).*(1-x3)
         (1+x1).*(1+x2).*(1-x3)
         (1-x1).*(1+x2).*(1-x3)
         (1-x1).*(1-x2).*(1+x3)
         (1+x1).*(1-x2).*(1+x3)
         (1+x1).*(1+x2).*(1+x3)
         (1-x1).*(1+x2).*(1+x3)];
      dNdxi = 1/8*[...
            -(x2- 1)*(x3- 1), -(x1- 1)*(x3- 1), -(x1- 1)*(x2- 1)
            +(x2- 1)*(x3- 1),  (x1+ 1)*(x3- 1),  (x1+ 1)*(x2- 1)
            -(x2+ 1)*(x3- 1), -(x1+ 1)*(x3- 1), -(x1+ 1)*(x2+ 1)
            +(x2+ 1)*(x3- 1),  (x1- 1)*(x3- 1),  (x1- 1)*(x2+ 1)
            +(x2- 1)*(x3+ 1),  (x1- 1)*(x3+ 1),  (x1- 1)*(x2- 1)
            -(x2- 1)*(x3+ 1), -(x1+ 1)*(x3+ 1), -(x1+ 1)*(x2- 1)
            +(x2+ 1)*(x3+ 1),  (x1+ 1)*(x3+ 1),  (x1+ 1)*(x2+ 1)
            -(x2+ 1)*(x3+ 1), -(x1- 1)*(x3+ 1), -(x1- 1)*(x2+ 1)];
   otherwise
      error('unimplemented shape function');
end
