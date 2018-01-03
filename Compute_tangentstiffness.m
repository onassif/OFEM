function [D, ctan] = Compute_tangentstiffness(gp, props)

   c = zeros(3,3,3,3);
   I = eye(3);
   mu = props(1);
   lambda = props(2);
   J = gp.J;
   
   for i=1:3
       for j=1:3
           for k=1:3
               for l=1:3
                   c(i,j,k,l) = mu/J*( I(i,k)*I(j,l) + I(i,l)*I(j,k) )+...
                       lambda*(  (2*J-1)*I(i,j)*I(k,l) - (J-1)*( I(i,k)*I(j,l) + I(i,l)*I(j,k) )  );
               end
           end
       end
   end
   
 D = [ c(1,1,1,1) c(1,1,2,2) c(1,1,1,2)
       c(2,2,1,1) c(2,2,2,2) c(2,2,1,2)
       c(1,2,1,1) c(1,2,2,2) c(1,2,1,2)];
 ctan = c;
end
