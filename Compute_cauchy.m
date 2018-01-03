function sigma_voigt = Compute_cauchy(ndof,gp,props)
   I = eye(ndof);  
   b = gp.b;
   J = gp.J;
   mu = props(1);
   lambda = props(2);

   sigma  = mu/J* (b-I) + lambda*(J-1)*I;
   sigma_voigt = [sigma(1,1); sigma(2,2); sigma(1,2)];
end
