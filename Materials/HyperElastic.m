classdef HyperElastic
   %HyperNeo Computes Cauchy stress based on NeoHookean model
   
   properties (SetAccess = private)
      ndm;
      ndof;
      ngp;
      finiteDisp = 1;
      
      lambda;
      shear;
      tau;
      
      I
      I4_dev
      I4_bulk
   end
   %%
   methods
      %% Construct
      function ob = HyperElastic(num, props, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.ngp  = num.gp;
         
         if strcmp(props{1,1} ,'nu') && strcmp(props{2,1} ,'E')
            nu  = props{1,2};
            E  = props{2,2};
         elseif strcmp(props{1,1} ,'E') && strcmp(props{2,1} ,'nu')
            E  = props{1,2};
            nu  = props{2,2};
         else
            error(['You''ve chosen Neo Hookean material but specified',...
               ' incompatible material properties, I''m disapponted']);
         end
         ob.lambda = nu*E/((1+nu)*(1-2*nu));
         ob.shear  = E/(2*(1+nu));
         
         ob.I      = eye(ob.ndm);
         ob.I4_dev = identity.I4_dev;
         ob.I4_bulk= identity.I4_bulk;
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, gp, el, ~)
         eps = gp.B * el.Uvc;
      end
      %% Sigma
      function [sigma_voigt, ob] = computeCauchy(ob, gp, ~, ~)
         lam = ob.lambda;
         mu  = ob.shear;
         
         JxX = det(gp.F);
         sigma = mu*(gp.b - ob.I) + JxX*lam*(JxX-1)*ob.I;
         
         if ob.ndm == 2
            sigma_voigt = [sigma(1,1) sigma(2,2) sigma(1,2)]';
         elseif ob.ndm == 3
            sigma_voigt = [sigma(1,1) sigma(2,2) sigma(3,3) sigma(1,2) sigma(2,3) sigma(3,1)]';
         end
      end
      %% Tangential stiffness
      function [D, ctan, ob] = computeTangentStiffness(ob, gp, ~, ~)
         lam = ob.lambda;
         mu  = ob.shear;
         
         matE = diag([2,2,2,1,1,1]);
         JxX = det(gp.F);
         
         D = mu*matE + JxX*lam*( (2*JxX-1)*ob.I4_bulk - (JxX-1)*matE );
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if ob.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         % Definitions
         B=gp.Bf;
         D=gp.D;
         S = gp.sigma;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         if ob.ndm == 2
            Dmat = [D, zeros(3,1); zeros(1,4)];
            Dgeo = [...
               +S(1)      0        1/2*S(3)         1/2*S(3)
               +0         S(2)     1/2*S(3)        -1/2*S(3)
               +1/2*S(3)  1/2*S(3) 1/4*(S(1)+S(2))  1/4*(S(2)-S(1))
               +1/2*S(3) -1/2*S(3) 1/4*(S(2)-S(1))  1/4*(S(1)+S(2))];
         elseif ob.ndm == 3
            Dmat = [D, zeros(6,3); zeros(3,9)];
            Dgeo = [...
               +S(1)      0         0         1/2*S(4)         0                1/2*S(6)         1/2*S(4)         0               -1/2*S(6)
               +0         S(2)      0         1/2*S(4)         1/2*S(5)         0               -1/2*S(4)         1/2*S(5)         0
               +0         0         S(3)      0                1/2*S(5)         1/2*S(6)         0               -1/2*S(5)         1/2*S(6)
               +1/2*S(4)  1/2*S(4)  0         1/4*(S(1)+S(2))  1/4*S(6)         1/4*S(5)         1/4*(S(2)-S(1))  1/4*S(6)        -1/4*S(5)
               +0         1/2*S(5)  1/2*S(5)  1/4*S(6)         1/4*(S(2)+S(3))  1/4*S(4)        -1/4*S(6)         1/4*(S(3)-S(2))  1/4*S(4)
               +1/2*S(6)  0         1/2*S(6)  1/4*S(5)         1/4*S(4)         1/4*(S(1)+S(3))  1/4*S(5)        -1/4*S(4)         1/4*(S(1)-S(3))
               +1/2*S(4) -1/2*S(4)  0         1/4*(S(2)-S(1)) -1/4*S(6)         1/4*S(5)         1/4*(S(1)+S(2)) -1/4*S(6)        -1/4*S(5)
               +0         1/2*S(5) -1/2*S(5)  1/4*S(6)         1/4*(S(3)-S(2)) -1/4*S(4)        -1/4*S(6)         1/4*(S(2)+S(3)) -1/4*S(4)
               -1/2*S(6)  0         1/2*S(6) -1/4*S(5)         1/4*S(4)         1/4*(S(1)-S(3)) -1/4*S(5)        -1/4*S(4)         1/4*(S(1)+S(3))];
         end
         if gp.i == 1
            Kel = B'*(Dmat+Dgeo)*B *gp.J*gp.w;
         else
            Kel = el.K + B'*(Dmat+Dgeo)*B *gp.J*gp.w;
         end
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, ~)
         if gp.i == 1
            Fint = (gp.B'*gp.sigma) *gp.J *gp.w;
         else
            Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
         end
      end
   end
end

