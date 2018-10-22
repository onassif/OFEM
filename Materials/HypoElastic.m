classdef HypoElastic
   %Simple Hypo elastic routine
   %   Currently only if ndm=2, ndof=2
   
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      ngp;
      finiteDisp = 0;
      
      sig0;
      eps0;
      nExp;
      nu;
      K;
      
      I
      I4_dev
      I4_bulk
      
      e
   end
   %%
   methods
      %% Construct
      function ob = HypoElastic(num, props, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.ngp  = num.gp;
         
         for i=1:length(props)
            switch props{i,1}
               case 'sig_0'
                  ob.sig0 = props{i,2};
               case 'eps_0'
                  ob.eps0 = props{i,2};
               case 'n'
                  ob.nExp = props{i,2};
               case 'nu'
                  ob.nu   = props{i,2};
               otherwise
                  error(['You''ve chosen Hypo Elastic material but specified',...
                     ' incompatible material properties, I''m disapponted']);
            end
         end
         
         ob.K = (ob.nExp*ob.sig0/ob.eps0) / (3*(1-2*ob.nu));
         
         ob.I      = eye(num.ndm);
         ob.I4_dev = identity.I4_dev;
         ob.I4_bulk= identity.I4_bulk;
         
      end
      %% Epsilon
      function [eps, ob] = Strain(ob, gp, el, ~)
         eps  = gp.B * el.Uvc;
         ob.e = T1T2(eps,2) - (1/3)*trace(T1T2(eps,2))*ob.I;
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, ~, ~)
         K    = ob.K;
         sig0 = ob.sig0;
         eps0 = ob.eps0;
         nExp = ob.nExp;
         
         e = ob.e;
         eEff = sqrt( (2/3)*sum(dot(e,e)) );
         sEff = ob.effectivestress(eEff, sig0, eps0, nExp);
         dsde = ob.hardeningslope (eEff, sig0, eps0, nExp);
         
         if     ob.ndm == 2
            e = [e(1,1) e(2,2) 0      e(1,2) 0      0     ]';
         elseif ob.ndm == 3
            e = [e(1,1) e(2,2) e(3,3) e(1,2) e(2,3) e(1,3)]';
         end
         
         if (eEff > 0 )
            D = 2/3*(sEff/eEff)*ob.I4_dev + 4/9*(dsde-(sEff/eEff))/(eEff^2)*(e*e') + K/3*ob.I4_bulk;
            S = 2/3*sEff * ob.e/eEff;
         else
            D = 2/3*dsde*ob.I4_dev + (K/3)*ob.I4_bulk;
            S = zeros(ob.ndm);
         end
         
         sigma   = S + (1/3)*K*sum(gp.eps(1:ob.ndm))*ob.I;
         sigma_v = T2T1(sigma,1);
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, ~)
         if (ob.ndm == 2); gp.D = gp.D([1,2,4],[1,2,4]); end
         
         Kel = el.K + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, ~)
         Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
      end
   end
   
   methods (Static)
      function se = effectivestress(ee,s0,e0,n)
         if (ee < e0)
            if (n-1<10^(-12))
               se = s0*ee/e0;
            else
               se = s0*(sqrt( (1+n^2)/(n-1)^2 - (n/(n-1)-ee/e0)^2 )-1/(n-1));
            end
         else
            se = s0*( (ee/e0)^(1/n)  );
         end
         
      end
      function dsde = hardeningslope(ee,s0,e0,n)
         if (ee < e0)
            if (n-1<10^(-12))
               dsde = s0/e0;
            else
               dsde = s0*(n/(n-1)-ee/e0)/e0;
               dsde = dsde/sqrt( (1+n^2)/(n-1)^2 - (n/(n-1)-ee/e0)^2 );
            end
         else
            dsde = s0*( (ee/e0)^(1/n)  )/(n*ee);
         end
         
      end
   end
end