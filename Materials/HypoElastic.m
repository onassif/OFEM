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
   end
   
   properties (SetAccess = private)
      sEff
      eEff
      e
   end
   %%
   methods
      %% Construct
      function obj = HypoElastic(num, props, identity)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.ngp  = num.gp;
         
         for i=1:length(props)
            switch props{i,1}
               case 'sig_0'
                  obj.sig0 = props{i,2};
               case 'eps_0'
                  obj.eps0 = props{i,2};
               case 'n'
                  obj.nExp = props{i,2};
               case 'nu'
                  obj.nu   = props{i,2};
               otherwise
                  error(['You''ve chosen Hypo Elastic material but specified',...
                     ' incompatible material properties, I''m disapponted']);
            end
         end
         
         obj.K = (obj.nExp*obj.sig0/obj.eps0) / (3*(1-2*obj.nu));
         
         obj.I      = eye(num.ndm);
         obj.I4_dev = identity.I4_dev;
         obj.I4_bulk= identity.I4_bulk;
         
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         I = obj.I;
         
         eps = gp.B * el.Uvc;
         
         if     obj.ndm == 2
            e=[1.0*eps(1) 0.5*eps(3)
               0.5*eps(3) 1.0*eps(2)];
         elseif obj.ndm == 3
            e=[1.0*eps(1) 0.5*eps(4) 0.5*eps(6)
               0.5*eps(4) 1.0*eps(2) 0.5*eps(5)
               0.5*eps(6) 0.5*eps(5) 1.0*eps(3)];
         end
         e = e - (1/3)*trace(e)*I;
         obj.e = e;
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~, ~)
         K    = obj.K;
         I    = obj.I;
         
         e    = obj.e;
         eEff = obj.eEff;
         sEff = obj.sEff;
         
         if eEff>0
            S = (2/3)*sEff * e/eEff;
         else
            S = zeros(obj.ndm);
         end
         sigma = S + (1/3)*K*sum(gp.eps(1:obj.ndm))*I;
         
         sigma_voigt = obj.voigtize(sigma,'col', obj.ndm);
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, ~, ~, ~)
         K    = obj.K;
         sig0 = obj.sig0;
         eps0 = obj.eps0;
         nExp = obj.nExp;
         
         e = obj.e;
         eEff = sqrt( (2/3)*sum(dot(e,e)) );
         sEff = obj.effectivestress(eEff, sig0, eps0, nExp);
         dsde = obj.hardeningslope (eEff, sig0, eps0, nExp);
         
         obj.eEff = eEff;
         obj.sEff = sEff;
         
         I4_dev  = obj.I4_dev;
         I4_bulk = obj.I4_bulk;
         
         if     obj.ndm == 2
            e = [e(1,1) e(2,2) 0      e(1,2) 0      0     ]';
         elseif obj.ndm == 3
            e = [e(1,1) e(2,2) e(3,3) e(1,2) e(2,3) e(1,3)]';
         end
         
         if (eEff > 0 )
            D = 2/3*(sEff/eEff)*I4_dev + 4/9*(dsde-(sEff/eEff))/(eEff^2)*(e*e') + K/3*I4_bulk;
         else
            D = 2/3*dsde*I4_dev + (K/3)*I4_bulk;
         end
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         
         if     obj.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(~, gp, el, ~)
         if gp.i == 1
            Kel = (gp.B'*gp.D*gp.B) *gp.J *gp.w;
         else
            Kel = el.K + (gp.B'*gp.D*gp.B) *gp.J *gp.w;
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
      function vec = voigtize(mat, orientation, ndm)
         switch orientation
            case {'Col', 'col', 'Column', 'column'}
               if     ndm == 2
                  vec = [mat(1,1) mat(2,2) mat(1,2)]';
               elseif ndm == 3
                  vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)]';
               end
            case {'Row', 'row'}
               if     ndm == 2
                  vec = [mat(1,1) mat(2,2) mat(1,2)];
               elseif ndm == 3
                  vec = [mat(1,1) mat(2,2) mat(3,3) mat(1,2) mat(2,3) mat(1,3)];
               end
            otherwise
               error('Unknown input');
         end
      end
   end
end

