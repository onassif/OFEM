classdef HypoElastic
  % Simple Hypo elastic routine
  % Currently only if ndm=2, ndof=2
   
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
      
    I;
    I4_dev;
    I4_bulk;
      
    e;
  end

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
    
    %% Compute gp K, F and in case of dynamics: M
    function [gp, el, ob] = computeKFM(ob, gp, el, ~)
      gp.eps  = gp.B * el.Uvc;
      ob.e = T1T2(gp.eps,2) - (1/3)*trace(T1T2(gp.eps,2))*ob.I;

         
      eEff = sqrt( (2/3)*sum(dot(ob.e,ob.e)) );
      sEff = ob.effectivestress(eEff, ob.sig0, ob.eps0, ob.nExp);
      dsde = ob.hardeningslope (eEff, ob.sig0, ob.eps0, ob.nExp);
         
      if (ob.ndm == 2)
        e_vec = [
          ob.e(1,1) 
          ob.e(2,2) 
          0.0      
          ob.e(1,2) 
          0.0      
          0.0     
        ];
      elseif (ob.ndm == 3)
        e_vec = [
          ob.e(1,1) 
          ob.e(2,2) 
          ob.e(3,3) 
          ob.e(1,2) 
          ob.e(2,3) 
          ob.e(1,3)
        ];
      end
         
      if (eEff > 0 )
        gp.D =... 
          2/3*(sEff/eEff)*ob.I4_dev +... 
          4/9*(dsde-(sEff/eEff))/(eEff^2)*(e_vec*e_vec') +...
          ob.K/3*ob.I4_bulk;
        S = 2/3*sEff * ob.e/eEff;
      else
        gp.D = 2/3*dsde*ob.I4_dev + (ob.K/3)*ob.I4_bulk;
        S = zeros(ob.ndm);
      end
         
      sigma   = S + (1/3)*ob.K*sum(gp.eps(1:ob.ndm))*ob.I;
      gp.sigma = T2T1(sigma,1);

      if (ob.ndm == 2)
        D = gp.D([1,2,4],[1,2,4]);
      elseif (ob.ndm == 3)
        D = gp.D;
      end
         
      el.K    = el.K    + (gp.B'*D*gp.B)   *gp.J *gp.w;
      el.Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
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