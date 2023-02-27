classdef PlasticRI
  % Rate Independent Plastic class
     
  properties (Hidden, SetAccess = private)
    ndm;
    ndof;
    dNdX;
    ngp;
    finiteDisp = 0;
      
    I;
    I4_dev;
    I4_bulk;
    linear = false;
  end
   
  properties (SetAccess = private)
    E;
    v;
    Kp;
    Hp;
    Y;
    G;
    K;
      
    eEff;
    e;
    ep;
    theta;
    bkStrss;
  end

  methods
    %% Construct
    function ob = PlasticRI(num, props, identity)
      ob.ndm  = num.ndm;
      ob.ndof = num.ndof;
      ob.ngp  = num.gp;
      ob.eEff    = zeros(                 num.gp, num.el, num.steps+1);
      ob.bkStrss = zeros(num.ndm,num.ndm, num.gp, num.el, num.steps+1);
      ob.ep      = zeros(num.ndm,num.ndm, num.gp, num.el, num.steps+1);
         
      [ob.E, ob.v, ob.Kp, ob.Hp, ob.Y] = ob.getProps(props);
         
      ob.G =  0.5*ob.E/(1+  ob.v);
      ob.K =(1/3)*ob.E/(1-2*ob.v);
         
      ob.I      = eye(ob.ndm);
      ob.I4_dev = identity.I4_dev;
      ob.I4_bulk= identity.I4_bulk;
    end
      
    %% Compute gp K, F and in case of dynamics: M
    function [gp, el, ob] = computeKFM(ob, gp, el, step)
      gp.eps  = gp.B * el.Uvc;
      ob.e = T1T2(gp.eps,2) - (1/3)*trace(T1T2(gp.eps,2))*ob.I;
         
      % Values from previous step
      n.ep      = ob.ep     (:,:, gp.i, gp.iel, step);
      n.eEff    = ob.eEff   (     gp.i, gp.iel, step);
      n.bkStrss = ob.bkStrss(:,:, gp.i, gp.iel, step);
         
      % Trial equivalent stress
      np1.e    = ob.e;
      np1.S_ht = 2*ob.G*(np1.e - n.ep) - n.bkStrss;
      np1.sEqv = sqrt( (3/2)*sum(dot(np1.S_ht, np1.S_ht)) );
         
         
      np1.deEff = ob.getEffPlasStrnInc( ...
        n.eEff, ...
        ob.G, ...
        ob.Y, ...
        ob.Kp, ...
        ob.Hp, ...
        np1.sEqv);
      ob.eEff(gp.i, gp.iel, step+1) =  n.eEff + np1.deEff;
         
      if (np1.sEqv*np1.deEff>0)
        np1.dbkStrss = np1.deEff * ob.Hp * (np1.S_ht/np1.sEqv);
        np1.dep      = 3/2*np1.deEff* (np1.S_ht/np1.sEqv);
            
        np1.theta = 1 - 3*ob.G* (np1.deEff/np1.sEqv);
        c2        = (1-np1.theta) - 1/(1+ (ob.Kp+ob.Hp)/(3*ob.G)) ;
        c3        = 3/2*c2/(np1.sEqv)^2;
      else
        np1.dbkStrss = zeros(size(np1.S_ht));
        np1.dep      = 0;
            
        np1.theta    = 1;
        c3 = 0;
      end
      np1.ep = n.ep + np1.dep;
         
      % Save states
      ob.bkStrss(:,:, gp.i, gp.iel, step+1) = n.bkStrss + np1.dbkStrss;
      ob.ep     (:,:, gp.i, gp.iel, step+1) = np1.ep;
         
      % Compute D
      if (ob.ndm == 2)
        S=[...
          np1.S_ht(1,1) np1.S_ht(2,2) 0 ...
          np1.S_ht(1,2) 0             0];
      
      elseif (ob.ndm == 3)
        S=[...
          np1.S_ht(1,1) np1.S_ht(2,2) np1.S_ht(3,3)...
          np1.S_ht(1,2) np1.S_ht(2,3) np1.S_ht(1,3)];
      end
         
      gp.D =...
        ob.K*ob.I4_bulk              +... % Elastic bulk
        np1.theta*(2*ob.G)*ob.I4_dev +... % Elastic deviatoric
        c3*(2*ob.G)*(S'*S);               % Plastic
         
      sigma    = 2*ob.G*(np1.e - np1.ep) + ob.K*sum(gp.eps(1:ob.ndm))*ob.I;
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
    function deEff = getEffPlasStrnInc(eEff_n, G, Y, Kp,Hp, sEquiv)
      deEff = (sEquiv - Y - eEff_n*Kp) / (Kp + Hp + 3*G);
      if (deEff < 0)
        deEff = 0;
      end
    end

    function [E, v, K, H, Y] = getProps(props)
      for i=1:length(props)
        switch props{i,1}
          case 'E'
            E = props{i,2};
          case 'v'
            v = props{i,2};
          case 'K'
            K = props{i,2};
          case 'H'
            H = props{i,2};
          case 'Y'
            Y = props{i,2};
          otherwise
            error([...
              'You''ve chosen Rate Independent Plastic material ' ...
              'but specified incompatible material properties, ' ...
              'I''m disapponted']);
        end
      end
    end
  end
end

