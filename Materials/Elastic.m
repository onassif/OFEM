classdef Elastic
  %Elastic 3D class
  properties (SetAccess = private)
    ndm;
    ndof;
    ngp;
    finiteDisp = 0;
    C0;
    dyn;
    rho;
  end
  %%
  methods
    %% Construct
    function ob = Elastic(num, props, identity, dynamic)
      ob.ndm  = num.ndm;
      ob.ndof = num.ndof;
      ob.ngp  = num.gp;
         
      ob.dyn = dynamic; 
         
      [E, nu, ob.rho] = ob.getProps(props);
         
      G =   0.5*E/(1+  nu);
      K = (1/3)*E/(1-2*nu);
         
      ob.C0 = K*identity.I4_bulk + 2*G*identity.I4_dev;
    end
    
    %% Compute gp K, F and in case of dynamics: M
    function [gp, el, ob] = computeKFM(ob, gp, el, ~)
         
      % Strain
      gp.eps = gp.B * el.Uvc;
         
      % Tangential Stiffness
      gp.D = ob.C0;
      if (ob.ndm == 2)
        D = ob.C0([1,2,4],[1,2,4]); 
      elseif (ob.ndm == 3) 
        D = ob.C0;
      end
         
      % Stress
      gp.sigma = D*gp.eps;

      % K
      el.K = el.K + (gp.B'*D*gp.B) *gp.J *gp.w;

      % F
      el.Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
         
      % M
      if ob.dyn == 1
        el.M = el.M + ob.rho*gp.NN   *gp.J*gp.w; 
      elseif ob.dyn == 2
        el.M = el.M + ob.rho*gp.NN_2 *gp.J*gp.w;
      end
    end
  end
  methods (Static)
    % Loop over different material properties
    function [E, nu, rho] = getProps(props)
      for j = 1:length(props)
        switch props{j,1}
          case 'E'
            E = props{j,2};
          case 'nu'
            nu = props{j,2};
          case 'rho'
            rho = props{j,2};
        end
      end
      
      % TODO: change this workaround to something more meaningful
      if ~exist('rho', 'var')
        rho = 1;
      end
    end
  end
end