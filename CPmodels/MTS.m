classdef MTS
   properties
      permut = [...
         0  0  0  0  0  1  0 -1  0
         0  0 -1  0  0  0  1  0  0
         0  1  0 -1  0  0  0  0  0]';
      pr = struct(...
         'tau_ht_y' ,[], 'g0_y'   ,[], 'q_y'        ,[], 'p_y'        ,[], 'eps0_dt_y',[],...
         'tau_ht_v' ,[], 'g0_v'   ,[], 'q_v'        ,[], 'p_v'        ,[], 'eps0_dt_v',[],...
         'tau_a'    ,[], 'k0'     ,[], 'boltz'      ,[], 'b'          ,[], 'theta0'   ,[],...
         'mu0'      ,[], 'D0'     ,[], 't0'         ,[], 'hardenExp'  ,[], 'E'        ,[],...
         'nu'       ,[], 'voche_m',[], 'elasticType',[], 'numCrystals',[], 'C0'       ,[]);
      
      tauT_n = 0;
      list;
      
      step;
      iter;
      iel;
      gp;
      gradFeinv;
      nel;
      nsteps;
      
      slip = struct('n',[]);
      gRot
      
      tau_l
      
      ms;
      qs;
      q_cr;
      
      deEff
      S
      tauT
      rss
      gamm
      de
      de_mod
      xi
   end
   
   methods
      %% Construct
      function ob = MTS(cpProps,gRot,slip,numel,numsteps)
         ob.pr = ob.saveProps(cpProps);

         ob.gradFeinv = zeros(3,9,numel,numsteps+1);
         ob.slip.n    = slip.n;
         ob.gRot      = gRot;
      end
      %% dbarp
      function dbarp = compute_dbarp(ob)
         dbarp = ob.ms'*ob.gamm;
      end
      %% wbarp
      function wbarp = compute_wbarp(ob)
         wp = ob.qs'*ob.gamm;
         wbarp = [...
            0       wp(3) wp(2)
            -wp(3)  0     wp(1)
            -wp(2) -wp(1) 0];
      end
      %% Wp
      function Wp = compute_Wp(ob)
         Wp = ob.q_cr'*ob.gamm;
      end
      %% dgamma/dtau
      function dgammdtau = compute_dgammdtau(ob)
         n = ob.pr.nExp;
         dgammdtau = (n*ob.deEff)*abs(ob.rss).^(n-1)/(ob.tauT^n);
      end
      %% dgamma/dtau_tilde  
      function dgammdtauT = compute_dgammdtauT(ob)
         n = ob.pr.nExp;
         dgammdtauT = (-n/ob.tauT)*ob.gamm;
      end
      %% deltaTau_tilde 
      function dtauT = compute_dtauT(ob)
         voche_m= ob.pr.voche_m;
         dtauT =  ob.pr.theta0*sum(abs(ob.xi.^voche_m .* ob.gamm));
      end
      %% dR2/dtau
      function dR2dtau = compute_dR2dtau(ob)
         theta0 = ob.pr.theta0;
         voche_m= ob.pr.voche_m;

         n = ob.pr.nExp;
         dgammadtau = (n*ob.deEff)*abs(ob.rss).^(n-1)/(ob.tauT^n);
         dR2dtau = -theta0*(ob.xi.^voche_m .* sign(ob.rss))'*(dgammadtau.*ob.ms);
      end
      %% dR2/dtauT
      function dR2dtauT = compute_dR2dtauT(ob)
         theta0 = ob.pr.theta0;
         tau_v  = ob.pr.tau_ht_v;
         m= ob.pr.voche_m;
         ct    = (ob.tauT-ob.pr.tau_a) - ob.pr.tau_ht_y;
         n = ob.pr.nExp;
         dR2dtauT = 1+theta0*sum(abs((m./ob.xi.*(1/tau_v+ob.tau_l/ct^2)+n/ob.tauT).*ob.xi.^(m).*ob.gamm));
         if (imag(dR2dtauT)~=0)
            dR2dtauT = NaN;
         end
      end
      %% dgamma/de
      function dgammde = compute_dgammde(ob)
         dgammde = (2/3)*(1/ob.deEff^2)*ob.gamm*ob.de_mod';
      end
      %% dS/de
      function dSde = compute_dSde(ob, dgammde)
         dSde = ob.pr.theta0*(abs(dgammde')*ob.xi).*sign(ob.de_mod);
      end
      %% Computations at the element level
      function ob = endGPComp(ob)
         ob.gradFeinv(:,:, ob.iel, ob.step+1) = ob.calcGradFeInv(...
            ob.list.R(:,:,:,ob.iel,ob.step+1),ob.list.Rp(:,:,:,ob.iel,ob.step),ob.gp.dXdxi_list(:,:,ob.iel), ob.gp.xi');
      end
      %% Get functions
      function val = get.tauT_n(ob)
         n.tauT  = ob.list.tauT(ob.gp.i, ob.iel, ob.step);
         tau_y   = ob.pr.tau_ht_y;
         tau_a   = ob.pr.tau_a;
         val = (n.tauT>0)*n.tauT +(n.tauT<=0)*(tau_a+tau_y+0.1);
      end
      
      function val = get.tau_l(ob)
         R  = ob.list.R( :,:,ob.gp.i,ob.iel,ob.step);
         Rp = ob.list.Rp(:,:,ob.gp.i,ob.iel,ob.step);
         n_cr = R*Rp'*ob.gRot'*ob.slip.n';
         tm = ob.gradFeinv(:,:, ob.iel, ob.step)*ob.permut*n_cr;
         
         val = ob.pr.k0*ob.pr.b*ob.pr.mu0^2*sqrt(diag(tm'*tm))/(18*ob.pr.theta0);
      end
      
      function val = get.rss(ob)
         val = ob.ms*ob.S;
      end
      
      function val = get.gamm(ob)
         val  = ob.deEff * (ob.rss./ob.tauT).^ob.pr.nExp .* sign(ob.rss);
      end
      
      function val = get.de_mod(ob)
         val  = [ob.de(1:3); 0.5*ob.de(4:6)];
      end
      
      function val = get.xi(ob)
         ct     = (ob.tauT-ob.pr.tau_a) - ob.pr.tau_ht_y;
         val    = 1 - ct/ob.pr.tau_ht_v + ob.tau_l/ct;
      end
      
   end
   %% Static methods
   methods (Static)
      function gradFeInv = calcGradFeInv(R, Rp, dXdxi, intermat)
         ngp = size(R,3);
         Rt  = zeros(ngp,3,3);
         
         % Get R components and stick in the right place
         jacinv = inv(dXdxi);
         for i = 1:ngp
            Rt(i,:,:) = jacinv*Rp(:,:,i)*R(:,:,i)';
         end
         
         %       Vector:
         LHS = reshape(Rt,ngp,9)';
         RHS = LHS/intermat;
         
         grads = reshape(RHS,3,3,3);
         
         gradFeInv = reshape(grads(1:3,1:3,1:3), 3, 9);
      end
      
      function props = saveProps(cpProps)
         for i=1:length(cpProps)
            switch cpProps{i,1}
               case 'tau_ht_y'
                  props.tau_ht_y    = cpProps{i,2};
               case 'g0_y'
                  props.g0_y        = cpProps{i,2};
               case 'q_y'
                  props.q_y         = cpProps{i,2};
               case 'p_y'
                  props.p_y         = cpProps{i,2};
               case 'eps0_dt_y'
                  props.eps0_dt_y   = cpProps{i,2};
               case 'tau_ht_v'
                  props.tau_ht_v    = cpProps{i,2};
               case 'g0_v'
                  props.g0_v        = cpProps{i,2};
               case 'q_v'
                  props.q_v         = cpProps{i,2};
               case 'p_v'
                  props.p_v         = cpProps{i,2};
               case 'eps0_dt_v'
                  props.eps0_dt_v   = cpProps{i,2};
               case 'tau_a'
                  props.tau_a       = cpProps{i,2};
               case 'k0'
                  props.k0          = cpProps{i,2};
               case 'boltz'
                  props.boltz       = cpProps{i,2};
               case 'b'
                  props.b           = cpProps{i,2};
               case 'theta0'
                  props.theta0      = cpProps{i,2};
               case 'mu0'
                  props.mu0         = cpProps{i,2};
               case 'D0'
                  props.D0          = cpProps{i,2};
               case 't0'
                  props.t0          = cpProps{i,2};
               case 'hardenExp'
                  props.nExp        = cpProps{i,2};
               case 'E'
                  props.E           = cpProps{i,2};
               case 'nu'
                  props.nu          = cpProps{i,2};
               case 'voche_m'
                  props.voche_m     = cpProps{i,2};
               case 'elasticType'
                  props.elasticType = cpProps{i,2};
               case 'numCrystals'
                  props.numCrystals = cpProps{i,2};
               otherwise
                  error("Unknown MTS hardening parameter");
            end
         end
      end
   end
end

