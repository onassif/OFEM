classdef MTS
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      nen;
      ngp;
      dNdX;
      finiteDisp = 1;
      
      I
      I2 = eye(3);
      I4_dev
      I4_bulk
      linear = false;
      
      hard
      
      name = 'MTS';
      
      schmidt
      nslip
      ms
      qs
      q_cr
      
      permut = [...
         0  0  0  0  0  1  0 -1  0
         0  0 -1  0  0  0  1  0  0
         0  1  0 -1  0  0  0  0  0]';
   end
   
   properties (SetAccess = private)
      props     = struct('E',[], 'nu',[], 'beta',[], 'Hp',[], 'Y',[], 'G',[], 'Bulk',[]);
      hardProps = struct(...
         'tau_ht_y' ,[], 'g0_y'   ,[], 'q_y'        ,[], 'p_y'        ,[], 'eps0_dt_y',[],...
         'tau_ht_v' ,[], 'g0_v'   ,[], 'q_v'        ,[], 'p_v'        ,[], 'eps0_dt_v',[],...
         'tau_a'    ,[], 'k0'     ,[], 'boltz'      ,[], 'b'          ,[], 'theta0'   ,[],...
         'mu0'      ,[], 'D0'     ,[], 't0'         ,[], 'hardenExp'  ,[], 'E'        ,[],...
         'nu'       ,[], 'voche_m',[], 'elasticType',[], 'numCrystals',[], 'C0'       ,[],...
         'angles'   ,struct('conv',[], 'type'       ,[], 'val'        ,[]));
      
      DList
      qList
      deList
      de
      ep
      S
      strss
      eEff
      sEquiv
      tauTilde
      gradFeinv
      Rp
      RpList
      R
      RList
      Uvc_n
      
      plastic = false;
      dir;
      
      dt;
      
      slip;
      gRot;
      
      mu_harden;
   end
   properties
      atol1 = 1.0e-6;
      rtol1 = 1.0e-12;
      atol2 = 1e-7;
      rtol2 = 1.0e-12;
      atol  = 1e-9;
      miter = 50;
      mmin  = 1;
      c     = 1.0e-4;
      red   = 0.5;
      mls   = 4;
      Jtol1 = 1e-4;
      Jtol2 = 1e-4;
      mcuts = 8;
   end
   %%
   methods
      %% Construct
      function ob = MTS(num, props, hardProps, slip, time, identity)
         ob.ndm   = num.ndm;
         ob.ndof  = num.ndof;
         ob.nen   = num.nen;
         ob.ngp   = num.gp;
         ob.eEff  = zeros(                num.gp, num.el, num.steps+1);
         ob.ep    = zeros(ob.ndm, ob.ndm, num.gp, num.el, num.steps+1);
         ob.DList = zeros(num.str, num.str, num.gp, num.el, num.steps+1);
         ob.qList = zeros(num.str, num.str, num.gp, num.el, num.steps+1);
         ob.deList= zeros(num.str, num.gp, num.el, num.steps+1); 
         ob.strss = zeros(num.str, num.gp, num.el, num.steps+1);
         ob.tauTilde  = zeros(num.gp, num.el, num.steps+1);
         ob.gradFeinv = zeros(3, 9, num.gp, num.el, num.steps+1);
         ob.RpList    = zeros(3, 3, num.gp, num.el, num.steps+1);
         ob.RList     = zeros(3, 3, num.gp, num.el, num.steps+1);
         ob.RpList(1,1,:,:,1) = 1;     ob.RpList(2,2,:,:,1) = 1;     ob.RpList(3,3,:,:,1) = 1;
         ob.RList(1,1,:,:,1)  = 1;     ob.RList(2,2,:,:,1)  = 1;     ob.RList(3,3,:,:,1)  = 1;
         
         ob.Uvc_n = zeros(num.nen*num.ndof, num.el, num.steps+1);
         
         for i=1:length(props)
            switch props{i,1}
               case 'E'
                  ob.props.E     = props{i,2};
               case 'nu'
                  ob.props.nu    = props{i,2};
               case 'beta'
                  ob.props.beta  = props{i,2};
               case 'Hp'
                  ob.props.Hp    = props{i,2};
               case 'Y'
                  ob.props.Y     = props{i,2};
               otherwise
                  error( ['You''ve chosen mm10 material but specified incompatible ',...
                     'material properties, I''m disapponted']);
            end
         end
         for i=1:length(hardProps)
            switch hardProps{i,1}
               case 'tau_ht_y'
                  ob.hardProps.tau_ht_y    = hardProps{i,2};
               case 'g0_y'
                  ob.hardProps.g0_y        = hardProps{i,2};
               case 'q_y'
                  ob.hardProps.q_y         = hardProps{i,2};
               case 'p_y'
                  ob.hardProps.p_y         = hardProps{i,2};
               case 'eps0_dt_y'
                  ob.hardProps.eps0_dt_y   = hardProps{i,2};
               case 'tau_ht_v'
                  ob.hardProps.tau_ht_v    = hardProps{i,2};
               case 'g0_v'
                  ob.hardProps.g0_v        = hardProps{i,2};
               case 'q_v'
                  ob.hardProps.q_v         = hardProps{i,2};
               case 'p_v'
                  ob.hardProps.p_v         = hardProps{i,2};
               case 'eps0_dt_v'
                  ob.hardProps.eps0_dt_v   = hardProps{i,2};
               case 'tau_a'
                  ob.hardProps.tau_a       = hardProps{i,2};
               case 'k0'
                  ob.hardProps.k0          = hardProps{i,2};
               case 'boltz'
                  ob.hardProps.boltz       = hardProps{i,2};
               case 'b'
                  ob.hardProps.b           = hardProps{i,2};
               case 'theta0'
                  ob.hardProps.theta0      = hardProps{i,2};
               case 'mu0'
                  ob.hardProps.mu0         = hardProps{i,2};
               case 'D0'
                  ob.hardProps.D0          = hardProps{i,2};
               case 't0'
                  ob.hardProps.t0          = hardProps{i,2};
               case 'hardenExp'
                  ob.hardProps.hardenExp   = hardProps{i,2};
               case 'E'
                  ob.hardProps.E           = hardProps{i,2};
               case 'nu'
                  ob.hardProps.nu          = hardProps{i,2};
               case 'voche_m'
                  ob.hardProps.voche_m     = hardProps{i,2};
               case 'elasticType'
                  ob.hardProps.elasticType = hardProps{i,2};
               case 'numCrystals'
                  ob.hardProps.numCrystals = hardProps{i,2};
               case 'angleConv'
                  ob.hardProps.angles.conv = hardProps{i,2};
               case 'angleType'
                  ob.hardProps.angles.type = hardProps{i,2};
               case 'angles'
                  ob.hardProps.angles.val  = hardProps{i,2};
               otherwise
                  error("Unknown MTS hardening parameter");
            end
         end
         
         ob.nslip = slip.nslip;
         ob.slip.n  = slip.n;
         [ob.gRot, ob.schmidt]= ob.compute_gRot_schmidt(ob.hardProps.angles.val, slip.b, slip.n, slip.nslip);
         ob.props.G    =   0.5*ob.hardProps.E/(1+  ob.hardProps.nu);
         ob.props.Bulk = (1/3)*ob.hardProps.E/(1-2*ob.hardProps.nu);
         
         ob.I       = eye(num.ndm);
         ob.I4_dev  = identity.I4_dev;
         ob.I4_bulk = identity.I4_bulk;
         
         ob.hardProps.C0 = ob.props.Bulk*ob.I4_bulk + 2*ob.props.G*ob.I4_dev;
         
         ob.dt = zeros(sum(time(:,2),1), 1);
         for i=1:size(time,1)
            for j = 1:time(i,2)
               if i == 1
                  ob.dt(j) = time(i,1)/time(i,2);
               else
                  ob.dt(  j+sum( time(1:i-1,2) )  ) = (time(i,1)-time(i-1,1))/time(i,2);
               end
            end
         end
         
      end
      %% Epsilon
      function [eps, ob] = computeStrain(ob, gp, el, step)
         eps   = gp.q_hf' * gp.B_hf(1:6,:) *  el.Uvc;
         if el.iter ==0 && step>1
            ob.de = gp.q_hf' * gp.B_hf(1:6,:) * (el.Uvc - ob.Uvc_n(:,el.i,step-1));
         else
            ob.de = gp.q_hf' * gp.B_hf(1:6,:) * (el.Uvc - el.Uvc_n);
         end
         ob.qList(:, :, gp.i, el.i, step+1) = gp.q;
         ob.deList(  :, gp.i, el.i, step+1) = eps;
         ob.Uvc_n(:,el.i, step+1) = el.Uvc ;
      end
      %% Sigma
      function [sigma_voigt, ob] = computeCauchy(ob, gp, el, step)
         if el.iter == 0 && step > 1
            sigma_voigt = ob.qList(:, :, gp.i, el.i, step)*(ob.strss(:,gp.i,el.i, step) + ob.DList( :,:,gp.i,el.i, step)*ob.de);
            ob.S = sigma_voigt;
         else
            sigma_voigt = ob.S;
         end
      end
      %% Tangential stiffness
      function [D, ctan, ob] = computeTangentStiffness(ob, gp, el, step)
         if el.iter == 0 && step > 1
            D = ob.DList( :,:,gp.i,el.i, step);
            
            gradFeinv   = ob.calcGradFeInv(ob.RList(:,:,:,el.i,step),ob.RpList(:,:,:,el.i,step-1),...
               gp.dXdxi_list(:,:,el.i), gp.xi');
            ob.gradFeinv( :,:,:,el.i, step) = reshape(gradFeinv,3,9,8);
         else
            dt   = ob.dt(step);
            % Identities
            ob.Rp = ob.RpList(:, :, gp.i, el.i, step);
            ob.R  = gp.R;
            % Material properties
            mu_harden = ob.hardProps.mu0;
            tau_v     = ob.hardProps.tau_ht_v;
            tau_y     = ob.hardProps.tau_ht_y;
            tau_a     = ob.hardProps.tau_a;
            k0        = ob.hardProps.k0;
            theta0    = ob.hardProps.theta0;
            b         = ob.hardProps.b;
            nExp      = ob.hardProps.hardenExp;
            C0        = ob.hardProps.C0;
            voche_m   = ob.hardProps.voche_m;
            % Gauss-Point properties
            ms        = ob.ms;
            qs        = ob.qs;
            q_cr      = ob.q_cr;
            %          ob.calc_grads(8,1,gp.Rn_list(:,:,:,gp.iel),
            %          tau_l     =
            % Values from previous step
            
            np1.deEff = sqrt( 2/3* (ob.de(1:3)'*ob.de(1:3) + 0.5*ob.de(4:6)'*ob.de(4:6)) );
            n.tauTilde  = ob.tauTilde(gp.i, gp.iel, step);
            n.tauTilde  = (n.tauTilde>0)*n.tauTilde +(n.tauTilde<=0)*(np1.deEff~=0)*(tau_a+tau_y+0.1);
            n.S         = ob.strss(:,gp.i, gp.iel, step);
            n.sol       = [n.S;n.tauTilde];
            R_n         = ob.RList( :,:,gp.i,el.i, step);
            gradFeinv   = ob.gradFeinv(:,:,gp.i, gp.iel, step);
            n_cr        = R_n*ob.Rp'*ob.gRot'*ob.slip.n';
            
            tm = gradFeinv*ob.permut*n_cr;
            tau_l = k0*b*mu_harden^2*sqrt(diag(tm'*tm))/(18*theta0); %% Change
            
            
            frac = 0;
            stp = 1;
            S = n.S;
            tauT = n.tauTilde;
            ostress = S;
            ott = tauT;
            cuts = 0;
            fail = 0;
            while (frac<1)
               iter = 0;
               R  = reshape(gp.R, 9,1);
               de = ob.de*(stp+frac);
               deEff = np1.deEff*(stp+frac);
               dt = 1; %%%%%%%%%%% fix later
               rss   = ms*S;
               if tauT ==0
                  gamm = zeros(size(rss));
               else
                  gamm  = deEff * (rss./tauT).^nExp .* sign(rss);
               end
               dbarp = ms'*gamm;
               Wp    = q_cr'*gamm;
               np1.R1 = S - n.S - C0*(de - dbarp) + ob.twoSymm(S,Wp);
               ct = (tauT-tau_a) - tau_y;
               xi = 1 - ct/tau_v + tau_l/ct;
               h = n.tauTilde + theta0*sum(abs(xi.^voche_m .* gamm));
               np1.R2 = tauT - h;
               np1.R = [np1.R1;np1.R2];
               nrm.R1 = norm(np1.R1); nrm.iR1 = norm(np1.R1);
               nrm.R2 = norm(np1.R2); nrm.iR2 = norm(np1.R2);
               while (nrm.R1>ob.atol1) && (nrm.R1/nrm.iR1>ob.rtol1) ||...
                     (nrm.R2>ob.atol2) && (nrm.R2/nrm.iR2>ob.rtol2)
                  dgammdrss  = (nExp*deEff)*abs(rss).^(nExp-1)/(tauT^nExp);
                  dgammdtauT = (-nExp/tauT)*gamm;
                  
                  corr = ob.twoSymm(S,q_cr');
                  J11 = (C0*ms'+corr)*(dgammdrss.*ms)+ob.IW(Wp)+eye(6);
                  J12 = (C0*ms'+corr)*dgammdtauT;
                  J21 = -theta0*(xi.^voche_m .* sign(rss))'*(dgammdrss.*ms);
                  J22 = 1+theta0*sum(abs((voche_m./xi.*(1/tau_v+tau_l/ct^2)+nExp/tauT).*xi.^(voche_m).*gamm));
                  if (imag(J22)~=0)
                     J22 = NaN;
                  end
                  J = [J11 J12;J21 J22];
                  dx = J\np1.R;
                  
                  alpha = 1;
                  ls  = 0;
                  x = [S;tauT];
                  while ls < ob.mls
                     ls = ls + 1;
                     nlsx = (0.5 - ob.c*alpha) *norm(np1.R)^2;
                     xnew = x - alpha*dx;
                     S    = xnew(1:6);
                     tauT = xnew(7:end);
                     
                     rss   = ms*S;
                     gamm  = deEff * (rss./tauT).^nExp .* sign(rss);
                     dbarp = ms'*gamm;
                     Wp    = q_cr'*gamm;
                     np1.R1 = S - n.S - C0*(de - dbarp) + ob.twoSymm(S,Wp);
                     ct = (tauT-tau_a) - tau_y;
                     xi = 1 - ct/tau_v + tau_l/ct;
                     h = n.tauTilde + theta0*sum(abs(xi.^voche_m .* gamm));
                     np1.R2 = tauT - h;
                     np1.R = [np1.R1;np1.R2];
                     nrm.R1 = norm(np1.R1);
                     nrm.R2 = norm(np1.R2);
                     nRs = 0.5*norm(np1.R)^2;
                     if ((nRs <= nlsx) || (ls >= ob.mls)) % || (cuts > 6 && iter <= 3)
                        break
                     else
                        alpha = ob.red*alpha;
                     end
                  end
                  
                  iter = iter + 1;
                  % c           Increment and check for failure
                  if ((iter > ob.miter) || any(isnan(x)) || any(imag(x)))
                     if  (cuts > 5 && nR1 < 5 && nR2 < 5e4)
                        S    = n.S;
                        tauT = n.tauTilde;
                        break
                     else
                        fail = 1;
                        break;
                     end
                  end
               end
               if (fail)
                  ob.S = ostress;
                  tauT = ott;
                  stp = stp * 0.5;
                  cuts = cuts + 1;
                  if (cuts > ob.mcuts), break
                  end
                  fail = 0;
               else
                  ostress = S;
                  ott = tauT;
                  frac = frac + stp;
               end
               ob.S = S;
            end
            
            de_mod  = [de(1:3); 0.5*de(4:6)];
            if deEff == 0
               dgammde = zeros(size(gamm,1),size(de_mod,1));
               J11 = eye(size(de_mod,1));
               J12 = zeros(size(de_mod,1),1);
               J21 = zeros(1,size(de_mod,1));
               J22 = eye(size(tauT,1));
            else
               dgammde = (2/3)*(1/deEff^2)*gamm*de_mod';
            end
            dSde = theta0*(abs(dgammde')*xi).*sign(de_mod);
            corr = ob.twoSymm(S,q_cr');
            
            JA = (C0*ms'+corr)*dgammde;
            JB = J12/J22*dSde';
            
            JJ = J11 - J12/J22*J21;
            JR = C0 - JA - JB;
            
            D = JJ\JR;
            D = 0.5*(D+D');
            
            wp = qs'*gamm;
            wbarp = [...
               0      wp(3) wp(2)
               -wp(3)  0     wp(1)
               -wp(2) -wp(1) 0];
            
            ob.DList( :,:,gp.i,el.i, step+1) = D;
            ob.RpList(:,:,gp.i,el.i, step+1) = expm(wbarp)*ob.RpList(:,:,gp.i,el.i, step);
            ob.RList( :,:,gp.i,el.i, step+1) = gp.R;
            ob.tauTilde(  gp.i,el.i, step+1) = tauT;
            ob.strss(   :,gp.i,el.i, step+1) = ob.S;
         end
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         if ob.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, step)
         % Definitions
         if el.iter == 0 && step > 1
            S = ob.strss(   :, gp.i, el.i, step);
            q = ob.qList(:, :, gp.i, el.i, step);
            S = q*S;
            gp.U = (gp.U + gp.dU)';
         else
            S = gp.sigma_un;
            q = gp.q;
         end
         B=gp.B;
         D = [q*gp.D*q' zeros(6,3); zeros(3,9)] - [...
            +S(1)      0         0         0.5*S(4)          0                 0.50*S(6)        -0.5*S(4)          0                 0.50*S(6)
            +0         S(2)      0         0.5*S(4)          0.50*S(5)         0                 0.5*S(4)         -0.5*S(5)          0
            +0         0         S(3)      0                 0.50*S(5)         0.50*S(6)         0                 0.5*S(5)         -0.50*S(6)
            +0.5*S(4)  0.5*S(4)  0         0.25*(S(2)+S(1))  0.25*S(6)         0.25*S(5)        -0.25*(S(2)-S(1)) -0.25*S(6)         0.25*S(5)
            +0         0.5*S(5)  0.5*S(5)  0.25*S(6)         0.25*(S(3)+S(2))  0.25*S(4)         0.25*S(6)        -0.25*(S(3)-S(2)) -0.25*S(4)
            +0.5*S(6)  0         0.5*S(6)  0.25*S(5)         0.25*S(4)         0.25*(S(1)+S(3)) -0.25*S(5)         0.25*S(4)        -0.25*(S(1)-S(3))
            -0.5*S(4)  0.5*S(4)  0        -0.25*(S(2)-S(1))  0.25*S(6),       -0.25*S(5)        -0.25*(S(1)+S(2))  0.25*S(6)         0.25*S(5)
            +0        -0.5*S(5)  0.5*S(5) -0.25*S(6)        -0.25*(S(3)-S(2))  0.25*S(4)         0.25*S(6)        -0.25*(S(2)+S(3))  0.25*S(4)
            +0.5*S(6)  0        -0.5*S(6)  0.25*S(5),       -0.25*S(4)        -0.25*(S(1)-S(3))  0.25*S(5)         0.25*S(4)        -0.25*(S(1)+S(3))];
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if gp.i == 1
            Kel = (B'*D*B)*gp.j*gp.w;
         else
            Kel = el.K + (B'*D*B)*gp.j*gp.w;
         end
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, step)
         if el.iter == 0 && step > 1
            sigma = [gp.sigma;0;0;0];
            gp.U = (gp.U + gp.dU)';
         else
            sigma = [gp.sigma_un;0;0;0];
         end
         if gp.i == 1
            Fint = (gp.B'*sigma) *gp.j *gp.w;
         else
            Fint = el.Fint + (gp.B'*sigma) *gp.j *gp.w;
         end
      end
      %% m(slip)
      function value = get.ms(ob)
         value = zeros(ob.nslip,6);
         RE = ob.RT2RVE(ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)+ob.schmidt(:,:,s)');
            value(s,:) = [A(1,1) A(2,2) A(3,3) 2*A(1,2) 2*A(2,3) 2*A(1,3)]*RE';
         end
      end
      %% q(slip)
      function value = get.qs(ob)
         value = zeros(ob.nslip,3);
         RW  = ob.RT2RVW(ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)-ob.schmidt(:,:,s)');
            value(s,:) = [A(2,3) A(1,3) A(1,2)]*RW';
         end
      end
      %% q(slip) in current coordinates
      function value = get.q_cr(ob)
         value = zeros(ob.nslip,3);
         RWC  = ob.RT2RVW(ob.R*ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)-ob.schmidt(:,:,s)');
            value(s,:) = [A(2,3) A(1,3) A(1,2)]*RWC';
         end
      end
      
   end
   methods (Static)
      function [gRot, schmidt] = compute_gRot_schmidt(angles, bi, ni, nslip)
         r = angles(1)*pi()/180;
         s = angles(2)*pi()/180;
         t = angles(3)*pi()/180;
         
         gRot = [...
            -sin(r)*sin(t)-cos(r)*cos(t)*cos(s), cos(r)*sin(t)-sin(r)*cos(t)*cos(s), cos(t)*sin(s)
            +sin(r)*cos(t)-cos(r)*sin(t)*cos(s),-cos(r)*cos(t)-sin(r)*sin(t)*cos(s), sin(t)*sin(s)
            +cos(r)*sin(s)                     , sin(r)*sin(s)                     , cos(s)       ];
         
         bs = bi * gRot;
         ns = ni * gRot;
         schmidt = zeros(3,3,nslip);
         for s = 1:nslip
            schmidt(:,:,s) = bs(s,:)'*ns(s,:);
         end
      end
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
         
         gradFeInv = repmat(reshape(grads(1:3,1:3,1:3), 27, 1),1,ngp);
      end
      function res = twoSymm(S,Wp)
         res = [...
            2*S(4)*Wp(3,:) - 2*S(6)*Wp(2,:)
            2*S(4)*Wp(3,:) - 2*S(5)*Wp(1,:)
            2*S(6)*Wp(2,:) + 2*S(5)*Wp(1,:)
            Wp(3,:)*(S(1)-S(2)) + Wp(1,:)*S(6) - Wp(2,:)*S(5)
            Wp(1,:)*(S(2)-S(3)) + Wp(2,:)*S(4) + Wp(3,:)*S(6)
            Wp(2,:)*(S(1)-S(3)) + Wp(1,:)*S(4) - Wp(3,:)*S(5)];
      end
      function res = IW(W)
         res = [...
            0     0     0     2*W(3)   0       -2*W(2)
            0     0     0     2*W(3)   -2*W(1)  0
            0     0     0     0         2*W(1)  2*W(2)
            W(3)  -W(3) 0     0        -1*W(2)  1*W(1)
            0     W(1)  -W(1) 1*W(2)     0      1*W(3)
            W(2)  0     -W(2) 1*W(1)   -1*W(3)  0.0000];
      end
      function RV  = RT2RVE(R)
         RV = [...
            R(1,1)^2      R(1,2)^2      R(1,3)^2      2*R(1,1)*R(1,2)             2*R(1,3)*R(1,2)             2*R(1,1)*R(1,3)
            R(2,1)^2      R(2,2)^2      R(2,3)^2      2*R(2,1)*R(2,2)             2*R(2,3)*R(2,2)             2*R(2,1)*R(2,3)
            R(3,1)^2      R(3,2)^2      R(3,3)^2      2*R(3,1)*R(3,2)             2*R(3,3)*R(3,2)             2*R(3,1)*R(3,3)
            R(1,1)*R(2,1) R(1,2)*R(2,2) R(1,3)*R(2,3) R(1,1)*R(2,2)+R(2,1)*R(1,2) R(1,2)*R(2,3)+R(1,3)*R(2,2) R(1,1)*R(2,3)+R(1,3)*R(2,1)
            R(2,1)*R(3,1) R(3,2)*R(2,2) R(2,3)*R(3,3) R(2,1)*R(3,2)+R(2,2)*R(3,1) R(2,2)*R(3,3)+R(3,2)*R(2,3) R(2,1)*R(3,3)+R(2,3)*R(3,1)
            R(1,1)*R(3,1) R(1,2)*R(3,2) R(1,3)*R(3,3) R(1,1)*R(3,2)+R(1,2)*R(3,1) R(1,2)*R(3,3)+R(1,3)*R(3,2) R(1,1)*R(3,3)+R(3,1)*R(1,3)];
      end
      function RV = RT2RVW(R)
         RV = [...
         R(2,2)*R(3,3)-R(2,3)*R(3,2) R(2,1)*R(3,3)-R(2,3)*R(3,1) R(2,1)*R(3,2)-R(2,2)*R(3,1)
         R(1,2)*R(3,3)-R(1,3)*R(3,2) R(1,1)*R(3,3)-R(1,3)*R(3,1) R(1,1)*R(3,2)-R(1,2)*R(3,1)
         R(1,2)*R(2,3)-R(1,3)*R(2,2) R(1,1)*R(2,3)-R(1,3)*R(2,1) R(1,1)*R(2,2)-R(1,2)*R(2,1)];
      end

   end
end
