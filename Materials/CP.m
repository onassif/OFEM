classdef CP
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   properties
      de
      list = struct('D',[], 'Rp',[], 'R', [], 'S',[], 'tauT', []);
   end
   properties (SetAccess = private)
      ndm;
      ngp;
      nstr;
      Q
      finiteDisp = 1;
      
      schmidt
      nslip
      ms
      qs
      q_cr
      
      angles    = struct('angleConv',[], 'angleType',[], 'val',[]);

      S

      Rp
      R
      Uvc_n
      C0
      
      dt;
      
      slip;
      gRot;
      
      cpM
   end
   properties
      atol1 = 1e-6;	atol2 = 1e-7;
      rtol1 = 1e-12;	rtol2 = 1e-12;
      Jtol1 = 1e-4;  Jtol2 = 1e-4;
      miter = 50;
      mmin  = 1;
      c     = 1.0e-4;
      red   = 0.5;
      mls   = 4;
      mcuts = 8;
   end
   %%
   methods
      %% Construct
      function ob = CP(num, props, cpType, cpProps, angles, slip, time, identity)
         ob.ndm   = num.ndm;
         ob.nstr  = num.str;
         ob.ngp   = num.gp;
         ob.list.D = zeros(   6, 6, num.gp, num.el2, num.steps+1);
         ob.list.S = zeros(      6, num.gp, num.el2, num.steps+1);
         ob.list.tauT = zeros(      num.gp, num.el2, num.steps+1);
         ob.list.Rp   = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.list.R    = zeros(3, 3, num.gp, num.el2, num.steps+1);
         ob.list.Rp(1,1,:,:,1) = 1;     ob.list.Rp(2,2,:,:,1) = 1;     ob.list.Rp(3,3,:,:,1) = 1;
         ob.list.R( 1,1,:,:,1) = 1;     ob.list.R( 2,2,:,:,1) = 1;     ob.list.R( 3,3,:,:,1) = 1;
         
         ob.Uvc_n = zeros(num.nen*num.ndof, num.el, num.steps+1);
         
         [Bulk, G] = ob.getProps(props);
         ob.C0     = Bulk*identity.I4_bulk + 2*G*identity.I4_dev;
         
         for i=1:length(angles)
            switch angles{i,1}
               case 'angleConv'
                  ob.angles.conv = angles{i,2};
               case 'angleType'
                  ob.angles.type = angles{i,2};
               case 'angles'
                  ob.angles.val  = angles{i,2};
               otherwise
                  error( ['You''ve chosen mm10 material but specified incompatible ',...
                     'material properties, I''m disapponted']);
            end
         end
         ob.nslip = slip.nslip;
         ob.slip.n  = slip.n;
         [ob.gRot, ob.schmidt] = gRotSchmidt(ob.angles.val, slip.b, slip.n, slip.nslip);
         
         switch cpType
            case 'MTS'
               ob.cpM = MTS(cpProps, ob.gRot, slip, num.el2, num.steps);
         end
         
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
      function [eps, ob] = Strain(ob, gp, el, step)
         gp.U  = (gp.U_n + 1/2*gp.dU); % n + 1/2
         if ob.ndm == 2
            Q = Qmat([gp.R zeros(2,1);zeros(1,2) 1]);
            Q = Q([1,2,4],[1,2,4]);
            ob.de = Q'*gp.B*el.UresVc;
            ob.de = [ob.de(1) ob.de(2) 0 ob.de(3) 0 0]';
         elseif ob.ndm == 3
            Q     = Qmat(gp.R);
            ob.de = Q'*gp.B*el.UresVc;
         end
         
         eps = Q'*gp.B*el.Uvc;
         ob.Uvc_n(:,el.i, step+1) = el.Uvc;
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, el, step)
         ob.cpM.list = ob.list;
         ob.cpM.step = step;
         ob.cpM.iter = el.iter;
         ob.cpM.iel  = el.i;
         ob.cpM.gp = gp;
         
         if el.iter == 0 && step > 1
            D = ob.list.D( :,:,gp.i,el.i, step);
            Q = Qmat(check2D(gp.R));
            sigma_v = Q*(ob.list.S(:,gp.i,el.i, step) + D*ob.de);
            ob.S = sigma_v;
            qbar = computeqbar(Q*ob.list.S(:,gp.i,el.i, step));
         else
%             dt   = ob.dt(step);
            % Identities
            ob.Rp = ob.list.Rp(:, :, gp.i, el.i, step); %eye(3)
            ob.R  = check2D(gp.R);

            
            % Material properties
            C0   = ob.C0;
            % Gauss-Point properties
            ms   = ob.ms;
            q_cr = ob.q_cr;
            ob.cpM.ms   = ob.ms;
            ob.cpM.qs   = ob.qs;
            ob.cpM.q_cr = ob.q_cr;
            % Values from previous step
            if ob.ndm == 2
               np1.deEff = sqrt( 2/3* (ob.de(1:2)'*ob.de(1:2) + 0.5*ob.de(3)'*ob.de(3)) );              
            elseif ob.ndm == 3
               np1.deEff = sqrt( 2/3* (ob.de(1:3)'*ob.de(1:3) + 0.5*ob.de(4:6)'*ob.de(4:6)) );
            end
            
            n.tauTilde  = ob.cpM.tauT_n;
            n.S         = ob.list.S(:,gp.i, gp.iel, step); %zeros(6,1);
            
            frac = 0;
            stp = 1;
            S = n.S;
            tauT = n.tauTilde;
            ob.cpM.S    = n.S;
            ob.cpM.tauT = n.tauTilde;
            
            ostress = S;
            ott = tauT;
            cuts = 0;
            fail = 0;
            while (frac<1)
               iter = 0;
               de    = ob.de*(stp+frac);        ob.cpM.de    = ob.de;
               deEff = np1.deEff*(stp+frac);    ob.cpM.deEff = deEff;
               
               dbarp = ob.cpM.compute_dbarp;
               Wp    = ob.cpM.compute_Wp;
               dtauT = ob.cpM.compute_dtauT;
               
               R1 = S - n.S - C0*(de - dbarp) + twoSymm(S,Wp);
               R2 = tauT - n.tauTilde - dtauT;
               R = [R1;R2];
               nrm.iR1 = norm(R1); nrm.iR2 = norm(R2);
               while (norm(R1)>ob.atol1) && (norm(R1)/nrm.iR1>ob.rtol1) || (norm(R2)>ob.atol2) && (norm(R2)/nrm.iR2>ob.rtol2)
                  dgammdtau  = ob.cpM.compute_dgammdtau;
                  dgammdtauT = ob.cpM.compute_dgammdtauT;
                  
                  corr = twoSymm(S,q_cr');
                  J11 = (C0*ms'+corr)*(dgammdtau.*ms) + IW(Wp) + eye(6);
                  J12 = (C0*ms'+corr)*dgammdtauT;
                  J21 = ob.cpM.compute_dR2dtau;
                  J22 = ob.cpM.compute_dR2dtauT;
                  
                  J = [J11 J12;J21 J22];
                  dx = J\R;
                  
                  alpha = 1;
                  ls  = 0;
                  x = [S;tauT];
                  while ls < ob.mls
                     ls = ls + 1;
                     nlsx = (0.5 - ob.c*alpha) *norm(R)^2;
                     xnew = x - alpha*dx;
                     S    = xnew(1:6);
                     tauT = xnew(7:end);
                     ob.cpM.S    = xnew(1:6);
                     ob.cpM.tauT = xnew(7:end);
                     
                     dbarp = ob.cpM.compute_dbarp;
                     Wp    = ob.cpM.compute_Wp;
                     R1 = S - n.S - C0*(de - dbarp) + twoSymm(S,Wp);
                     
                     dtauT = ob.cpM.compute_dtauT;
                     R2 = tauT -n.tauTilde - dtauT;
                     R = [R1;R2];
                     nRs = 0.5*norm(R)^2;
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
            Q = Qmat(check2D(gp.R));

            sigma_v = Q*ob.S;
            qbar = computeqbar(sigma_v);
            if deEff == 0
               D = C0;
            else
               dgammde = ob.cpM.compute_dgammde;
               
               dSde    = ob.cpM.compute_dSde(dgammde);
               JA = ( C0*ms' + twoSymm(S,q_cr') )*dgammde;
               JB = J12/J22*dSde';
               JJ = J11 - J12/J22*J21;
               JR = C0 - JA - JB;
               D = JJ\JR;
               D = 0.5*(D+D');
            end
            
            wbarp = ob.cpM.compute_wbarp;
            
            ob.list.D( :,:,gp.i,el.i, step+1) = D; % [116770.243255197,57513.7019018134,57513.7019018134,0,0,0[29628.2706766917]]
            ob.list.Rp(:,:,gp.i,el.i, step+1) = expm(wbarp)*ob.list.Rp(:,:,gp.i,el.i, step); %eye(3)
            ob.list.R( :,:,gp.i,el.i, step+1) = check2D(gp.R); %eye(3)
            ob.list.tauT(  gp.i,el.i, step+1) = tauT; %155.1
            ob.list.S(   :,gp.i,el.i, step+1) = ob.S; %zeros(6,1)
         end
         D  = Q*D*Q' - qbar;
         
         if gp.i == ob.ngp
            ob.cpM.list = ob.list;
            ob.cpM      = ob.cpM.endGPComp();
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, step)
         % Definitions
         Q = Qmat(check2D(gp.R));
         if el.iter == 0 && step > 1
            S = Q*ob.list.S(:, gp.i, el.i, step);
            gp.U = (gp.U + gp.dU);
         else
            S = gp.sigma;
         end
         B=gp.Bf;

         D = formCombD(S, gp.D, gp.finiteDisp);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         if ob.ndm == 2
            D = D([1,2,4,7],[1,2,4,7]);
         end
         Kel = el.K + gp.j*gp.w* (B'*D*B);
      end
      %% Element Fint
      function Fint = computeFint(ob, gp, el, step)
         sigma = [gp.sigma;0;0;0]; % already rotated
         if el.iter == 0 && step > 1
            gp.U = (gp.U + gp.dU);
         end
         if ob.ndm == 2
            sigma = sigma([1,2,4,7]);
         end
         Fint = el.Fint + (gp.Bf'*sigma) *gp.j *gp.w;
      end
      %% m(slip)
      function value = get.ms(ob)
         value = zeros(ob.nslip,6);
         RE = RT2RVE(ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)+ob.schmidt(:,:,s)');
            value(s,:) = [A(1,1) A(2,2) A(3,3) 2*A(1,2) 2*A(2,3) 2*A(1,3)]*RE';
         end
      end
      %% q(slip)
      function value = get.qs(ob)
         value = zeros(ob.nslip,3);
         RW  = RT2RVW(ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)-ob.schmidt(:,:,s)');
            value(s,:) = [A(2,3) A(1,3) A(1,2)]*RW';
         end
      end
      %% q(slip) in current coordinates
      function value = get.q_cr(ob)
         value = zeros(ob.nslip,3);
         RWC  = RT2RVW(ob.R*ob.Rp');
         for s = 1:ob.nslip
            A = 0.5*(ob.schmidt(:,:,s)-ob.schmidt(:,:,s)');
            value(s,:) = [A(2,3) A(1,3) A(1,2)]*RWC';
         end
      end
   end
   
   methods (Static)
      function [Bulk, G] = getProps(props)
         for j = 1:length(props)
            switch props{j,1}
               case 'E'
                  E  = props{j,2};
               case 'nu'
                  nu = props{j,2};
            end
         end
         G    = 1/2*E/(1+  nu);
         Bulk = 1/3*E/(1-2*nu);
      end
   end
end