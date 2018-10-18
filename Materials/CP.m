classdef CP
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (SetAccess = private)
      ndm;
      ndof;
      nen;
      ngp;
      nstr;
      dNdX;
      Q
      finiteDisp = 1;
      
      I
      linear = false;
      
      hard
      
      name = 'CP';
      
      schmidt
      nslip
      ms
      qs
      q_cr
      
      
      props     = struct('E',[], 'nu',[], 'G',[], 'Bulk',[]);
      angles    = struct('angleConv',[], 'angleType',[], 'val',[]);
      
      list = struct('D',[], 'Rp',[], 'R', [], 'S',[], 'tauT', []);
      de
      S
      
      
      gradFeinv
      Rp
      R
      Uvc_n
      C0
      
      plastic = false;
      dir;
      
      dt;
      
      slip;
      gRot;
      
      mu_harden;
      
      cpM
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
      function ob = CP(num, props, cpType, cpProps, angles, slip, time, identity)
         ob.ndm   = num.ndm;
         ob.nstr  = num.str;
         ob.ndof  = num.ndof;
         ob.nen   = num.nen;
         ob.ngp   = num.gp;
         ob.list.D = zeros(num.str, num.str, num.gp, num.el, num.steps+1);
         ob.list.S = zeros(num.str, num.gp, num.el, num.steps+1);
         ob.list.tauT = zeros(num.gp, num.el, num.steps+1);
         ob.gradFeinv = zeros(3, 9, num.el, num.steps+1);
         ob.list.Rp   = zeros(3, 3, num.gp, num.el, num.steps+1);
         ob.list.R    = zeros(3, 3, num.gp, num.el, num.steps+1);
         ob.list.Rp(1,1,:,:,1) = 1;     ob.list.Rp(2,2,:,:,1) = 1;     ob.list.Rp(3,3,:,:,1) = 1;
         ob.list.R(1,1,:,:,1)  = 1;     ob.list.R(2,2,:,:,1)  = 1;     ob.list.R(3,3,:,:,1)  = 1;
         
         ob.Uvc_n = zeros(num.nen*num.ndof, num.el, num.steps+1);
         
         for i=1:length(props)
            switch props{i,1}
               case 'E'
                  ob.props.E     = props{i,2};
               case 'nu'
                  ob.props.nu    = props{i,2};
               otherwise
                  error( ['You''ve chosen mm10 material but specified incompatible ',...
                     'material properties, I''m disapponted']);
            end
         end
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
         [ob.gRot, ob.schmidt] = ob.compute_gRot_schmidt(ob.angles.val, slip.b, slip.n, slip.nslip);
         ob.props.G    =   0.5*ob.props.E/(1+  ob.props.nu);
         ob.props.Bulk = (1/3)*ob.props.E/(1-2*ob.props.nu);
         
         switch cpType
            case 'MTS'
               ob.cpM = MTS(cpProps, ob.gRot, slip, num.el, num.steps);
         end
         
         ob.C0 = ob.props.Bulk*identity.I4_bulk + 2*ob.props.G*identity.I4_dev;
         
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
         gp.U  = (gp.U_n + (1/2)*gp.dU)'; % n + 1/2
         Q     = ob.Qmat(gp.R);
         if (el.iter==0&&step>1)
            Un = ob.Uvc_n(:,el.i,step-1);
         else
            Un = ob.Uvc_n(:,el.i,step);
         end
         ob.de = Q'*gp.B*(el.Uvc - Un);
         
         eps = Q'*gp.B*el.Uvc;
         ob.Uvc_n(:,el.i, step+1) = el.Uvc ;
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
            sigma_v = ob.Qmat(gp.R)*(ob.list.S(:,gp.i,el.i, step) + D*ob.de);
            ob.S = sigma_v;
         else
            dt   = ob.dt(step);
            % Identities
            ob.Rp = ob.list.Rp(:, :, gp.i, el.i, step);
            ob.R  = gp.R;
            % Material properties
            
            C0        = ob.C0;
            % Gauss-Point properties
            ms        = ob.ms;
            q_cr      = ob.q_cr;
            ob.cpM.ms   = ob.ms;
            ob.cpM.qs   = ob.qs;
            ob.cpM.q_cr = ob.q_cr;
            % Values from previous step
            
            np1.deEff = sqrt( 2/3* (ob.de(1:3)'*ob.de(1:3) + 0.5*ob.de(4:6)'*ob.de(4:6)) );
            
            n.tauTilde  = ob.cpM.tauT_n;
            n.S         = ob.list.S(:,gp.i, gp.iel, step);
            
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
               
               R1 = S - n.S - C0*(de - dbarp) + ob.twoSymm(S,Wp);
               R2 = tauT - n.tauTilde - dtauT;
               R = [R1;R2];
               nrm.iR1 = norm(R1); nrm.iR2 = norm(R2);
               while (norm(R1)>ob.atol1) && (norm(R1)/nrm.iR1>ob.rtol1) || (norm(R2)>ob.atol2) && (norm(R2)/nrm.iR2>ob.rtol2)
                  dgammdtau  = ob.cpM.compute_dgammdtau;
                  dgammdtauT = ob.cpM.compute_dgammdtauT;
                  
                  corr = ob.twoSymm(S,q_cr');
                  J11 = (C0*ms'+corr)*(dgammdtau.*ms) + ob.IW(Wp) + eye(ob.nstr);
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
                     R1 = S - n.S - C0*(de - dbarp) + ob.twoSymm(S,Wp);
                     
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
            
            if deEff == 0
               D = C0;
            else
               dgammde = ob.cpM.compute_dgammde;
               
               dSde    = ob.cpM.compute_dSde(dgammde);
               JA = ( C0*ms' + ob.twoSymm(S,q_cr') )*dgammde;
               JB = J12/J22*dSde';
               JJ = J11 - J12/J22*J21;
               JR = C0 - JA - JB;
               D = JJ\JR;
               D = 0.5*(D+D');
            end
            
            wbarp = ob.cpM.compute_wbarp;
            
            ob.list.D( :,:,gp.i,el.i, step+1) = D;
            ob.list.Rp(:,:,gp.i,el.i, step+1) = expm(wbarp)*ob.list.Rp(:,:,gp.i,el.i, step);
            ob.list.R( :,:,gp.i,el.i, step+1) = gp.R;
            ob.list.tauT(  gp.i,el.i, step+1) = tauT;
            ob.list.S(   :,gp.i,el.i, step+1) = ob.S;
            
            sigma_v = ob.Qmat(gp.R)*ob.S;
         end
         
         if gp.i == ob.ngp
            ob.cpM.list = ob.list;
            ob.cpM = ob.cpM.endGPComp();
         end
      end
      %% Element K
      function Kel = computeK_el(ob, gp, el, step)
         % Definitions
         Q = ob.Qmat(gp.R);
         if el.iter == 0 && step > 1
            S = Q*ob.list.S(:, gp.i, el.i, step);
            gp.U = (gp.U + gp.dU)';
         else
            S = gp.sigma;
         end
         B=gp.Bf;
         cep  = Q*gp.j*gp.w*gp.D*Q';
         qbar = gp.j*gp.w*[...
            2*S(1) 0      0      S(4)            0               S(6)
            0      2*S(2) 0      S(4)            S(5)            0
            0      0      2*S(3) 0               S(5)            S(6)
            S(4)   S(4)   0      1/2*(S(1)+S(2)) 1/2*S(6)        1/2*S(5)
            0      S(5)   S(5)   1/2*S(6)        1/2*(S(2)+S(3)) 1/2*S(4)
            S(6)   0      S(6)   1/2*S(5)        1/2*S(4)        1/2*(S(1)+S(3))];
         Kmat = [cep - qbar, zeros(6,3); zeros(3,9)];
         Kgeo = gp.j*gp.w*formGeo(S);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         Kel = el.K + (B'*(Kmat+Kgeo)*B);
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el, step)
         sigma = [gp.sigma;0;0;0]; % already rotated
         if el.iter == 0 && step > 1
            gp.U = (gp.U + gp.dU)';
         end
         
         Fint = el.Fint + (gp.Bf'*sigma) *gp.j *gp.w;
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
      function value = Qmat(R)
         R11 = R(1,1); R12 = R(1,2); R13 = R(1,3);
         R21 = R(2,1); R22 = R(2,2); R23 = R(2,3);
         R31 = R(3,1); R32 = R(3,2); R33 = R(3,3);
         value = [...
            R11^2   	R12^2   	R13^2   	2*R11*R12        	2*R12*R13        	2*R11*R13
            R21^2   	R22^2   	R23^2   	2*R21*R22        	2*R22*R23        	2*R21*R23
            R31^2    R32^2   	R33^2   	2*R31*R32        	2*R32*R33        	2*R31*R33
            R11*R21  R12*R22  R13*R23  R11*R22+R12*R21	R12*R23+R13*R22   R11*R23+R13*R21
            R21*R31  R32*R22  R23*R33  R21*R32+R22*R31	R22*R33+R23*R32   R21*R33+R23*R31
            R11*R31  R12*R32  R13*R33  R11*R32+R12*R31	R12*R33+R13*R32   R11*R33+R13*R31];
      end
   end
end