classdef MTS
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      nen;
      dNdX;
      finiteDisp = 1;
      
      I
      I2 = eye(3);
      I4_dev
      I4_bulk
      linear = false;
      
      hard
      
      name = 'mm10';
   end
   
   properties (SetAccess = private)
      props     = struct('E',[], 'nu',[], 'beta',[], 'Hp',[], 'Y',[], 'G',[], 'Bulk',[]);
      hardProps = struct(...
         'tau_ht_y' ,[], 'g0_y'   ,[], 'q_y'        ,[], 'p_y'        ,[], 'eps0_dt_y',[],...
         'tau_ht_v' ,[], 'g0_v'   ,[], 'q_v'        ,[], 'p_v'        ,[], 'eps0_dt_v',[],...
         'tau_a'    ,[], 'k0'     ,[], 'boltz'      ,[], 'b'          ,[], 'theta0'   ,[],...
         'mu0'      ,[], 'D0'     ,[], 't0'         ,[], 'hardenExp'  ,[], 'E'        ,[],...
         'nu'       ,[], 'voche_m',[], 'elasticType',[], 'numCrystals',[],...
         'angles'   ,struct('conv',[], 'type'       ,[], 'val'        ,[]));
      
      e
      ep
      eEff
      sEquiv
      
      plastic = false;
      dir;
      
      dt;
      
      slip;
      gRot;
      
      mu_harden;
   end
   %%
   methods
      %% Construct
      function obj = MTS(num, props, hardProps, slip, time, identity)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.nen  = num.nen;
         obj.eEff = zeros(                 num.gp, num.el, num.steps+1);
         obj.ep   = zeros(obj.ndm,obj.ndm, num.gp, num.el, num.steps+1);
         
         for i=1:length(props)
            switch props{i,1}
               case 'E'
                  obj.props.E     = props{i,2};
               case 'nu'
                  obj.props.nu    = props{i,2};
               case 'beta'
                  obj.props.beta  = props{i,2};
               case 'Hp'
                  obj.props.Hp    = props{i,2};
               case 'Y'
                  obj.props.Y     = props{i,2};
               otherwise
                  error( ['You''ve chosen mm10 material but specified incompatible ',...
                     'material properties, I''m disapponted']);
            end
         end
         for i=1:length(hardProps)
            switch hardProps{i,1}
               case 'tau_ht_y'
                  obj.hardProps.tau_ht_y    = hardProps{i,2};
               case 'g0_y'
                  obj.hardProps.g0_y        = hardProps{i,2};
               case 'q_y'
                  obj.hardProps.q_y         = hardProps{i,2};
               case 'p_y'
                  obj.hardProps.p_y         = hardProps{i,2};
               case 'eps0_dt_y'
                  obj.hardProps.eps0_dt_y   = hardProps{i,2};
               case 'tau_ht_v'
                  obj.hardProps.tau_ht_v    = hardProps{i,2};
               case 'g0_v'
                  obj.hardProps.g0_v        = hardProps{i,2};
               case 'q_v'
                  obj.hardProps.q_v         = hardProps{i,2};
               case 'p_v'
                  obj.hardProps.p_v         = hardProps{i,2};
               case 'eps0_dt_v'
                  obj.hardProps.eps0_dt_v   = hardProps{i,2};
               case 'tau_a'
                  obj.hardProps.tau_a       = hardProps{i,2};
               case 'k0'
                  obj.hardProps.k0          = hardProps{i,2};
               case 'boltz'
                  obj.hardProps.boltz       = hardProps{i,2};
               case 'b'
                  obj.hardProps.b           = hardProps{i,2};
               case 'theta0'
                  obj.hardProps.theta0      = hardProps{i,2};
               case 'mu0'
                  obj.hardProps.mu0         = hardProps{i,2};
               case 'D0'
                  obj.hardProps.D0          = hardProps{i,2};
               case 't0'
                  obj.hardProps.t0          = hardProps{i,2};
               case 'hardenExp'
                  obj.hardProps.hardenExp   = hardProps{i,2};
               case 'E'
                  obj.hardProps.E           = hardProps{i,2};
               case 'nu'
                  obj.hardProps.nu          = hardProps{i,2};
               case 'voche_m'
                  obj.hardProps.voche_m     = hardProps{i,2};
               case 'elasticType'
                  obj.hardProps.elasticType = hardProps{i,2};
               case 'numCrystals'
                  obj.hardProps.numCrystals = hardProps{i,2};
               case 'angleConv'
                  obj.hardProps.angles.conv = hardProps{i,2};
               case 'angleType'
                  obj.hardProps.angles.type = hardProps{i,2};
               case 'angles'
                  obj.hardProps.angles.val  = hardProps{i,2};
               otherwise
                  error("Unknown MTS hardening parameter");
            end
         end
         obj.slip = slip;
         obj.props.G    =   0.5*obj.props.E/(1+  obj.props.nu);
         obj.props.Bulk = (1/3)*obj.props.E/(1-2*obj.props.nu);
         
         obj.I       = eye(num.ndm);
         obj.I4_dev  = identity.I4_dev;
         obj.I4_bulk = identity.I4_bulk;
         
         obj.dt = zeros(sum(time(:,2),1), 1);
         for i=1:size(time,1)
            for j = 1:time(i,2)
               if i == 1
                  obj.dt(j) = time(i,1)/time(i,2);
               else
                  obj.dt(  j+sum( time(1:i-1,2) )  ) = (time(i,1)-time(i-1,1))/time(i,2);
               end
            end
         end
         
      end
      %% Epsilon
      function [eps, obj] = computeStrain(obj, gp, el, ~)
         I = obj.I;
         
         eps= gp.q * gp.B(1:6,:) * el.Uvc;
         
         e=[1.0*eps(1) 0.5*eps(4) 0.5*eps(6)
            0.5*eps(4) 1.0*eps(2) 0.5*eps(5)
            0.5*eps(6) 0.5*eps(5) 1.0*eps(3)];
         e = e - (1/3)*trace(e)*I;
         obj.e = e;
         obj.eEff = sqrt( (2/3)*sum( dot(e,e) ));
      end
      %% Sigma
      function [sigma_voigt, obj] = computeCauchy(obj, gp, step)
         I = obj.I;
         G = obj.G;
         K = obj.Bulk;
         
         np1.e  = obj.e;
         np1.ep = obj.ep(:,:, gp.i, gp.iel, step+1);
         
         sigma       = 2*G* (np1.e - np1.ep) + K*sum(gp.eps(1:obj.ndm))*I;
         sigma_voigt = obj.voigtize(sigma,'col', obj.ndm);
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, step)
         dt   = obj.dt(step);
         % Identities
         I4_dev  = obj.I4_dev;
         I4_bulk = obj.I4_bulk;
         % Material properties
         mu_harden = obj.hardProps.mu0;
         tau_v     = obj.hardProps.tau_ht_v;
         tau_y     = obj.hardProps.tau_ht_y;
         % Gauss-Point properties
         ms        = gp.ms;
         qs        = gp.qs;
         q_cr      = gp.q_cr;
         obj.calc_grads(8,1,gp.Rn_list(:,:,:,gp.iel),
         %          tau_l     =
         % Values from previous step
         n.ep    = obj.ep   (:,:, gp.i, gp.iel, step);
         n.eEff  = obj.eEff (     gp.i, gp.iel, step);
         
         % Trial equivalent stress
         np1.e      = obj.e;
         np1.S      = 2*G* (np1.e - n.ep); % Assuming no plasticity
         np1.sEquiv = sqrt( (3/2)*sum( dot(np1.S, np1.S) ));% Assuming no plasticity
         
         % Effective plastic strain increment
         np1.deEff =...
            obj.getEffPlasStrnInc(G, Y, n.eEff, np1.sEquiv, e0, nExp, dt, edot0, m);
         np1.eEff  = n.eEff + np1.deEff;
         
         % Decide whether solving plastic or elastic tangential
         if (np1.sEquiv*np1.deEff>0)
            np1.dep = (3/2)*np1.deEff* (np1.S/np1.sEquiv);
            
            c1 = 1 - 3*G*np1.deEff/np1.sEquiv;
            c2 =...
               3*G/np1.sEquiv + c1*( 1/(nExp*(e0+n.eEff+np1.deEff)) + 1/(m*np1.deEff) );
            c3 = 4.5*G* (np1.deEff - 1/c2)/(c1*(np1.sEquiv)^3);
         else
            np1.dep = (3/2)*np1.deEff* (np1.S/1);
            
            c1 = 1;
            c3 = 0;
         end
         np1.ep = n.ep + np1.dep;
         
         % Save states as they'll be sent to stress function
         obj.eEff(     gp.i, gp.iel, step+1)= np1.eEff;
         obj.ep  (:,:, gp.i, gp.iel, step+1)= np1.ep;
         
         if     obj.ndm == 2
            S  = [np1.S(1,1) np1.S(2,2) 0          np1.S(1,2) 0          0         ];
         elseif obj.ndm == 3
            S  = [np1.S(1,1) np1.S(2,2) np1.S(3,3) np1.S(1,2) np1.S(2,3) np1.S(1,3)];
         end
         % Compute D
         D =K*I4_bulk          +... % Elastic bulk
            c1*   (2*G)*I4_dev +... % Elastic deviatoric
            c1*c3*(2*G)*(S'*S);     % Plastic
         
         ctan = reshape(D([1,4,6,4,2,5,6,5,3],[1,4,6,4,2,5,6,5,3]),3,3,3,3);
         if obj.ndm == 2
            D =D([1,2,4],[1,2,4]);
         end
      end
      %% Element K
      function Kel = computeK_el(~, Kel, gp, ~)
         % Definitions
         B=gp.B;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         Kel = Kel + (B'*gp.D*B) *gp.J *gp.w;
      end
      %% Element Fint
      function Fint = computeFint(~, gp, el)
         Fint = el.Fint + (gp.B'*gp.sigma) *gp.J *gp.w;
      end
      
      
   end
   
   methods (Static)
      
      function deEff = getEffPlasStrnInc(G, Y, eEff, sEquiv, e0, nExp, dt, edot0, m)
         deEff = 10^(-15);
         err   = Y;
         tol   = 10^(-06)*Y;
         
         if (sEquiv*edot0 == 0)
            deEff = 0;
         else
            while (err>tol)
               c = (1+(eEff+deEff)/e0)^(1/nExp) * (deEff/(dt*edot0))^(1/m);
               f = sEquiv/Y - 3*G*(deEff/Y)- c;
               dfde = -3*G/Y - c*(1/(nExp*(eEff+deEff+e0)) + 1/(m*deEff));
               enew = deEff - f/dfde;
               if (enew<0)
                  %        e must be >0, so if new approx to e <0 the solution
                  %        must lie between current approx to e and zero.
                  deEff = deEff/10;
               else
                  deEff = enew;
               end
               err = abs(f);
            end
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
      
      function gradFes = mm10_calc_grads(ngp, geonl, rot_blk, jac, Rps, ~)         
         Rt = zeros(ngp,3,3);
         
         % Get R components and stick in the right place
         if (geonl)
            jacinv = inv(jac);
            for i = 1:ngp
               Rt(i,1:3,1:3) = jacinv*reshape(Rps(1:9,i),3,3)*reshape(rot_blk(1:9,i),3,3)';
            end
         else
            Rt = reshape(Rps',ngp, 3,3);
         end
         
         intermat = 1/sqrt(3)*[...
            -1 +1 -1 +1 -1 +1 -1 +1
            -1 -1 +1 +1 -1 -1 +1 +1
            -1 -1 -1 -1 +1 +1 +1 +1];
         %       Vector:
         LHS = reshape(Rt,ngp,9)';
         RHS = LHS/intermat;
         
         grads = reshape(RHS,3,3,3);
         
         gradFes = repmat(reshape(grads(1:3,1:3,1:3), 27, 1),1,ngp); 
      end
   end
end  