classdef ViscoPlastic
   %Elastic 3D elastic class
   %   Detailed explanation goes here
   
   properties (Hidden, SetAccess = private)
      ndm;
      ndof;
      ngp;
      dNdX;
      finiteDisp = 0;
      
      I
      I2 = eye(3);
      I4_dev
      I4_bulk
   end
   
   properties (SetAccess = private)
      E
      nu
      Y
      e0
      nExp
      edot0
      m
      G
      K
      
      e
      ep
      eEff
      sEquiv
      
      plastic = false;
      dir;
      
      dt;
   end
   %%
   methods
      %% Construct
      function ob = ViscoPlastic(num, props, time, identity)
         ob.ndm  = num.ndm;
         ob.ndof = num.ndof;
         ob.ngp  = num.gp;
         ob.eEff = zeros(               num.gp, num.el, num.steps+1);
         ob.ep   = zeros(ob.ndm,ob.ndm, num.gp, num.el, num.steps+1);
         
         [ob.E, ob.nu, ob.Y, ob.e0, ob.nExp, ob.edot0, ob.m] = ob.getProps(props);
         
         ob.G =   0.5*ob.E/(1+  ob.nu);
         ob.K = (1/3)*ob.E/(1-2*ob.nu);
         
         ob.I       = eye(num.ndm);
         ob.I4_dev  = identity.I4_dev;
         ob.I4_bulk = identity.I4_bulk;
         
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
      function [eps, ob] = Strain(ob, gp, el, ~)
         eps  = gp.B * el.Uvc;
         ob.e = T1T2(eps,2) - (1/3)*trace(T1T2(eps,2))*ob.I;
      end
      %% Sigma & Tangential stiffness
      function [sigma_v, D, ob] = SigmaCmat(ob, gp, ~, step)
         dt   = ob.dt(step);
         % Material properties
         e0   = ob.e0;
         nExp = ob.nExp;
         edot0= ob.edot0;
         m    = ob.m;
         
         % Values from previous step
         n.ep    = ob.ep   (:,:, gp.i, gp.iel, step);
         n.eEff  = ob.eEff (     gp.i, gp.iel, step);
         
         % Trial equivalent stress
         np1.e      = ob.e;
         np1.S      = 2*ob.G* (np1.e - n.ep); % Assuming no plasticity
         np1.sEquiv = sqrt( (3/2)*sum( dot(np1.S, np1.S) ));% Assuming no plasticity
         
         % Effective plastic strain increment
         np1.deEff = ob.getEffPlasStrnInc(ob.G,ob.Y,n.eEff,np1.sEquiv,e0,nExp,dt,edot0,m);
         np1.eEff  = n.eEff + np1.deEff;
         
         % Decide whether solving plastic or elastic tangential
         if (np1.sEquiv*np1.deEff>0)
            np1.dep = (3/2)*np1.deEff* (np1.S/np1.sEquiv);
            
            c1 = 1 - 3*ob.G*np1.deEff/np1.sEquiv;
            c2 = 3*ob.G/np1.sEquiv + c1*( 1/(nExp*(e0+n.eEff+np1.deEff)) + 1/(m*np1.deEff) );
            c3 = 4.5*ob.G* (np1.deEff - 1/c2)/(c1*(np1.sEquiv)^3);
         else
            np1.dep = (3/2)*np1.deEff* (np1.S/1);
            
            c1 = 1;
            c3 = 0;
         end
         np1.ep = n.ep + np1.dep;
         
         % Save states as they'll be sent to stress function
         ob.eEff(     gp.i, gp.iel, step+1)= np1.eEff;
         ob.ep  (:,:, gp.i, gp.iel, step+1)= np1.ep;
         
         if     ob.ndm == 2
            S  = [np1.S(1,1) np1.S(2,2) 0          np1.S(1,2) 0          0         ];
         elseif ob.ndm == 3
            S  = [np1.S(1,1) np1.S(2,2) np1.S(3,3) np1.S(1,2) np1.S(2,3) np1.S(1,3)];
         end
         % Compute D
         D =ob.K*ob.I4_bulk          +... % Elastic bulk
            c1*   (2*ob.G)*ob.I4_dev +... % Elastic deviatoric
            c1*c3*(2*ob.G)*(S'*S);        % Plastic
         
         % Compute Sigma
         sigma   = 2*ob.G*(np1.e-np1.ep) + ob.K*sum(gp.eps(1:ob.ndm))*ob.I;
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
      function [E, nu, Y, e0, nExp, edot0, m] = getProps(props)
         for i=1:length(props)
            switch props{i,1}
               case 'E'
                  E     = props{i,2};
               case 'nu'
                  nu    = props{i,2};
               case 'Y'
                  Y     = props{i,2};
               case 'e0'
                  e0    = props{i,2};
               case 'n'
                  nExp  = props{i,2};
               case 'edot0'
                  edot0 = props{i,2};
               case 'm'
                  m     = props{i,2};
               otherwise
                  error(['You''ve chosen Rate Independent Plastic material but ',...
                     'specified incompatible material properties, I''m disapponted']);
            end
         end
      end
   end
end