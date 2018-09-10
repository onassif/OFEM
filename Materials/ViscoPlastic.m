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
      linear = false;
      
      name = 'ViscoPlastic';
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
      Bulk
      
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
      function obj = ViscoPlastic(num, props, time, identity)
         obj.ndm  = num.ndm;
         obj.ndof = num.ndof;
         obj.ngp  = num.gp;
         obj.eEff = zeros(                 num.gp, num.el, num.steps+1);
         obj.ep   = zeros(obj.ndm,obj.ndm, num.gp, num.el, num.steps+1);
         
         for i=1:length(props)
            switch props{i,1}
               case 'E'
                  obj.E     = props{i,2};
               case 'nu'
                  obj.nu    = props{i,2};
               case 'Y'
                  obj.Y     = props{i,2};
               case 'e0'
                  obj.e0    = props{i,2};
               case 'n'
                  obj.nExp  = props{i,2};
               case 'edot0'
                  obj.edot0 = props{i,2};
               case 'm'
                  obj.m     = props{i,2};
               otherwise
                  error(['You''ve chosen Rate Independent Plastic material but ',...
                     'specified incompatible material properties, I''m disapponted']);
            end
         end
         
         obj.G    =   0.5*obj.E/(1+  obj.nu);
         obj.Bulk = (1/3)*obj.E/(1-2*obj.nu);
         
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
         
         eps= gp.B * el.Uvc;
         if obj.ndm == 2
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
      function [sigma_voigt, obj] = computeCauchy(obj, gp, ~, step)
         I = obj.I;
         G = obj.G;
         K = obj.Bulk;
         
         np1.e  = obj.e;
         np1.ep = obj.ep(:,:, gp.i, gp.iel, step+1);
         
         sigma       = 2*G* (np1.e - np1.ep) + K*sum(gp.eps(1:obj.ndm))*I;
         sigma_voigt = obj.voigtize(sigma,'col', obj.ndm);
      end
      %% Tangential stiffness
      function [D, ctan, obj] = computeTangentStiffness(obj, gp, ~, step)
         dt   = obj.dt(step);
         % Identities
         I4_dev  = obj.I4_dev;
         I4_bulk = obj.I4_bulk;
         % Material properties
         Y    = obj.Y;
         e0   = obj.e0;
         nExp = obj.nExp;
         edot0= obj.edot0;
         m    = obj.m;
         G    = obj.G;
         K    = obj.Bulk;
         
         % Values from previous step
         n.ep    = obj.ep   (:,:, gp.i, gp.iel, step);
         n.eEff  = obj.eEff (     gp.i, gp.iel, step);
         
         % Trial equivalent stress
         np1.e      = obj.e;
         np1.S      = 2*G* (np1.e - n.ep); % Assuming no plasticity
         np1.sEquiv = sqrt( (3/2)*sum( dot(np1.S, np1.S) ));% Assuming no plasticity
         
         % Effective plastic strain increment
         np1.deEff = obj.getEffPlasStrnInc(G,Y,n.eEff,np1.sEquiv,e0,nExp,dt,edot0,m);
         np1.eEff  = n.eEff + np1.deEff;
         
         % Decide whether solving plastic or elastic tangential
         if (np1.sEquiv*np1.deEff>0)
            np1.dep = (3/2)*np1.deEff* (np1.S/np1.sEquiv);
            
            c1 = 1 - 3*G*np1.deEff/np1.sEquiv;
            c2 = 3*G/np1.sEquiv + c1*( 1/(nExp*(e0+n.eEff+np1.deEff)) + 1/(m*np1.deEff) );
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
      function Kel = computeK_el(~, gp, ~, ~)
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
      
   end
end