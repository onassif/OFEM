classdef NewtonRaphson
   %NewtonRaphson Some parameters used for NR
   %   Detailed explanation goes here
   properties
      tol;
      max_iter;
      mult
      iter
      residual
      correction
      dyn
      step=1
      
      a;
      v;
      u;
   end
   properties (SetAccess=private)
      numsteps
      time=0;
      fctr=0;
      dt
      dyn_steps
   end
   
   methods
      %% Construct
      function ob = NewtonRaphson(NR_tol, max_iter, num, dynamic, varargin)
         ob.tol        = NR_tol;
         ob.max_iter   = max_iter;
         ob.numsteps   = num.steps;
         ob.mult       = zeros(num.mlLen,2);
         ob.dyn        = dynamic;
         
         if dynamic > 0
            ob.dt = num.dyn_time / num.dyn_steps;
            ob.dyn_steps = num.dyn_steps;
            
            ob.a = zeros(num.redEq, num.dyn_steps+1);
            ob.v = zeros(num.redEq, num.dyn_steps+1);
            ob.u = zeros(num.redEq, num.dyn_steps+1);
         end
         switch nargin
            case 5
               ob.time     =[zeros(1,size(varargin{1},2)); varargin{1}];
               ob.time(:,2)=cumsum(ob.time(:,2));
            case 6
               ob.time     =[zeros(1,size(varargin{1},2)); varargin{1}];
               ob.time(:,2)=cumsum(ob.time(:,2));
               ob.fctr=[zeros(num.mlLen,1) varargin{2}(:,:,1), zeros(num.mlLen,1), varargin{2}(:,:,2)];
               ob.fctr=reshape(ob.fctr, num.mlLen , size(varargin{2},2)+1,2);
         end
      end
      %% set the multiplier
      function ob = set.mult(ob,value)
         if (~ob.time | length(ob.fctr)==1)%if user didn't pass time or factors or Initialization
            if length(value)==1
               ob.mult(:,1) = value;
               ob.mult(:,2) = value;
            else % Initialization
               ob.mult(:,1) = value(:,2);
               ob.mult(:,2) = value(:,2);
            end
         else
            row = find(ob.step>ob.time(:,2),1,'last')+1;
            ncrmnt = squeeze(...
               (ob.fctr(:,row,:)-ob.fctr(:,row-1,:))/...
               (ob.time(  row,2)-ob.time(  row-1,2)) );
            
            ob.mult(:,1) = ob.mult(:,1) + ncrmnt;
            ob.mult(:,2) = ob.mult(:,2) + ncrmnt;
         end
      end
      %% Solve the system of linear equations
      function [dU, ob] = solveSystem(ob, globl, G, rmIndc, knwndU, extrapolate)
         dU = zeros(size(rmIndc));
         
         if ob.iter==0 && ob.step>1 && extrapolate
            dU(rmIndc) = globl.K\G + globl.U(rmIndc);
         else
            dU(rmIndc) = globl.K\G;
         end
         dU(~rmIndc) = knwndU;
         
         if ob.dyn >0
            beta1 = 0.5;
            beta2 = 0.0;
            
            if ob.dyn == 1
               invMK = inv(globl.M + 0.5*beta2*ob.dt*ob.dt*globl.K);
            elseif ob.dyn == 2
               invMK = 1./diag(globl.M + 0.5*beta2*ob.dt*ob.dt*globl.K);
            end
            
            ob.a(:,1) = globl.M \ (G - globl.K*ob.u(:,1));
            for i = 1:ob.dyn_steps
               if ob.dyn == 1
                  ob.a(:,i+1) = invMK *...
                     (G - globl.K*(ob.u(:,i) + ob.dt*ob.v(:,i) + 0.5*ob.dt^2*(1.-beta2)*ob.a(:,i)) );
               elseif ob.dyn == 2
                  ob.a(:,i+1) = invMK .*...
                     (G - globl.K*(ob.u(:,i) + ob.dt*ob.v(:,i) + 0.5*ob.dt^2*(1.-beta2)*ob.a(:,i)) );
               end
               ob.v(:,i+1) = ob.v(:,i) + ob.dt*(1.-beta1)*ob.a(:,i) + ob.dt*beta1*ob.a(:,i+1);
               ob.u(:,i+1) = ob.u(:,i) + ob.dt*ob.v(:,i) + (1.-beta2)*0.5*ob.dt^2*ob.a(:,i) +...
                  0.5*beta2*ob.dt^2*ob.a(:,i+1) ;
            end
         end
         
      end
   end
end
