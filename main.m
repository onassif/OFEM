clear; close all
% Dialog
run(ReadInput())
initialize;
mat = mat; % Becauce MATLAB is crazy

for step=1:num.steps % Steps loop
   %% Start time loop
   %%% 1n - Load multiplier
   NR.step = step;
   NR.mult = step/num.steps;
   NR.iter = 0; fprintf('\n');
   NR.correction = inf;
   el.iter = 0; % To save U_n when overwriting with U_n+1
   
   while (NR.correction > NR.tol)
      %% Start NR loop
      % Clear global K and Fint
      globl.K    = sparse(num.eq, num.eq);
      globl.Fint = zeros (num.eq,1);
      dU         = zeros (num.eq,1);
      % Send new U to the element object
      el.U_global = globl.U;
      for iel =1:num.el
         %% Start elements loop
         % Get information related to current element:
         el.i   = iel;        gp.iel = iel;        % element number
         coor   = el.coor;                         % element coordinates
         gp.U   = el.Umt;     gp.U_n = el.Umt_n;   gp.dU = el.w;% element unknowns (array form, n and n+1)
         props  = el.props;                        % element material properties
         if isa(mat{el.im},'DG') && firstInstance
            firstInstance = false;
            tempK = globl.K;
         end
         for igp = 1:mat{el.im}.ngp;     gp.i = igp;
            %%   Start Loop over Gauss points
            %%%   2gp. Strain tensor
            [gp.eps, mat{el.im}]        = mat{el.im}.computeStrain(gp, el, step);
            
            %%%   3gp. Tangential stifness
            [gp.D, gp.ctan, mat{el.im}] = mat{el.im}.computeTangentStiffness(gp, el, step);
            
            %%%   4gp. Stress
            [gp.sigma, mat{el.im}]      = mat{el.im}.computeCauchy(gp, el, step);
            
            %%%   5gp. K
            el.K    = mat{el.im}.computeK_el(gp, el, step);
            
            %%%   6gp. Fint
            el.Fint = mat{el.im}.computeFint(gp, el, step);
            
            % store states
            hist.eps (:,igp,iel)       = gp.eps;
            hist.stre(:,igp,iel)       = gp.sigma;
            hist.ctan(:,:,:,:,igp,iel) = gp.ctan;
         end
         
         % Assemble to global arrays
         [globl.K, globl.Fint] = Assemble(globl.K, globl.Fint, el);
      end
      
      %% Finished gauss points loop, back to NR loop
      %%%   7i. Fext and apply constrains
      [globl.K, Fext, globl.Fint, G, knwndU, rmIndc]  =  ApplyConstraints_and_Loads(...
         NR.mult, globl.K, Fext, globl.Fint, globl.U, tempK, inpt, num.ndm, step, NR.iter, finiteDisp);
      
      %%%   8i. dU and update Ui
      if el.iter == 0 && step > 1 && extrapolate
         dU( rmIndc) = (globl.K\G+globl.w(rmIndc));
      else
         dU( rmIndc) = globl.K\G;
      end
      dU(~rmIndc) = knwndU;
                        
      globl.w = globl.w*(NR.iter>0) + dU;
      globl.U = globl.U + dU;
      
      NR.correction = norm(dU)/norm(globl.w);
      NR.residual   = norm(G)/(num.np*num.ndof);
      el.w_global = globl.w;
                        
      % Print iteration information:
      fprintf('step: %4i\t iteration:%2i\t correction: %.10f\t residual: %.10f\n',...
         step, NR.iter, NR.correction, NR.residual);
      NR.iter = NR.iter + 1; 
      el.iter = NR.iter;
      
      %%%   9i. History arrays
      hist.resid(NR.iter) = NR.residual;
      % write to file
      if (NR.correction < NR.tol)
         WriteReadHistory(hist, num, step, Fext, globl.U + dU,...
            nodes+reshape(globl.U(1:num.ndm*num.np),num.ndm,num.np)');
      end
      
      if (NR.iter > NR.max_iter)
         error('>>>> Error!! Exceeded number of NR iterations permissible, exiting..');
      end
   end
   %% Finished NR loop, back to time step loop
end
%%%     10. Post Processing
[hist, num] = WriteReadHistory(hist, num, nodes);
PostProcess(hist, num, gp);