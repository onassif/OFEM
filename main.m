clear; close all
% Import project folders
addpath( ...
  'CPmodels\', ...
  'Hardenings\', ...
  'Materials\', ...
  'Shapes\', ...
  'Slips\', ...
  'Utilities\');

% Input Dialog
run(ReadInput())

initialize;
tic
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
    globl.M    = sparse(num.eq, num.eq);
    globl.Fint = zeros( num.eq,1);
    % Send new U to the element object
    el.Ures_glb = globl.U(1:num.np*num.ndm) - el.U_glb_n;
    if el.iter ==0
      el.U_glb_n  = el.U_global;
    end
    el.U_global = globl.U;
    for iel =1:num.el
      %% Start elements loop
      % Get information related to current element:
      el.i = iel;        gp.iel = iel;        % element number
      gp.U = el.Umt;     gp.U_n = el.Umt_n;   gp.dU = el.w;% element unknowns (array form, n and n+1)
      el.Fint = zeros(length(el.indices),1);
      el.K    = zeros(length(el.indices), length(el.indices));
      el.M    = zeros(length(el.indices), length(el.indices));

      for igp = 1:el.mat{el.im}.ngp;     gp.i = igp;
        %%   Start Loop over Gauss points
        %%%   Compute gp K, F and in case of dynamics: M
        [gp, el, el.mat{el.im}] = el.mat{el.im}.computeKFM(gp, el, step);
            
        % store states
        if igp <= num.gp
          hist.eps( :,igp,iel) = gp.eps;
          if size(hist.stre,1) == 3 && size(gp.sigma,1) == 6
            gp.sigma = gp.sigma([1,2,4]);
          end
          hist.stre(:,igp,iel) = gp.sigma;
          hist.D( :,:,igp,iel) = gp.D;
        end
      end
         
      % Assemble to global arrays
      [globl.K, globl.M, globl.Fint] = Assemble(globl.K, globl.M, globl.Fint, el);
    end
      
    %% Finished gauss points loop, back to NR loop
    %%%   7i. Fext and apply constrains
    [globl, Fext, G, knwndU, rmIndc]  =  ApplyConstraints_and_Loads(...
         NR.mult   ,...
         globl     ,...
         Fext      ,...
         globl.U   ,...
         inpt      ,...
         num.ndm   ,...
         step      ,...
         NR.iter   ,...
         finiteDisp,...
         mixed);
      
    %%%   8i. dU and update Ui
    [dU, NR] = NR.solveSystem(globl, G, rmIndc, knwndU, extrapolate);
                        
    globl.w = globl.w*(NR.iter>0) + dU;
    globl.U = globl.U + dU;
      
    NR.correction = norm(dU)/norm(globl.w);
    NR.residual   = norm(G)/(num.np*num.ndof);
    el.w_global   = globl.w;
                        
    % Print iteration information:
    fprintf(['' ...
      'step: %4i\t ' ...
      'iteration:%2i\t ' ...
      'correction: %.10f\t ' ...
      'residual: %.10f\n'],...
      step, ...
      NR.iter, ...
      NR.correction, ...
      NR.residual);

    NR.iter = NR.iter + 1;
    el.iter = NR.iter;
      
    %%%   9i. History arrays
    hist.resid(NR.iter) = NR.residual;
    % write to file
    if (NR.correction < NR.tol)
      WriteReadHistory(...
        hist        ,...
        num         ,...
        step        ,...
        Fext        ,...
        globl.U + dU,...
        nodes+reshape(globl.U(1:num.ndm*num.np),num.ndm,num.np)');
    end
      
    if (NR.iter > NR.max_iter)
      error('>>>> Error!! Exceeded number of NR iterations permissible, exiting..');
    end
  end
  %% Finished NR loop, back to time step loop
end
t = toc
%%%     10. Post Processing
[hist, num] = WriteReadHistory(hist, num, nodes);
PostProcess(hist, num, gp);