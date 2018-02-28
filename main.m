%     Program CE691
clear; close all

% Dialog
run(ReadInput())
initialize

for step=1:num.steps % Steps loop
    %% Start time loop
    
    %%% 1n - Load multiplier
    NR.step = step;
    NR.mult = step/num.steps;
    NR.iter = 0; fprintf('\n');
    NR.correction = inf;
    globl.w    = zeros(num.eq,1);
    
    while (NR.correction > NR.tol)
        %% Start NR loop
        % Clear global K and Fint
        globl.K    = sparse(num.eq, num.eq);
        globl.Fint = zeros(num.eq,1);
        
        % Send new U to the element object
%         el.U_global = globl.w;
        el.U_global = globl.U;
        for iel =1:num.el
            %% Start elements loop
            % Get information related to current element:
            el.i = iel;             % element number
            coor = el.coor;         % element coordinates
            Umt  = el.Umt;          % element unknowns (array form)
            props= el.props;        % element material properties
            
            % clear previous values of elemental K and Fint
            el.K    = 0;
            el.Fint = 0;
            
            for igp = 1:num.gp;     gp.i = igp;
                %%   Start Loop over Gauss points
                
                %%%   2gp. Gauss points geometry-related values
                gp = Compute_gp_info(gp, coor, Umt, iel, num);
                
                %%%   3gp. Strain tensor
                [gp.eps, mat] = mat.computeStrain(gp, el, step);
                
                %%%   4gp. Tangential stifness
                [gp.D, gp.ctan, mat] = mat.computeTangentStiffness(gp, step);
                
                %%%   5gp. Stress
                [gp.sigma, mat] = mat.computeCauchy(gp, step);
                
                %%%   6gp. K
                el.K            = mat.computeK_el(el.K, gp, num.gp);
                
                %%%   7gp. Fint
                el.Fint         = mat.computeFint(gp, el);
                
                % store states
                hist.eps (:,igp,iel) = gp.eps;
                hist.stre(:,igp,iel) = gp.sigma;
                hist.ctan(:,:,:,:,igp,iel) = gp.ctan;
            end
            
            % Assemble to global arrays
            [globl.K, globl.Fint] = Assemble(globl.K, globl.Fint, el);
        end
        
        %% Finished gauss points loop, back to NR loop
        
        %%%   8i. Fext and apply constrains
        [globl.K, Fext, globl.Fint, G]  =  ApplyConstraints_and_Loads(...
            NR.mult, globl.K, Fext, globl.Fint, globl.U, inpt, num.ndm);
        %             NR.mult, globl.K, Fext, globl.Fint, inpt, num.ndof);        
        
        %%%   10i. dU and update Ui
        dU      = globl.K\G .* ~(mat.linear==1 && NR.iter>0);
        globl.w = globl.w + dU;
        globl.U = globl.U + dU;
        
        NR.correction = norm(dU)/norm(globl.w); 
        NR.residual   = norm(G)/(num.np*num.ndof);
        
        %           Print iteration information:
        NR.iter = NR.iter + 1;
        fprintf('step: %4.0d\t iteration:%2.0d\t correction: %.10f\t residual: %.10f\n',...
            step, NR.iter, NR.correction, NR.residual);
        
        
        
        %%%   11i. History arrays
        hist.resid(NR.iter) = NR.residual;
        % write to file
        if (NR.correction < NR.tol)
%             globl.U = globl.U + globl.w;
            
            WriteReadHistory(hist, num, step, globl.Fint, globl.U + dU,...
                nodes+reshape(globl.U(1:num.ndm*num.np),num.ndm,num.np)');
        end
        
        
        if (NR.iter > NR.max_iter)
            error('>>>>> Error!! Exceeded the number of NR iterations permissible, exiting...');
        end
    end
    %% Finished NR loop, back to time step loop
end
%%%     12. Post Processing
[hist, num] = WriteReadHistory(hist, num, nodes);
PostProcess(hist, num, gp);