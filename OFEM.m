%     Program CE691
clear; close all

% Dialog
inputParameters
initialize

for istep=1:num.steps % Steps loop
    %% Start time loop
    
    %%% 1n - Load multiplier
    mult = istep/num.steps;
    
    iter = 0; fprintf('\n');
    residual = inf;
    while (residual > NR_tol)
        %% Start NR loop
        % Clear global K and Fint
        globl.K = sparse(num.eq, num.eq);          globl.Fint = zeros(num.eq,1);
        
        % Send new U to the element object
        el.U_global = globl.U;

        for iel =1:num.el          
            %% Start elements loop
            % Get information related to current element:
            el.i = iel;             % element number
            coor = el.coor;         % element coordinates
            Umt  = el.Umt;          % element unknowns (array form)
            Uvc  = el.Uvc;          % element unknowns (vector form)
            props= el.props;        % element material properties
            
            % clear previous values of elemental K and Fint
            el.clear_K_Fint();
            
            for igp = 1:num.gp;     gp.i = igp;
                %%   Start Loop over Gauss points
                
                %%%   2gp. Gauss points geometry-related values
                gp = Compute_gp_info(gp, coor, Umt, iel);
                
                %%%   3gp. Strain tensor
                gp.eps = gp.B * Uvc;
                
                %%%   4gp. Stress
                gp.sigma = Compute_cauchy(num.ndm, gp, props);
                
                %%%   5gp. Tangential stifness
                [gp.D, gp.ctan] = Compute_tangentstiffness(gp, props);
                
                %%%   6gp. K
                el.K = Compute_Kel(el.K, gp, num.ndof, num.gp);
                
                %%%   7gp. Fint
                el.Fint =  el.Fint + (gp.B'*gp.sigma)*gp.j *gp.w;
                
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
        [globl.K, Fext, globl.Fint]  =  ApplyConstraints_and_Loads(...
            mult, globl.K, Fext, globl.Fint, inpt, num.ndof);
        
        %%%   9i. Residual and normalized residual
        G = Fext - globl.Fint;
        residual = norm(G)/norm(Fext);

        %%%   10i. dU and update Ui
        dU = globl.K\G;
        globl.U = globl.U + dU;
        
        %   Print iteration information:
        iter = iter + 1;
        fprintf(['step: %d\t iteration: %d\t residual: ',...
            '[abs: %.10f\t rel: %.10f]\n'],istep, iter, norm(G),residual);...
        
        %%%   11i. History arrays
        hist.resid(iter) = residual;
        % write to file 
        if (residual < NR_tol)
            WriteReadHistory(hist, num, istep, Fext, globl.U + dU,...
                nodes+reshape(globl.U,num.ndof,num.np)');
        end
        
        
        if (iter > max_iter)
            error('>>>>> Error!! Exceeded the number of NR iterations permissible, exiting...');
        end
    end
    %% Finished NR loop, back to time step loop
end
%%%     12. Post Processing
[hist, num] = WriteReadHistory(hist, num, nodes);
PostProcess(hist, num, gp);