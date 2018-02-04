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
    NR.residual = inf;
    while (NR.residual > NR.tol)
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
            el.K    = 0;    
            el.Fint = 0;
            
            for igp = 1:num.gp;     gp.i = igp;
                %%   Start Loop over Gauss points
                
                %%%   2gp. Gauss points geometry-related values
                gp = Compute_gp_info(gp, coor, Umt, iel, num);
                 
                %%%   3gp. Strain tensor
                gp.eps = gp.B * Uvc;
                
                %%%   4gp. Tangential stifness
                [gp.D, gp.ctan, mat] = mat.Compute_tangentstiffness(gp);
                
                %%%   5gp. Stress
                gp.sigma        = mat.Compute_cauchy(gp);
                
                %%%   6gp. K
                el.K            = mat.Compute_Kel(el.K, gp, num.gp);
                
                %%%   7gp. Fint
                el.Fint         = mat.Compute_Fint(el.Fint, gp);
                
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
            NR.mult, globl.K, Fext, globl.Fint, inpt, num.ndof, NR.iter==0);
        
        %%%   9i. Residual and normalized residual
        G = Fext - globl.Fint;
        NR.residual = norm(G)/norm(Fext);

        %%%   10i. dU and update Ui
        dU = globl.K\G;
        globl.U = globl.U + dU;
        
        %   Print iteration information:
        NR.iter = NR.iter + 1;
        fprintf(['step: %d\t iteration: %d\t residual: ',...
            '[abs: %.10f\t rel: %.10f]\n'],step, NR.iter, norm(G),NR.residual);
        
        %%%   11i. History arrays
        hist.resid(NR.iter) = NR.residual;
        % write to file 
        if (NR.residual < NR.tol)
            WriteReadHistory(hist, num, step, Fext, globl.U + dU,...
                nodes+reshape(globl.U,num.ndof,num.np)');
        end
        
        
        if (NR.iter > NR.max_iter)
            error('>>>>> Error!! Exceeded the number of NR iterations permissible, exiting...');
        end
    end
    %% Finished NR loop, back to time step loop
end
%%%     12. Post Processing
[hist, num] = WriteReadHistory(hist, num, nodes);
% PostProcess(hist, num, gp);