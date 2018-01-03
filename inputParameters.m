    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input Start
    coor = [0 0
            1 0
            1 1
            0 1];
    
    BC_T = {'x', 0, 'u', 0;...
%             'node',1,'v',0};
            'y', 0, 'v', 0}; 
        
    FORCE_Tx = {'x', 1, 'u', 0.6*19.0652152};       
         
    NR_tol = 1e-9;%1e-11;
    max_iter = 100;%20
    n_steps = 5;
    props = [40, 40]; % props = [mu lambda];
    eltype = 'Q4';

    

    [nodes, elements, nen, ngp, numnp, numel, ndm, BC, FORCE] =...
        Generate_mesh2(eltype, coor, 3, 3, true, BC_T, FORCE_Tx);
    nummat = 1;
    material = 2; % Elastic
    ndof = 2;
    numeq = ndof*numnp;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% User input End