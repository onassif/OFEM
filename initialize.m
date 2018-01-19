%   voigt size
numstr = (ndm*ndm + ndm) /2;

%   Assign zero values to Arrays
dU              = zeros(numeq,1);
Fext            = zeros(numeq,1);
G               = zeros(numeq,1);
globl.U         = zeros(numeq,1);

step            = cell (n_steps+1,1);

hist.eps        = zeros(numstr, ngp, numel, 'single');
hist.stre       = zeros(numstr, ngp, numel, 'single');
hist.ctan       = zeros(3,3,3,3,ngp, numel, 'single');
hist.resid      = zeros(max_iter,1,'single');
hist.conn       = uint32(elements(:,1:nen));
hist.nconn      = zeros(numnp,nen+1);
for i=1:numnp
    count = 0;
    for j =1:numel
        if( find(hist.conn(j,:)==i) )
            count = count+1;
            hist.nconn(i,count) = j;
        end
    end
    hist.nconn(i,nen+1) = count;
end

num.el          = numel;
num.np          = numnp;
num.nen         = nen;
num.ndof        = ndof;
num.eq          = numeq;
num.gp          = ngp;
num.steps       = n_steps;
num.str         = numstr;
num.ndm         = ndm;


%%%%%%%%%%%%%%%%%%%%% Create objects:
% Material-related
if     (material == 1)
    mat = Elastic(num, props);
elseif (material == 2)
    mat = PlaneStrain(num, props);
elseif (material == 3)
    mat = HyperNeo(num, props);  
elseif (material == 4)
    mat = ClassicPlasticityRI(num, props);
    hist.dt          = n_steps/total_time;
    hist.ttime       = total_time;
end
    

% gp-related 
if     (eltype == 'Q4')
    gp = Q4(num, mat.finiteDisp);
elseif (eltype == 'Q9')
    gp = Q9;
elseif (eltype == 'T3')
    gp = T3(numel, mat.finiteDisp);
elseif (eltype == 'T6')
    gp = T6;
elseif (eltype == 'Q8')
    gp = Q8(num, mat.finiteDisp);
end

% Element-related
el = Elements(elements, nodes, num, props, globl.U);

% NR-related
NR = NewtonRaphson(NR_tol, max_iter, n_steps);

inpt.BC         = BC;
inpt.FORCE      = FORCE;

mon_str={'Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sep','Oct','Nov','Dec'};
current = fix(clock);
hist.start= sprintf('hist%shist_%d-%d-%s_%d:%d:%d',...
    filesep, current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start); 

clearvars -except dU el elements Fext globl gp hist inpt mat nodes NR num           