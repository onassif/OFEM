%   voigt size
numstr = (ndof*ndof + ndof) /2;

%   Assign zero values to Arrays
dU              = zeros(numeq,1);
Fext            = zeros(numeq,1);
G               = zeros(numeq,1);
globl.U         = zeros(numeq,1);

step            = cell (n_steps+1,1);

hist.eps         = zeros(numstr, ngp, numel, 'single');
hist.stre        = zeros(numstr, ngp, numel, 'single');
hist.ctan        = zeros(3,3,3,3,ngp, numel, 'single');
hist.resid       = zeros(max_iter,1,'single');
hist.conn        = uint32(elements(:,1:nen));

hist.nconn       = zeros(numnp,nen+1);
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
%%%%%%%%%%%%%%%%%%%%% Create objects:
% gp-related 
if(eltype == 'Q4')
    gp = Q4(numel);
elseif (eltype == 'Q9')
    gp = Q9;
elseif    (eltype == 'T3')
    gp = T3;
elseif    (eltype == 'T6')
    gp = T6;
end

% Element-related
el = Elements(elements, nodes, numel, numnp, ndof, ndof, nen, props, globl.U);

% Material-related
if (material == 1)
    mat = HyperNeo;  
end
    
inpt.BC         = BC;
inpt.FORCE      = FORCE;

num.el          = numel;
num.np          = numnp;
num.nen         = nen;
num.ndof        = ndof;
num.eq          = numeq;
num.gp          = ngp;
num.steps       = n_steps;
num.str         = numstr;
num.ndm         = ndm;

mon_str={'Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sep','Oct','Nov','Dec'};
current = fix(clock);
hist.start= sprintf('hist_%d-%d-%s_%d:%d:%d',...
    current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start); 

clear('BC','BC_T', 'FORCE', 'numel','numnp', 'nen', 'ndof', 'numeq', 'ngp',...
    'n_steps', 'props', 'numstr', 'ndm', 'count', 'mon_str','current', 'material')
            