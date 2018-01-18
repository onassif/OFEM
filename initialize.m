%   voigt size
numstr = (ndm*ndm + ndm) /2;

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

num.el          = numel;
num.np          = numnp;
num.nen         = nen;
num.ndof        = ndof;
num.eq          = numeq;
num.gp          = ngp;
num.steps       = n_steps;
num.str         = numstr;
num.ndm         = ndm;

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
% Material-related
if (material == 1)
    mat = Elastic(ndm, ndof);
elseif   (material == 2)
    mat = PlaneStrain(ndm, ndof);
elseif (material == 3)
    mat = HyperNeo(ndm, ndof);  
end

% gp-related 
if(eltype == 'Q4')
    gp = Q4(num, mat.finiteDisp);
elseif (eltype == 'Q9')
    gp = Q9;
elseif    (eltype == 'T3')
    gp = T3(numel, mat.finiteDisp);
elseif    (eltype == 'T6')
    gp = T6;
elseif    (eltype == 'Q8')
    gp = Q8(num, mat.finiteDisp);
end

% Element-related
el = Elements(elements, nodes, numel, numnp, ndm, ndof, nen, props, globl.U);

inpt.BC         = BC;
inpt.FORCE      = FORCE;

mon_str={'Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sep','Oct','Nov','Dec'};
current = fix(clock);
hist.start= sprintf('hist%shist_%d-%d-%s_%d:%d:%d',...
    filesep, current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start); 

clear('BC','BC_T', 'FORCE', 'numel','numnp', 'nen', 'ndof', 'numeq', 'ngp',...
    'n_steps', 'props', 'numstr', 'ndm', 'count', 'mon_str','current', 'material')
            