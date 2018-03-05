%   voigt size
numstr = (ndm*ndm + ndm) /2;

%   Assign zero values to Arrays
dU         = zeros(numeq,1);
Fext       = zeros(numeq,1);
G          = zeros(numeq,1);
globl.U    = zeros(numeq,1);
globl.w    = zeros(numeq,1);
globl.Fint = zeros(numeq,1);
step       = cell (n_steps+1,1);

hist.eps   = zeros(numstr, ngp, numel, 'single');
hist.stre  = zeros(numstr, ngp, numel, 'single');
hist.ctan  = zeros(3,3,3,3,ngp, numel, 'single');
hist.resid = zeros(max_iter,1,'single');
hist.conn  = uint32(elements(:,1:nen));
hist.nconn = zeros(numnp,nen+1);
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

num.el    = numel;
num.np    = numnp;
num.nen   = nen;
num.ndof  = ndof;
num.eq    = numeq;
num.gp    = ngp;
num.steps = n_steps;
num.str   = numstr;
num.ndm   = ndm;
num.BC    = BC(end,end);
num.FORCE = FORCE(end,end);

%%%%%%%%%%%%%%%%%%%%% Identities:
run('identities.m');

%%%%%%%%%%%%%%%%%%%%% Create objects:
% Hardening-related
% if exist('hardType','var')
%     if strcmp(hardType,'Mixed')
%         hard = ConstantMixedHardening(hardprops);
%     end
% end

% Material-related
switch material
    case 1
        mat = Elastic(num, props);
    case 2
        mat = HypoElastic(num, props, ident.threeD.second);
    case 3
        mat = HyperNeo(num, props);
    case 4
        mat = ViscoPlastic(num, props, time, ident.threeD.second);
    case 5
        mat = PlasticRI(num, props, ident.threeD.second);
    case 6
        mat = MixedElasticPlaneStrain(num, props);
    case 7 
%         mat = HypoElastic(num, props, ident.threeD.second);
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
if exist('time','var') && exist('fctr','var')
    NR = NewtonRaphson(NR_tol, max_iter, num, time, fctr);
elseif exist('time','var')
    NR = NewtonRaphson(NR_tol, max_iter, num, time);
else
    NR = NewtonRaphson(NR_tol, max_iter, num);
end

inpt.BC    = BC;
inpt.FORCE = FORCE;

mon_str={'Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sep','Oct','Nov','Dec'};
current = fix(clock);
hist.start= sprintf('hist%shist_%d-%d-%s_%d:%d:%d',...
    filesep, current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start);

clearvars -except dU el elements Fext globl gp hist inpt mat nodes NR num