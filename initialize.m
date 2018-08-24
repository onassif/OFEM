%   voigt size
numstr = (ndm*ndm + ndm) /2;

%   Assign zero values to Arrays
dU         = zeros(numeq,1);
Fext       = zeros(numeq,1);
G          = zeros(numeq,1);
globl.U    = zeros(numeq,1);
globl.w    = zeros(numeq,1);
globl.Fint = zeros(numeq,1);
tempK      = [];
step       = cell (n_steps+1,1);

hist.eps   = zeros(numstr, ngp, size(elements,1), 'single');
hist.stre  = zeros(numstr, ngp, size(elements,1), 'single');
hist.ctan  = zeros(3,3,3,3,ngp, size(elements,1), 'single');
hist.resid = zeros(max_iter,1,'single');
hist.conn  = uint32(elements(:,1:nen));
hist.nodes = nodes;
maxSharedNode = 0;
realnumnp = size(nodes,1);
for i =1: numnp
  occurances =  length(find(hist.conn==i));
  if occurances > maxSharedNode
     maxSharedNode = occurances;
  end
end
hist.nconn = zeros(numnp, maxSharedNode+1);
for i=1:numnp
   count = 0;
   for j =1:numel
      if( find(hist.conn(j,:)==i) )
         count = count+1;
         hist.nconn(i,count) = j;
      end
   end
   hist.nconn(i,end) = count;
end

num.el    = numel;
num.np    = numnp;
% num.rnp   = realnumnp;
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
if exist('hardType','var')   
   switch slipType
      case {'FCC', 'fcc'}
         slip = FCC();
   end
end

% Material-related
mat = cell(length(material),1);
for i=1:length(material)
   if length(material) == 1
      prps = props;
   else
      prps = props{i};
   end
   switch material(i)
      case 1
         mat{i} = Elastic(num, prps);
      case 2
         mat{i} = HypoElastic(num, prps, ident.threeD.second);
      case 3
         mat{i} = HyperNeo(num, prps);
      case 4
         mat{i} = ViscoPlastic(num, prps, time, ident.threeD.second);
      case 5
         mat{i} = PlasticRI(num, prps, ident.threeD.second);
      case 6
         mat{i} = MixedElasticPlaneStrain(num, prps);
      case 7
         mat{i} = MTS(num, prps, hardProps, slip, time, ident.threeD.second);
      case 8
         mat{i} = DG(num, prps, props); 
   end
end

% gp-related
switch eltype
   case 'Q4'
    gp = Q4(mat{1}.finiteDisp, hist);
   case 'Q9'
    gp = Q9(mat{1}.finiteDisp, hist);
   case 'T3'
    gp = T3(mat{1}.finiteDisp, hist);
   case 'T6'
    gp = T6(mat{1}.finiteDisp, hist);
   case 'Q8'
    gp = Q8(mat{1}.finiteDisp, hist);
   case 'Q8Crys'
    gp = Q8Crys(mat{1}.finiteDisp, mat{1}.hardProps.angles, slip, hist);  
end

% Element-related
el = Elements(elements, nodes, num, props, globl.U, hist);
el.w_global = globl.w;

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
el.ulres = getInitU(1/num.steps, inpt, el);

num.el = size(elements,1);

mon_str={'Jan','Feb','Mar','Apr','May','June','Jul','Aug','Sep','Oct','Nov','Dec'};
current = fix(clock);
hist.start= sprintf('%s%shist%shist_%d-%d-%s_%d.%d.%d',...
   pwd,filesep,filesep, current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start);
firstInstance = true;
clearvars -except dU el elements Fext globl gp hist inpt mat nodes NR num hard slip tempK firstInstance