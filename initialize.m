%   voigt size
numstr = (ndm*ndm + ndm) /2;

%   Assign zero values to Arrays
dU         = zeros(numeq,1);
Fext       = zeros(numeq,1);
G          = zeros(numeq,1);
globl.U    = zeros(numeq,1);
globl.w    = zeros(numeq,1);
globl.Fint = zeros(numeq,1);

hist.eps   = zeros(numstr, ngp, size(elements,1), 'single');
hist.stre  = zeros(numstr, ngp, size(elements,1), 'single');
hist.D     = zeros(   6,6, ngp, size(elements,1), 'single');
hist.resid = zeros(max_iter,1,'single');
hist.conn  = uint32(elements(:,1:end-1));
hist.nodes = nodes;
hist.eltype= eltype;
maxSharedNode = 0;
realnumnp = size(nodes,1);
for i =1: numnp
  occurances =  length(hist.conn(hist.conn==i));
  if occurances > maxSharedNode
     maxSharedNode = occurances;
  end
end
hist.nconn = zeros(numnp, maxSharedNode+1);
for i=1:numnp
   count = 0;
   for j =1:numel
      if( sum(hist.conn(j,:)==i)>0 )
         count = count+1;
         hist.nconn(i,count) = j;
      end
   end
   hist.nconn(i,end) = count;
end

num.el    = size(elements,1);
num.els   = numel;
num.el2   = size(elements,1);
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
if exist('cpType','var')   
   switch slipType
      case {'FCC', 'fcc'}
         slip = FCC();
   end
end

% Element-related
el = Elements(elements, nodes, num, globl.U, hist);

% Material-related
el.mat = cell(length(material),1);
finiteDisp = 0;
for i=1:length(material)
   if length(material) == 1
      prps = props;
   else
      prps = props{i};
   end
   switch material(i)
      case 1
         el.mat{i} = Elastic(num, prps, ident.threeD.second);
      case 2
         el.mat{i} = HypoElastic(num, prps, ident.threeD.second);
      case 3
         el.mat{i} = HyperElastic(num, prps,  ident.threeD.second);
      case 4
         el.mat{i} = ViscoPlastic(num, prps, time, ident.threeD.second);
      case 5
         el.mat{i} = PlasticRI(num, prps, ident.threeD.second);
      case 6
         el.mat{i} = MixedElasticPlaneStrain(num, prps, ident.threeD.second);
      case 7
         el.mat{i} = CP(num, prps, cpType, hardProps, angles, slip, time, ident.threeD.second);
      case 8
         el.mat{i} = DG(num, prps, props, ident.threeD.second); 
      case 9
         el.mat{i} = DGHyper(num, prps, props, ident.threeD.second); 
      case 10
         el.mat{i} = DGCP(num, prps, props, ident.threeD.second); 
      case 11
         el.mat{i} = CP2(num, prps, cpType, hardProps, angles, slip, time, ident.threeD.second);
   end
   if el.mat{i}.finiteDisp
      finiteDisp = 1;
   end
end

% gp-related
switch eltype
   case 'Q4'
    gp = Q4(el.mat{1}.finiteDisp);
   case 'Q9'
    gp = Q9(el.mat{1}.finiteDisp);
   case 'T3'
    gp = T3(el.mat{1}.finiteDisp);
   case 'T6'
    gp = T6(el.mat{1}.finiteDisp);
   case 'T4'
    gp = T4(el.mat{1}.finiteDisp);
   case 'T10'
    gp = T10(el.mat{1}.finiteDisp);
   case 'Q8'
    gp = Q8(el.mat{1}.finiteDisp);
end
gp = gp.shapeIso();
[gp.det_dXdxi_list, gp.dNdX_list, gp.dXdxi_list] = shapeRef(...
   nodes, elements(1:numel,1:nen), gp.dNdxi_list);

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
hist.start= sprintf('%s%shist%shist_%d-%d-%s_%d.%d.%d',...
   pwd,filesep,filesep, current([1,3]), mon_str{current(2)}, current(4:6));
mkdir(hist.start);
firstInstance = true;

if ~exist('extrapolate','var')
   extrapolate = 0;
end
clearvars -except dU el elements Fext globl gp hist inpt mat nodes NR num hard slip tempK ...
   firstInstance finiteDisp extrapolate