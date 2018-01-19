function varargout = WriteReadHistory(hist, num, varargin )
%WriteHistory Writes & reads history arrays to binary files

if(nargin==6)
    step  = varargin{1};
    Fext  = varargin{2};
    d     = varargin{3};
    coor  = varargin{4};
    
    strgF = sprintf('%s%sforc_%0.4d', hist.start, filesep, step);
    strgd = sprintf('%s%sdisp_%0.4d', hist.start, filesep, step);
    strgc = sprintf('%s%scoor_%0.4d', hist.start, filesep, step);
    strgr = sprintf('%s%sresi_%0.4d', hist.start, filesep, step);
    strge = sprintf('%s%sepsi_%0.4d', hist.start, filesep, step);
    strgs = sprintf('%s%sstre_%0.4d', hist.start, filesep, step);
    strgC = sprintf('%s%sctan_%0.4d', hist.start, filesep, step);

    eps  = reshape(hist.eps , num.str, num.gp*num.el);
    stre = reshape(hist.stre, num.str, num.gp*num.el);
    ctan = reshape(hist.ctan, 3*3*3*3, num.gp*num.el);
    
    fidF = fopen(strgF,'w');
    fidd = fopen(strgd,'w');
    fidc = fopen(strgc,'w');
    fidr = fopen(strgr,'w');
    fide = fopen(strge,'w');
    fids = fopen(strgs,'w');
    fidC = fopen(strgC,'w');
    
    fwrite(fidF, Fext       , 'single');
    fwrite(fidd, d          , 'single');
    fwrite(fidc, coor       , 'single');
    fwrite(fidr, hist.resid , 'single');
    fwrite(fide, eps        , 'single');
    fwrite(fids, stre       , 'single');
    fwrite(fidC, ctan       , 'single');
    
elseif(nargin==3)
    nodes       = varargin{1};
    num.steps   = num.steps+1;
    
    hist.disp        = zeros(num.eq, num.steps, 'single');
    hist.force       = zeros(num.eq, num.steps, 'single');
    hist.coor        = zeros(num.np, num.ndm, num.steps, 'single');
    hist.coor(:,:,1) = nodes;
    hist.eps         = zeros(num.str,num.gp, num.el, num.steps, 'single');
    hist.stre        = zeros(num.str,num.gp, num.el, num.steps, 'single');
    hist.ctan        = zeros(3,3,3,3,num.gp, num.el, num.steps, 'single');
    
    for i=1:num.steps-1
        strgF = sprintf('%s%sforc_%0.4d', hist.start, filesep, i);
        strgd = sprintf('%s%sdisp_%0.4d', hist.start, filesep, i);
        strgc = sprintf('%s%scoor_%0.4d', hist.start, filesep, i);
        strge = sprintf('%s%sepsi_%0.4d', hist.start, filesep, i);
        strgs = sprintf('%s%sstre_%0.4d', hist.start, filesep, i);
        strgC = sprintf('%s%sctan_%0.4d', hist.start, filesep, i);
    
        fidF = fopen(strgF,'r');
        fidd = fopen(strgd,'r');
        fidc = fopen(strgc,'r');
        fide = fopen(strge,'r');
        fids = fopen(strgs,'r');
        fidC = fopen(strgC,'r');
     
        hist.force (:,i+1)	= fread(fidF, [num.eq, 1]       , 'single');
        hist.disp  (:,i+1)	= fread(fidd, [num.eq, 1]       , 'single');
        hist.coor  (:,:,i+1)= fread(fidc, [num.np, num.ndm] , 'single');
        hist.eps (:,:,:,i+1)      = reshape (fread(fide, [num.str, num.gp*num.el], 'single'), num.str, num.gp,num.el);
        hist.stre(:,:,:,i+1)      = reshape (fread(fids, [num.str, num.gp*num.el], 'single'), num.str, num.gp,num.el);
        hist.ctan(:,:,:,:,:,:,i+1)= reshape (fread(fidC, [3*3*3*3, num.gp*num.el], 'single'), 3,3,3,3, num.gp,num.el);      
        
        fclose('all');    
    end
    hist.disp = permute(reshape(hist.disp ,num.ndof,num.np,num.steps),[2 1 3]);
    hist.force= permute(reshape(hist.force,num.ndof,num.np,num.steps),[2 1 3]); 
    varargout{1} = hist;
    varargout{2} = num;
    
end
fclose('all');

end

