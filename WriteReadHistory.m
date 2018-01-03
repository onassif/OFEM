function varargout = WriteReadHistory(hist, num, varargin )
%WriteHistory Writes & reads history arrays to binary files

if(length(varargin)==4)
    step  = varargin{1};
    Fext  = varargin{2};
    d     = varargin{3};
    coor  = varargin{4};
    
    strgF = sprintf('%s/forc_%0.4d', hist.start,step);
    strgd = sprintf('%s/disp_%0.4d', hist.start,step);
    strgc = sprintf('%s/coor_%0.4d', hist.start,step);
    strgr = sprintf('%s/resi_%0.4d', hist.start,step);
    strge = sprintf('%s/epsi_%0.4d', hist.start,step);
    strgs = sprintf('%s/stre_%0.4d', hist.start,step);
    strgC = sprintf('%s/ctan_%0.4d', hist.start,step);

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
    
elseif(length(varargin)==1)
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
        strgF = sprintf('%s/forc_%0.4d', hist.start,i);
        strgd = sprintf('%s/disp_%0.4d', hist.start,i);
        strgc = sprintf('%s/coor_%0.4d', hist.start,i);
        strge = sprintf('%s/epsi_%0.4d', hist.start,i);
        strgs = sprintf('%s/stre_%0.4d', hist.start,i);
        strgC = sprintf('%s/ctan_%0.4d', hist.start,i);
    
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
    varargout{1} = hist;
    varargout{2} = num;
    
end
fclose('all');

end

