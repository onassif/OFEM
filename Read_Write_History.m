function [] = Read_Write_History(hist, num, step, Fext, d, coor, rw_flag)
%WriteHistory Writes history arrays to binary files

strgF = sprintf('%s/forc_%0.4d', hist.start,step);
strgd = sprintf('%s/disp_%0.4d', hist.start,step);
strgc = sprintf('%s/coor_%0.4d', hist.start,step);
strgr = sprintf('%s/resi_%0.4d', hist.start,step);
strge = sprintf('%s/epsi_%0.4d', hist.start,step);
strgs = sprintf('%s/stre_%0.4d', hist.start,step);
strgC = sprintf('%s/ctan_%0.4d', hist.start,step);

if(rw_flag)
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
else
    
    
    
%     hist.eps(:,igp,iel) = gp.eps;
%     hist.stre(:,igp,iel) = gp.sigma;
%     hist.ctan(:,:,:,:,igp,iel) = gp.ctan;
end
fclose('all');
end

