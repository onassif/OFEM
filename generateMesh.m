function [nde, el, nen, ngp, numnp, numel, ndm]= ...
   generateMesh(elmtype,coor, plot, numelx, numely, numelz)

switch elmtype
   case  'Q4'
      %%
      numnpx = numelx+1;
      numnpy = numely+1;
      numnp = numnpx*numnpy;
      numel = numelx * numely;
      nen = 4;
      ngp = 4;
      ndm = 2;
      if (size(coor) == [4 2])
      elseif (size(coor) == [2 4])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx;
      incy = (coor(4,2) - coor(2,2))/numely;
      nde = zeros(numnp,2);
      el = ones(numel,nen+1);
      
      ydir=coor(2,2);
      for j=1:numnpy
         xdir = coor(1,1);
         for i=1:numnpx
            nde(i+(j-1)*numnpx,:) = [xdir ydir];
            xdir = xdir+incx;
         end
         ydir = ydir+incy;
      end
      
      sz = [numnpx numnpy];
      for j=1:numely
         for i=1:numelx
            el(i + numelx*(j-1), 1:nen) = sub2ind(sz, [i i+1 i+1 i], [j j j+1 j+1]);
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 'Q9'
      %%
      numnpx = 2*numelx+1;
      numnpy = 2*numely+1;
      numnp = numnpx*numnpy;
      numel = numelx * numely;
      nen = 9;
      ngp = 9;
      ndm = 2;
      if (size(coor) == [4 2])
      elseif (size(coor) == [2 4])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx/2;
      incy = (coor(4,2) - coor(2,2))/numely/2;
      nde = zeros(numnp,2);
      el = ones(numel,nen+1,'uint32');
      
      ydir = coor(2,2);
      for j=1:numnpy
         xdir = coor(1,1);
         for i=1:numnpx
            nde(i+(j-1)*numnpx,:) = [xdir ydir];
            xdir = xdir+incx;
         end
         ydir = ydir+incy;
      end
      
      sz = [numnpx numnpy];
      for j=1:numely
         for i=1:numelx
            L = 2*i-1; R = 2*i+1; Mi = 2*i;
            B = 2*j-1; T = 2*j+1; Mj = 2*j;
            
            el(i + numelx*(j-1), 1:nen)=...
               sub2ind(sz, [L R R L Mi R Mi L Mi], [B B T T B Mj T Mj Mj]);
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 'T3'
      %%
      numnpx = numelx+1;  
      numnpy = numely+1;  
      numnp = numnpx*numnpy;
      nen = 3;    
      ngp = 1;    
      ndm = 2;
      numel = numelx * numely*2;
      
      if (size(coor) == [4 2])
      elseif (size(coor) == [2 4])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx;
      incy = (coor(4,2) - coor(2,2))/numely;
      nde = zeros(numnp,2);
      el = ones(numel,nen+1,'uint32');
      
      ydir=coor(2,2);
      for j=1:numnpy
         xdir = coor(1,1);
         for i=1:numnpx
            nde(i+(j-1)*numnpx,:) = [xdir ydir];
            xdir = xdir+incx;
         end
         ydir = ydir+incy;
      end
      
      sz = [numnpx numnpy];
      for j=1:numely
         for i=1:numelx
            el([2*i-1 2*i] + 2*numelx*(j-1), 1:nen) = [...
               sub2ind(sz, [i   i+1 i  ], [j   j   j+1])
               sub2ind(sz, [i+1 i   i+1], [j+1 j+1 j  ])];
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 'T6'
      %%
      numnpx = 2*numelx+1;
      numnpy = 2*numely+1;
      numnpz = 2*numelz+1;
      numnp = numnpx*numnpy;
      numel = numelx * numely * 2;
      nen = 6;
      ngp = 3;
      ndm = 2;
      if (size(coor) == [4 2])
      elseif (size(coor) == [2 4])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx/2;
      incy = (coor(4,2) - coor(4,1))/numely/2;
      nde = zeros(numnp,2);
      el = ones(numel,nen+1,'uint32');
      
      ydir=0;
      for j=1:numnpy
         xdir = 0;
         for i=1:numnpx
            nde(i+(j-1)*numnpx,:) = [xdir ydir];
            xdir = xdir+incx;
         end
         ydir = ydir+incy;
      end
      
      sz = [numnpx numnpy];
      for j=1:numely
         for i=1:numelx
            L = 2*i-1; R = 2*i+1; Mi = 2*i;
            B = 2*j-1; T = 2*j+1; Mj = 2*j;
            
            el([2*i-1 2*i] + 2*numelx*(j-1), 1:nen) = [...
               sub2ind(sz, [L R L Mi Mi L], [B B T B Mj Mj])
               sub2ind(sz, [R L R Mi Mi R], [T T B T Mj Mj])];
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case {'T4'}
      %%      
      numnpx = numelx+1;
      numnpy = numely+1;
      numnpz = numelz+1;
      numnp = numnpx * numnpy * numnpz;
      numel = numelx * numely * numelz*6;
      nen = 4;
      ngp = 1;
      ndm = 3;
      if (size(coor) == [8 3])
      elseif (size(coor) == [3 8])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx;
      incy = (coor(4,2) - coor(2,2))/numely;
      incz = (coor(6,3) - coor(3,3))/numelz;
      
      nde = zeros(numnp, ndm);
      el0 = ones(numel/6,9);
      el  = ones(numel,nen+1);
      
      zdir = coor(3,3);
      for k=1:numnpz
         ydir = coor(2,2);
         for j=1:numnpy
            xdir = coor(1,1);
            for i=1:numnpx
               nde(i+(j-1)*numnpx+(k-1)*numnpx*numnpy, :) = [xdir ydir zdir];
               xdir = xdir+incx;
            end
            ydir = ydir+incy;
         end
         zdir = zdir+incz;
      end
      
      sn = [numnpx numnpy numnpz];
      se = [numelx numely numelz];
      for k=1:numelz
         for j=1:numely
            for i=1:numelx
               el0(i + se(1)*(j-1) + se(1)*se(2)*(k-1), 1:8)=[...
                  sub2ind(sn, [i i+1 i+1 i], [j j j+1 j+1], [k   k   k   k  ]),...
                  sub2ind(sn, [i i+1 i+1 i], [j j j+1 j+1], [k+1 k+1 k+1 k+1])];
            end
         end
      end
      
      for i = 1:numel/6
         el((i-1)*6+(1:6),:) = [...
            el0(i,1) el0(i,2) el0(i,4) el0(i,8) el0(i,9)
            el0(i,3) el0(i,4) el0(i,2) el0(i,8) el0(i,9)
            el0(i,1) el0(i,5) el0(i,2) el0(i,8) el0(i,9)
            el0(i,6) el0(i,2) el0(i,5) el0(i,8) el0(i,9)
            el0(i,3) el0(i,2) el0(i,7) el0(i,8) el0(i,9)
            el0(i,6) el0(i,7) el0(i,2) el0(i,8) el0(i,9)];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case {'T10'}
      %%
      numnpx = 2*numelx+1;
      numnpy = 2*numely+1;
      numnpz = 2*numelz+1;
      
      numnp = numnpx * numnpy * numnpz;
      numel = numelx * numely * numelz*6;
      nen = 10;
      ngp = 4;
      ndm = 3;
      if (size(coor) == [8 3])
      elseif (size(coor) == [3 8])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx/2;
      incy = (coor(4,2) - coor(2,2))/numely/2;
      incz = (coor(6,3) - coor(3,3))/numelz/2;
      
      nde = zeros(numnp, ndm);
      el0 = ones(numel/6,28);
      el  = ones(numel,nen+1);
      
      zdir = coor(3,3);
      for k=1:numnpz
         ydir = coor(2,2);
         for j=1:numnpy
            xdir = coor(1,1);
            for i=1:numnpx
               nde(i+(j-1)*numnpx+(k-1)*numnpx*numnpy, :) = [xdir ydir zdir];
               xdir = xdir+incx;
            end
            ydir = ydir+incy;
         end
         zdir = zdir+incz;
      end
      
      sn = [numnpx numnpy numnpz];
      se = [numelx numely numelz];
      for k=1:numelz
         for j=1:numely
            for i=1:numelx
            Lx = 2*i-1; Ux = 2*i+1; Mx = 2*i;
            Ly = 2*j-1; Uy = 2*j+1; My = 2*j;
            Lz = 2*k-1; Uz = 2*k+1; Mz = 2*k;
               el0(i + se(1)*(j-1) + se(1)*se(2)*(k-1), 1:27)=[...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Lz Lz Lz Lz Lz Lz Lz Lz Lz]),...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Mz Mz Mz Mz Mz Mz Mz Mz Mz]),...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Uz Uz Uz Uz Uz Uz Uz Uz Uz])];
            end
         end
      end
      
      for i = 1:numel/6
         el((i-1)*6+(1:6),:) = [...
            el0(i, 1) el0(i, 7) el0(i, 3) el0(i,25) el0(i, 4) el0(i, 5) el0(i, 2) el0(i,13) el0(i,16) el0(i,14) el0(i,28)
            el0(i, 9) el0(i, 3) el0(i, 7) el0(i,25) el0(i, 6) el0(i, 5) el0(i, 8) el0(i,17) el0(i,14) el0(i,16) el0(i,28)
            el0(i, 1) el0(i, 3) el0(i,19) el0(i,25) el0(i, 2) el0(i,11) el0(i,10) el0(i,13) el0(i,14) el0(i,22) el0(i,28)
            el0(i,21) el0(i,19) el0(i, 3) el0(i,25) el0(i,20) el0(i,11) el0(i,12) el0(i,23) el0(i,22) el0(i,14) el0(i,28)
            el0(i, 9) el0(i,27) el0(i, 3) el0(i,25) el0(i,18) el0(i,15) el0(i, 6) el0(i,17) el0(i,26) el0(i,14) el0(i,28)
            el0(i,21) el0(i, 3) el0(i,27) el0(i,25) el0(i,12) el0(i,15) el0(i,24) el0(i,23) el0(i,14) el0(i,26) el0(i,28)];
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 'Q8'
      %%
      numnpx = numelx+1;
      numnpy = numely+1;
      numnpz = numelz+1;
      
      numnp = numnpx*numnpy*numnpz;
      numel = numelx * numely * numelz;
      nen = 8;
      ngp = 8;
      ndm = 3;
      if (size(coor) == [8 3])
      elseif (size(coor) == [3 8])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx;
      incy = (coor(4,2) - coor(2,2))/numely;
      incz = (coor(6,3) - coor(3,3))/numelz;
      
      nde = zeros(numnp, ndm);
      el = ones(numel,nen+1);
      
      zdir = coor(3,3);
      for k=1:numnpz
         ydir = coor(2,2);
         for j=1:numnpy
            xdir = coor(1,1);
            for i=1:numnpx
               nde(i+(j-1)*numnpx+(k-1)*numnpx*numnpy, :) = [xdir ydir zdir];
               xdir = xdir+incx;
            end
            ydir = ydir+incy;
         end
         zdir = zdir+incz;
      end
      
      sn = [numnpx numnpy numnpz];
      se = [numelx numely numelz];
      for k=1:numelz
         for j=1:numely
            for i=1:numelx
               el(i + se(1)*(j-1) + se(1)*se(2)*(k-1), 1:nen)=[...
                  sub2ind(sn, [i i+1 i+1 i], [j j j+1 j+1], [k   k   k   k  ]),...
                  sub2ind(sn, [i i+1 i+1 i], [j j j+1 j+1], [k+1 k+1 k+1 k+1])];
            end
         end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case {'Q27'}
      %%
      numnpx = 2*numelx+1;
      numnpy = 2*numely+1;
      numnpz = 2*numelz+1;
      
      numnp = numnpx*numnpy*numnpz;
      numel = numelx * numely * numelz;
      nen = 27;
      ngp = 27;
      ndm = 3;
      if (size(coor) == [8 3])
      elseif (size(coor) == [3 8])
         coor = coor';
      else
         error('wrong coordinates');
      end
      
      incx = (coor(2,1) - coor(1,1))/numelx/2;
      incy = (coor(4,2) - coor(2,2))/numely/2;
      incz = (coor(6,3) - coor(3,3))/numelz/2;
      
      nde = zeros(numnp, ndm);
      el = ones(numel,nen+1);
      
      zdir = coor(3,3);
      for k=1:numnpz
         ydir = coor(2,2);
         for j=1:numnpy
            xdir = coor(1,1);
            for i=1:numnpx
               nde(i+(j-1)*numnpx+(k-1)*numnpx*numnpy, :) = [xdir ydir zdir];
               xdir = xdir+incx;
            end
            ydir = ydir+incy;
         end
         zdir = zdir+incz;
      end
      
      sn = [numnpx numnpy numnpz];
      se = [numelx numely numelz];
      for k=1:numelz
         for j=1:numely
            for i=1:numelx
            Lx = 2*i-1; Ux = 2*i+1; Mx = 2*i;
            Ly = 2*j-1; Uy = 2*j+1; My = 2*j;
            Lz = 2*k-1; Uz = 2*k+1; Mz = 2*k;
               el(i + se(1)*(j-1) + se(1)*se(2)*(k-1), 1:nen)=[...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Lz Lz Lz Lz Lz Lz Lz Lz Lz]),...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Mz Mz Mz Mz Mz Mz Mz Mz Mz]),...
                  sub2ind(sn, [Lx Mx Ux Lx Mx Ux Lx Mx Ux], [Ly Ly Ly My My My Uy Uy Uy], [Uz Uz Uz Uz Uz Uz Uz Uz Uz])];
            end
         end
      end
end
%%%%%%%%%%% start plotting %%%%%%%%%%%
if plot
   PlotMesh(elmtype, nde, el, nen, numel)
end
end