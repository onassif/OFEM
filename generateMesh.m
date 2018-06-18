function [nde, el, nen, ngp, numnp, numel, ndm]= ...
   generateMesh(elmtype,coor, plot, numelx, numely, varargin)

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
      numnpx = numelx+1;  numnpy = numely+1;  numnp = numnpx*numnpy;
      nen = 3;    ngp = 1;    ndm = 2;
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
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case {'Q8','Q8Crys'}
      %%
      if (isempty(varargin))
         error("specified 3D element type but didn't specify increments in z direction");
      end
      numelz = varargin{1};
      
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
end
%%%%%%%%%%% start plotting %%%%%%%%%%%
if plot
   PlotMesh(elmtype, nde, el, nen, numel)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% BC
% new_BC = 0;
% for i=1:size(old_BC,1)
%    switch  old_BC{i,1}
%       case 'x'
%          dir = 1;
%          batch = true;
%       case 'y'
%          dir = 2;
%          batch = true;
%       case 'z'
%          dir = 3;
%          batch = true;
%       case 'node'
%          affected_node = old_BC{i,2};
%          starts = new_BC(end,1)+1;
%          ends   = new_BC(end,1)+1;
%          new_BC(starts, 1) = starts;
%          new_BC(starts, 2) = affected_node;
%          batch = false;
%       otherwise
%          error('unrecognized input!!');
%    end
%    if batch
%       
%       lmt = old_BC{i,2};
%       if sign(lmt) >=0
%          affected_nodes = find( nde(:,dir)<=lmt*1.001 & nde(:,dir)>=lmt*0.999 );
%       else
%          affected_nodes = find( nde(:,dir)>=lmt*1.001 & nde(:,dir)<=lmt*0.999 );
%       end
%       
% %       nBC = [affected_nodes , dir*ones(size(affected_nodes)), zeros(size(affected_nodes))];
% %       fc.nodes = affected_nodes;
% %       for ifc = 1:size(fc.faces, 1)
% %          fc.ifc = ifc;
% %          fc.coor = nde(fc.faces(ifc,:), (1:ndm ~= dir));
% %          for igp = 1:fc.numGP
% %             fc.gp.i = igp;
% %             if igp == 1 % No need to compute J for each gp
% %                fc.gp.dXdxi = (fc.gp.dNdxi*fc.coor)';
% %                J = det(fc.gp.dXdxi);
% %             end
% %             nBC(fc.indices, 3) = nBC(fc.indices, 3) + fc.gp.N'*old_BC{i,4}*J*fc.gp.w;
% %          end
% %       end
%       starts = new_BC(end,1)+1;
%       ends   = new_BC(end,1)+length(affected_nodes);
%       new_BC(starts:ends, 1) = starts:ends;
%       new_BC(starts:ends, 2) = affected_nodes;
%    end
%    
%    switch old_BC{i,3}
%       case 'u'
%          new_BC(starts:ends, 3) = 1;
%       case 'v'
%          new_BC(starts:ends, 3) = 2;
%       case 'w'
%          new_BC(starts:ends, 3) = 3;
%       otherwise
%          error('unrecognized input!!');
%    end
%    new_BC(starts:ends, 4) = old_BC{i,4};
%    new_BC(starts:ends, 5) = i;
%    
% end
% new_BC = new_BC(:,2:5);
% rowsToDelete = [];
% for i = 1:numnp
%    targetNode = find(new_BC(:,1)==i);
%    if length(targetNode)>1
%       for j = 1:ndm
%          sameDir = find(new_BC(targetNode,2)==j);
%          targetHomoBC    = find(new_BC(targetNode(sameDir),3) == 0);
%          targetNoHomoBC  = find(new_BC(targetNode(sameDir),3) ~= 0);
%          if isempty(targetHomoBC) && length(targetNoHomoBC)>1
%             error("One node, one direction but two non-homogeneous BCs");
%          elseif length(targetHomoBC)>=1 && length(targetNoHomoBC)>=1
%             for k=1:length(targetNoHomoBC)
%                rowsToDelete = [rowsToDelete;targetNode(sameDir(targetNoHomoBC))];
%             end
%          end
%       end
%    end
% end
% rowsToKeep = logical(1:size(new_BC,1))';
% rowsToKeep(rowsToDelete) = 0;
% new_BC = new_BC(rowsToKeep,:);
% %%
% % FORCE
% new_FORCE = [];
% for i=1:size(old_FORCE,1)
%    switch  old_FORCE{i,1}
%       case 'x'
%          dir = 1;
%          batch = true;
%       case 'y'
%          dir = 2;
%          batch = true;
%       case 'z'
%          dir = 3;
%          batch = true;
%       case 'node'
%          if isempty(new_FORCE)
%             new_FORCE = 0;
%          end
%          affected_node = old_FORCE{i,2};
%          starts = new_FORCE(end,1)+1;
%          ends   = new_FORCE(end,1)+1;
%          new_FORCE(starts, 1) = starts;
%          new_FORCE(starts, 2) = affected_node;
%          batch = false;
%       otherwise
%          error('unrecognized input!!');
%    end
%    if batch
%       lmt = old_FORCE{i,2};
%       if sign(lmt) >=0
%          affected_nodes = find( nde(:,dir)<=lmt*1.001 & nde(:,dir)>=lmt*0.999 );
%       else
%          affected_nodes = find( nde(:,dir)>=lmt*1.001 & nde(:,dir)<=lmt*0.999 );
%       end
%       
%       lnth = size(affected_nodes);
%       nF = [affected_nodes , dir*ones(lnth), zeros(lnth), i*ones(lnth)];
%       fc.nodes = affected_nodes;
%       for ifc = 1:size(fc.faces, 1)
%          fc.ifc = ifc;
%          fc.gp.mesh = struct('nodes', nde(:,std(nde(fc.faces(1,:),:))>1e-8), 'conn', fc.faces);
%          for igp = 1:fc.numGP
%             fc.gp.i = igp; fc.gp.iel = ifc;
%             
%             nF(fc.indices, 3) = nF(fc.indices, 3) + fc.gp.N'*old_FORCE{i,4}*abs(fc.gp.J)*fc.gp.w;
%          end
%       end
%       new_FORCE = [new_FORCE; nF];
%    else
%    
%    switch old_FORCE{i,3}
%       case 'u'
%          new_FORCE(starts+1:starts+j, 3) = 1;
%       case 'v'
%          new_FORCE(starts+1:starts+j, 3) = 2;
%       case 'w'
%          new_FORCE(starts+1:starts+j, 3) = 3;
%       otherwise
%          error('unrecognized input!!');
%    end
%    
%    new_FORCE(starts+1:starts+j, 5) = i;
%    end 
% end
% if ~batch
%    new_FORCE = new_FORCE(:,2:5);
% end
end