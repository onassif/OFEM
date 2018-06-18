function [new_BC, new_FORCE] = generateBC(old_BC, old_FORCE, nde, el, elmtype, ndm, numnp, numel)
%% Generate Faces 
fc =Faces(elmtype, el, numel);

%% BC
new_BC = 0;
for i=1:size(old_BC,1)
   switch  old_BC{i,1}
      case 'x'
         dir = 1;
         batch = true;
      case 'y'
         dir = 2;
         batch = true;
      case 'z'
         dir = 3;
         batch = true;
      case 'node'
         affected_node = old_BC{i,2};
         starts = new_BC(end,1)+1;
         ends   = new_BC(end,1)+1;
         new_BC(starts, 1) = starts;
         new_BC(starts, 2) = affected_node;
         batch = false;
      otherwise
         error('unrecognized input!!');
   end
   if batch
      
      lmt = old_BC{i,2};
      if sign(lmt) >=0
         affected_nodes = find( nde(:,dir)<=lmt*1.001 & nde(:,dir)>=lmt*0.999 );
      else
         affected_nodes = find( nde(:,dir)>=lmt*1.001 & nde(:,dir)<=lmt*0.999 );
      end
      
%       nBC = [affected_nodes , dir*ones(size(affected_nodes)), zeros(size(affected_nodes))];
%       fc.nodes = affected_nodes;
%       for ifc = 1:size(fc.faces, 1)
%          fc.ifc = ifc;
%          fc.coor = nde(fc.faces(ifc,:), (1:ndm ~= dir));
%          for igp = 1:fc.numGP
%             fc.gp.i = igp;
%             if igp == 1 % No need to compute J for each gp
%                fc.gp.dXdxi = (fc.gp.dNdxi*fc.coor)';
%                J = det(fc.gp.dXdxi);
%             end
%             nBC(fc.indices, 3) = nBC(fc.indices, 3) + fc.gp.N'*old_BC{i,4}*J*fc.gp.w;
%          end
%       end
      starts = new_BC(end,1)+1;
      ends   = new_BC(end,1)+length(affected_nodes);
      new_BC(starts:ends, 1) = starts:ends;
      new_BC(starts:ends, 2) = affected_nodes;
   end
   
   switch old_BC{i,3}
      case 'u'
         new_BC(starts:ends, 3) = 1;
      case 'v'
         new_BC(starts:ends, 3) = 2;
      case 'w'
         new_BC(starts:ends, 3) = 3;
      otherwise
         error('unrecognized input!!');
   end
   new_BC(starts:ends, 4) = old_BC{i,4};
   new_BC(starts:ends, 5) = i;
   
end
new_BC = new_BC(:,2:5);
rowsToDelete = [];
for i = 1:numnp
   targetNode = find(new_BC(:,1)==i);
   if length(targetNode)>1
      for j = 1:ndm
         sameDir = find(new_BC(targetNode,2)==j);
         targetHomoBC    = find(new_BC(targetNode(sameDir),3) == 0);
         targetNoHomoBC  = find(new_BC(targetNode(sameDir),3) ~= 0);
         if isempty(targetHomoBC) && length(targetNoHomoBC)>1
            error("One node, one direction but two non-homogeneous BCs");
         elseif length(targetHomoBC)>=1 && length(targetNoHomoBC)>=1
            for k=1:length(targetNoHomoBC)
               rowsToDelete = [rowsToDelete;targetNode(sameDir(targetNoHomoBC))];
            end
         end
      end
   end
end
rowsToKeep = logical(1:size(new_BC,1))';
rowsToKeep(rowsToDelete) = 0;
new_BC = new_BC(rowsToKeep,:);

%% FORCE
new_FORCE = [];
for i=1:size(old_FORCE,1)
   switch  old_FORCE{i,1}
      case 'x'
         dir = 1;
         batch = true;
      case 'y'
         dir = 2;
         batch = true;
      case 'z'
         dir = 3;
         batch = true;
      case 'node'
         if isempty(new_FORCE)
            new_FORCE = 0;
         end
         affected_node = old_FORCE{i,2};
         starts = new_FORCE(end,1)+1;
         ends   = new_FORCE(end,1)+1;
         new_FORCE(starts, 1) = starts;
         new_FORCE(starts, 2) = affected_node;
         batch = false;
      otherwise
         error('unrecognized input!!');
   end
   if batch
      lmt = old_FORCE{i,2};
      if sign(lmt) >=0
         affected_nodes = find( nde(:,dir)<=lmt*1.001 & nde(:,dir)>=lmt*0.999 );
      else
         affected_nodes = find( nde(:,dir)>=lmt*1.001 & nde(:,dir)<=lmt*0.999 );
      end
      
      lnth = size(affected_nodes);
      nF = [affected_nodes , dir*ones(lnth), zeros(lnth), i*ones(lnth)];
      fc.nodes = affected_nodes;
      for ifc = 1:size(fc.faces, 1)
         fc.ifc = ifc;
         fc.gp.mesh = struct('nodes', nde(:,std(nde(fc.faces(1,:),:))>1e-8), 'conn', fc.faces);
         for igp = 1:fc.numGP
            fc.gp.i = igp; fc.gp.iel = ifc;
            
            nF(fc.indices, 3) = nF(fc.indices, 3) + fc.gp.N'*old_FORCE{i,4}*abs(fc.gp.J)*fc.gp.w;
         end
      end
      new_FORCE = [new_FORCE; nF];
   else
   
   switch old_FORCE{i,3}
      case 'u'
         new_FORCE(starts+1:starts+j, 3) = 1;
      case 'v'
         new_FORCE(starts+1:starts+j, 3) = 2;
      case 'w'
         new_FORCE(starts+1:starts+j, 3) = 3;
      otherwise
         error('unrecognized input!!');
   end
   
   new_FORCE(starts+1:starts+j, 5) = i;
   end 
end
if ~batch
   new_FORCE = new_FORCE(:,2:5);
end
end

