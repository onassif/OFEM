function [new_BC, new_FORCE] = generateBC(old_BC, old_FORCE, nde, el, elmtype, ndm, numnp, numel)
%% Generate Faces
fc =Faces(elmtype, el, numel);

%% BC
new_BC = [];
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
         starts = size(new_BC,1)+1;
         ends   = size(new_BC,1)+1;
         new_BC(starts, 1) = affected_node;
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
      
      starts = size(new_BC,1)+1;
      ends   = size(new_BC,1)+length(affected_nodes);
      new_BC(starts:ends, 1) = affected_nodes;
   end
   
   switch old_BC{i,3}
      case 'u'
         new_BC(starts:ends, 2:4) = repmat([1, old_BC{i,4}, i], ends-starts+1,1);
      case 'v'
         new_BC(starts:ends, 2:4) = repmat([2, old_BC{i,4}, i], ends-starts+1,1);
      case 'w'
         new_BC(starts:ends, 2:4) = repmat([3, old_BC{i,4}, i], ends-starts+1,1);
      otherwise
         error('unrecognized input!!');
   end  
end
% If two (or more) BC exist for the same node, remove the non-homogenous one 
rowsToDelete = false(size(new_BC,1),1);
for i = 1:numnp
   targetNode = find(new_BC(:,1)==i);
   if length(targetNode)>1
      for j = 1:ndm
         sameDir = find(new_BC(targetNode,2)==j);
         targetHomoBC    = new_BC(targetNode(sameDir),3) == 0;
         targetNoHomoBC  = new_BC(targetNode(sameDir),3) ~= 0;
         if sum(targetHomoBC)==0 && sum(targetNoHomoBC)>1
            error("One node, one direction but two non-homogeneous BCs");
         elseif sum(targetHomoBC) && sum(targetNoHomoBC)
            rowsToDelete(targetNode(sameDir(targetNoHomoBC))) = true;
         end
      end
   end
end
new_BC(rowsToDelete,:) = [];

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
         affected_node = old_FORCE{i,2};
         starts = size(new_FORCE,1)+1;
         ends   = size(new_FORCE,1)+1;
         new_FORCE(starts, :) = [affected_node,0,0,0];
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
         [fc.gp.det_dXdxi_list, fc.gp.dNdX_list, fc.gp.dXdxi_list] = ...
            shapeRef(nde(:,std(nde(fc.faces(1,:),:))>1e-8), fc.faces, fc.gp.dNdxi_list);
         for igp = 1:fc.numGP
            fc.gp.i = igp; fc.gp.iel = ifc;
            
            nF(fc.indices, 3) = nF(fc.indices, 3) + fc.gp.N*old_FORCE{i,4}*abs(fc.gp.J)*fc.gp.w;
         end
      end
      new_FORCE = [new_FORCE; nF];
   else
      switch old_FORCE{i,3}
         case 'u'
            new_FORCE(starts:ends, 2:4) = [1, old_FORCE{i,4}, i];
         case 'v'
            new_FORCE(starts:ends, 2:4) = [2, old_FORCE{i,4}, i];
         case 'w'
            new_FORCE(starts:ends, 2:4) = [3, old_FORCE{i,4}, i];
         otherwise
            error('unrecognized input!!');
      end
   end
end
end