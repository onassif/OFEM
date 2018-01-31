function [nde, el, nen, ngp, numnp, numel, ndm, new_BC, new_FORCE]= ...
    Generate_mesh(elmtype,coor, old_BC, old_FORCE, plot, numelx, numely, varargin)

if elmtype == 'Q4'
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
    
    %%%%%%%%%%% start plotting %%%%%%%%%%%
    if plot
        figure('Position',[100,50,1200,700]);
        patch(reshape(nde(el(:,1:nen)',1),nen,numel),...
            reshape(nde(el(:,1:nen)',2),nen,numel),'r',...
            'FaceColor','none','EdgeColor','r', 'Marker', 'o', 'MarkerFaceColor', 'r');
        axis equal
        hold on
        
        
        for i=1:numel
            text(mean([nde(el(i,1),1),nde(el(i,2),1)]),...
                mean([nde(el(i,2),2),nde(el(i,3),2)]),num2str(i));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Q9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif elmtype == 'Q9'
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
    
    
    if plot
        figure('Position',[100,50,1200,700]);
        el2 = el(:,[1 5 2 6 3 7 4 8 9]);
        patch('Vertices',nde(reshape(el2(:,1:nen-1)',numel*(nen-1),1),:),...
            'Faces',reshape(1:(nen-1)*numel,(nen-1),numel)',...
            'FaceColor','none','EdgeColor','r');
        axis equal
        hold on
        scatter(nde(:,1), nde(:,2), 'Marker', 'o', 'MarkerFaceColor', 'r');
        
        for i=1:numel
            text(mean([nde(el(i,1),1),nde(el(i,2),1)]),...
                mean([nde(el(i,3),2),nde(el(i,2),2)]),num2str(i));
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif elmtype == 'T3'
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
    
    if plot
        figure('Position',[100,50,1200,700]);
        patch('Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:),...
            'Faces',reshape(1:nen*numel,nen,numel)',...
            'FaceColor','none','EdgeColor','r', 'Marker', 'o', 'MarkerFaceColor', 'r');
        axis equal
        hold on
        
        for i=1:numel
            text(1/3*(nde(el(i,2),1) - nde(el(i,1),1))+ nde(el(i,1),1),...
                1/3*(nde(el(i,3),2) - nde(el(i,1),2))+ nde(el(i,1),2),num2str(i));
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif elmtype == 'T6'
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
    
    if plot
        figure('Position',[100,50,1200,700]);
        el2 = el(:,[1 4 2 5 3 6]);
        patch('Vertices',nde(reshape(el2(:,1:nen)',numel*(nen),1),:),...
            'Faces',reshape(1:(nen)*numel,(nen),numel)',...
            'FaceColor','none','EdgeColor','r');
        axis equal
        hold on
        scatter(nde(:,1), nde(:,2), 'Marker', 'o', 'MarkerFaceColor', 'r');
        
        for i=1:numel
            text(1/3*(nde(el(i,2),1) - nde(el(i,1),1))+ nde(el(i,1),1),...
                1/3*(nde(el(i,3),2) - nde(el(i,1),2))+ nde(el(i,1),2),num2str(i));
        end
        
    end
    
elseif elmtype == 'Q8'
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
    
    
    %%%%%%%%%%% start plotting %%%%%%%%%%%
    if plot
        i = [1 2 6 5 2 3 7 6 3 4 8 7 4 1 5 8 1 2 3 4 5 6 7 8];
        fac = reshape(1:nen*numel,nen,numel);
        fac = reshape(fac(i,:),nen/2,6*numel)';
        figure('Position',[100,50,1200,700]);
        patch('Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:), 'Faces',fac,...
            'FaceColor','none','EdgeColor','r', 'Marker', 'o', 'MarkerFaceColor', 'r');
        axis equal
        set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
        view (gca,[-0.8 -0.4 0.1])
        hold on
        for i=1:numel
            text(mean([nde(el(i,1),1), nde(el(i,2),1)]),...
                mean([nde(el(i,2),2), nde(el(i,3),2)]),...
                mean([nde(el(i,4),3), nde(el(i,6),3)]),...
                num2str(i));
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% BC
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
%%
% FORCE
new_FORCE = 0;
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
            starts = new_FORCE(end,1)+1;
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
        
        portions = 0;
        for j=1:length(affected_nodes)
            if elmtype == 'T3'
                portions = portions + sum(length(find(el(:,1:nen-1) == affected_nodes(j))));
            else
                portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
            end
        end
        portioned_load = old_FORCE{i,4}/portions;
        starts = new_FORCE(end,1);
        
        for j=1:length(affected_nodes)
            new_FORCE(starts+j, 1) = starts+j;
            new_FORCE(starts+j, 2) = affected_nodes(j);
            if elmtype == 'T3'
                new_FORCE(starts+j, 4) = portioned_load *...
                    sum(length(find(el(:,1:nen-1) == affected_nodes(j))));
            else
                new_FORCE(starts+j, 4) = portioned_load *...
                    sum(length(find(el(:,1:nen) == affected_nodes(j))));
            end
            
        end
    end
            
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
new_FORCE = new_FORCE(:,2:5);

end
