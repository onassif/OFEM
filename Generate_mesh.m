function [nde, el, nen, ngp, numnp, numel, ndm, new_BC, new_FORCE]= ...
    Generate_mesh2(elmtype,coor, old_BC, old_FORCE, plot, numelx, numely, varargin)

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
    incy = (coor(4,2) - coor(4,1))/numely;
    nde = zeros(numnp,2);
    el = ones(numel,nen+1);
    
    ydir=0;
    for j=1:numnpy
        xdir = 0;
        for i=1:numnpx
            nde(i+(j-1)*numnpx,:) = [xdir ydir];
            xdir = xdir+incx;
        end
        ydir = ydir+incy;
    end
    
    l = 0;
    sz = [numnpx numnpy];
    for j=1:numely
        for i=1:numelx
            l = l+1;
            el(l,1:nen)=[sub2ind(sz,i,j),sub2ind(sz,i+1,j),sub2ind(sz,i+1,j+1),sub2ind(sz,i,j+1)];
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
    
    l = 0;
    sz = [numnpx numnpy];
    for j=1:numely
        for i=1:numelx
            l = l+1;
            i1 = 2*i-1; i2 = 2*i+1; i3 = 2*i;
            j1 = 2*j-1; j2 = 2*j+1; j3 = 2*j;
            el(l,1:nen)=[sub2ind(sz, i1, j1),sub2ind(sz, i2, j1),sub2ind(sz, i2, j2),sub2ind(sz, i1, j2),...
                sub2ind(sz, i3, j1),sub2ind(sz, i2, j3),sub2ind(sz, i3, j2),sub2ind(sz, i1, j3),...
                sub2ind(sz, i3, j3)];
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
    numnpx = numelx+1;
    numnpy = numely+1;
    numnp = numnpx*numnpy;
    numel = numelx * numely*2;
    nen = 3;
    ngp = 1;
    ndm = 2;
    if (size(coor) == [4 2])
    elseif (size(coor) == [2 4])
        coor = coor';
    else
        error('wrong coordinates');
    end
    
    incx = (coor(2,1) - coor(1,1))/numelx;
    incy = (coor(4,2) - coor(4,1))/numely;
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
    
    ydir=0;
    for j=1:numely
        p1 = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999, numelx);
        p2 = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999, numelx, 'last');
        p3 = find( nde(:,2)<=(ydir+incy)*1.00001 & nde(:,2)>=(ydir+incy)*0.99999, numelx);
        p4 = find( nde(:,2)<=(ydir+incy)*1.00001 & nde(:,2)>=(ydir+incy)*0.99999, numelx, 'last');
        p5 = p3;
        p6 = p2;
        
        
        el([1:2:numelx*2]+(j-1)*numelx*2, 1:nen) = [p1 p2 p3];
        el([2:2:numelx*2]+(j-1)*numelx*2, 1:nen) = [p4 p5 p6];
        ydir = ydir + incy;
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
    
    ydir=0;
    for j=1:numely
        p1 = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999);
        p1 = p1(1:2:end-1);
        p2 = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999);
        p2 = p2(2:2:end);
        p3 = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999);
        p3 = p3(3:2:end);
        p4 = find( nde(:,2)<=(ydir+incy)*1.00001 & nde(:,2)>=(ydir+incy)*0.99999);
        p4 = p4(2:2:end);
        p5 = find( nde(:,2)<=(ydir+2*incy)*1.00001 & nde(:,2)>=(ydir+2*incy)*0.99999);
        p5 = p5(1:2:end-1);
        p6 = find( nde(:,2)<=(ydir+incy)*1.00001 & nde(:,2)>=(ydir+incy)*0.99999);
        p6 = p6(1:2:end-1);
        p7 = find( nde(:,2)<=(ydir+2*incy)*1.00001 & nde(:,2)>=(ydir+2*incy)*0.99999);
        p7 = p7(3:2:end);
        p8 = find( nde(:,2)<=(ydir+2*incy)*1.00001 & nde(:,2)>=(ydir+2*incy)*0.99999);
        p8 = p8(2:2:end);
        p9 = p5;
        p10= p4;
        p11= p3;
        p12= find( nde(:,2)<=(ydir+incy)*1.00001 & nde(:,2)>=(ydir+incy)*0.99999);
        p12= p12(3:2:end);
        
        el([1:2:numelx*2]+(j-1)*numelx*2, 1:nen) = [p1  p3  p5   p2  p4  p6];
        el([2:2:numelx*2]+(j-1)*numelx*2, 1:nen) = [p7  p9  p11  p8  p10 p12];
        
        ydir = ydir + 2*incy;
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
    incy = (coor(4,2) - coor(1,2))/numely;
    incz = (coor(6,3) - coor(1,3))/numelz;
    
    nde = zeros(numnp, ndm);
    el = ones(numel,nen+1);
    
    zdir = 0;
    for k=1:numnpz
        ydir = 0;
        for j=1:numnpy
            xdir = 0;
            for i=1:numnpx
                nde(i+(j-1)*numnpx+(k-1)*numnpx*numnpy, :) = [xdir ydir zdir];
                xdir = xdir+incx;
            end
            ydir = ydir+incy;
        end
        zdir = zdir+incz;
    end
    
    l = 0;
    sz = [numnpx numnpy numnpz];
    for k=1:numelz
        for j=1:numely
            for i=1:numelx
                l = l+1;
                el(l,1:nen)=[sub2ind(sz,i,j,k),  sub2ind(sz,i+1,j,k),  sub2ind(sz,i+1,j+1,k),  sub2ind(sz,i,j+1,k),...
                    sub2ind(sz,i,j,k+1) sub2ind(sz,i+1,j,k+1) sub2ind(sz,i+1,j+1,k+1) sub2ind(sz,i,j+1,k+1)];
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
if ndm == 2
    % BC
    new_BC = 0;
    for i=1:size(old_BC,1)
        switch  old_BC{i,1}
            case 'x'
                xdir = old_BC{i,2};
                affected_nodes = find( nde(:,1)<=xdir*1.00001 & nde(:,1)>=xdir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts:ends, 3) = 1;
                    case 'v'
                        new_BC(starts:ends, 3) = 2;
                    case 'w'
                        new_BC(starts:ends, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts:ends, 4) = 0;
                
            case 'y'
                ydir = old_BC{i,2};
                affected_nodes = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts:ends, 3) = 1;
                    case 'v'
                        new_BC(starts:ends, 3) = 2;
                    case 'w'
                        new_BC(starts:ends, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts:ends, 4) = 0;
                
            case 'z'
                zdir = old_BC{i,2};
                affected_nodes = find( nde(:,3)<=zdir*1.00001 & nde(:,3)>=zdir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts:ends, 3) = 1;
                    case 'v'
                        new_BC(starts:ends, 3) = 2;
                    case 'w'
                        new_BC(starts:ends, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts:ends, 4) = 0;
                
            case 'node'
                affected_node = old_BC{i,2};
                starts = new_BC(end,1)+1;
                new_BC(starts, 1) = starts;
                new_BC(starts, 2) = affected_node;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts, 3) = 1;
                    case 'v'
                        new_BC(starts, 3) = 2;
                    case 'w'
                        new_BC(starts, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts, 4) = 0;
            otherwise
                error('unrecognized input!!');
        end
    end
    new_BC = new_BC(:,2:4);
    
    % FORCE
    new_FORCE = 0;
    for i=1:size(old_FORCE,1)
        switch  old_FORCE{i,1}
            case 'x'
                xdir = old_FORCE{i,2};
                affected_nodes = find( nde(:,1)<=xdir*1.00001 &...
                    nde(:,1)>=xdir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'y'
                ydir = old_FORCE{i,2};
                affected_nodes = find( nde(:,2)<=ydir*1.00001 &...
                    nde(:,2)>=ydir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'z'
                zdir = old_FORCE{i,2};
                affected_nodes = find( nde(:,3)<=zdir*1.00001 &...
                    nde(:,3)>=zdir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'node'
                affected_node = old_FORCE{i,2};
                starts = new_FORCE(end,1)+1;
                new_FORCE(starts, 1) = starts;
                new_FORCE(starts, 2) = affected_node;
                switch old_BC{i,3}
                    case  'u'
                        new_FORCE(starts, 3) = 1;
                    case 'v'
                        new_FORCE(starts, 3) = 2;
                    case 'w'
                        new_FORCE(starts, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts, 4) = 0;
            otherwise
                error('unrecognized input!!');
        end
    end
    new_FORCE = new_FORCE(:,2:4);
    
elseif ndm == 3
    
     % BC
    new_BC = 0;
    for i=1:size(old_BC,1)
        switch  old_BC{i,1}
            case 'x'
                xdir = old_BC{i,2};
                affected_nodes = find( nde(:,1)<=xdir*1.00001 & nde(:,1)>=xdir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
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
                new_BC(starts:ends, 4) = 0;
                
            case 'y'
                ydir = old_BC{i,2};
                affected_nodes = find( nde(:,2)<=ydir*1.00001 & nde(:,2)>=ydir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts:ends, 3) = 1;
                    case 'v'
                        new_BC(starts:ends, 3) = 2;
                    case 'w'
                        new_BC(starts:ends, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts:ends, 4) = 0;
                
            case 'z'
                zdir = old_BC{i,2};
                affected_nodes = find( nde(:,3)<=zdir*1.00001 & nde(:,3)>=zdir*0.99999);
                starts = new_BC(end,1)+1;
                ends   = new_BC(end,1)+length(affected_nodes);
                new_BC(starts:ends, 1) = starts:ends;
                new_BC(starts:ends, 2) = affected_nodes;
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
                new_BC(starts:ends, 4) = 0;
                
            case 'node'
                affected_node = old_BC{i,2};
                starts = new_BC(end,1)+1;
                new_BC(starts, 1) = starts;
                new_BC(starts, 2) = affected_node;
                switch old_BC{i,3}
                    case  'u'
                        new_BC(starts, 3) = 1;
                    case 'v'
                        new_BC(starts, 3) = 2;
                    case 'w'
                        new_BC(starts, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts, 4) = 0;
            otherwise
                error('unrecognized input!!');
        end
    end
    new_BC = new_BC(:,2:4);
    
    % FORCE
    new_FORCE = 0;
    for i=1:size(old_FORCE,1)
        switch  old_FORCE{i,1}
            case 'x'
                xdir = old_FORCE{i,2};
                affected_nodes = find( nde(:,1)<=xdir*1.00001 & nde(:,1)>=xdir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'y'
                ydir = old_FORCE{i,2};
                affected_nodes = find( nde(:,2)<=ydir*1.00001 &...
                    nde(:,2)>=ydir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'z'
                zdir = old_FORCE{i,2};
                affected_nodes = find( nde(:,3)<=zdir*1.00001 &...
                    nde(:,3)>=zdir*0.99999);
                portions = 0;
                for j=1:length(affected_nodes)
                    portions = portions + sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                portioned_load = old_FORCE{i,4}/portions;
                starts = new_FORCE(end,1);
                
                for j=1:length(affected_nodes)
                    new_FORCE(starts+j, 1) = starts+j;
                    new_FORCE(starts+j, 2) = affected_nodes(j);
                    new_FORCE(starts+j, 4) = portioned_load *...
                        sum(length(find(el(:,1:nen) == affected_nodes(j))));
                end
                
                switch old_FORCE{i,3}
                    case  'u'
                        new_FORCE(starts+1:starts+j, 3) = 1;
                    case 'v'
                        new_FORCE(starts+1:starts+j, 3) = 2;
                    case 'w'
                        new_FORCE(starts+1:starts+j, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
            case 'node'
                affected_node = old_FORCE{i,2};
                starts = new_FORCE(end,1)+1;
                new_FORCE(starts, 1) = starts;
                new_FORCE(starts, 2) = affected_node;
                switch old_BC{i,3}
                    case  'u'
                        new_FORCE(starts, 3) = 1;
                    case 'v'
                        new_FORCE(starts, 3) = 2;
                    case 'w'
                        new_FORCE(starts, 3) = 3;
                    otherwise
                        error('unrecognized input!!');
                end
                new_BC(starts, 4) = 0;
            otherwise
                error('unrecognized input!!');
        end
    end
    new_FORCE = new_FORCE(:,2:4);   
    

end
end
