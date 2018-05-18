function [] = PlotMesh(elmtype, nde, el, nen, numel)
        figure('Position',[100,50,1200,700]);
switch elmtype
    case 'Q4'
        patch(reshape(nde(el(:,1:nen)',1),nen,numel),...
            reshape(nde(el(:,1:nen)',2),nen,numel),'r',...
            'FaceColor','none','EdgeColor','r', 'Marker', 'o', 'MarkerFaceColor', 'r');
        axis equal
        hold on
        
        
        for i=1:numel
            text(mean([nde(el(i,1),1),nde(el(i,2),1)]),...
                mean([nde(el(i,2),2),nde(el(i,3),2)]),num2str(i));
        end
        
    case 'Q9'
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
        
    case 'T3'
        patch('Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:),...
            'Faces',reshape(1:nen*numel,nen,numel)',...
            'FaceColor','none','EdgeColor','r', 'Marker', 'o', 'MarkerFaceColor', 'r');
        axis equal
        hold on
        
        for i=1:numel
            text(1/3*(nde(el(i,2),1) - nde(el(i,1),1))+ nde(el(i,1),1),...
                1/3*(nde(el(i,3),2) - nde(el(i,1),2))+ nde(el(i,1),2),num2str(i));
        end
        
    case 'T6'
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
        
    case {'Q8','Q8Crys'}
        i = [1 2 6 5 2 3 7 6 3 4 8 7 4 1 5 8 1 2 3 4 5 6 7 8];
        fac = reshape(1:nen*numel,nen,numel);
        fac = reshape(fac(i,:),nen/2,6*numel)';
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

