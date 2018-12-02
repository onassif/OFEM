function [] = PlotMesh(elmtype, nde, el, nen, numel, color, ax)
if ~exist('ax','var')
   figure('Position',[100,50,1200,700]);
   ax = gca;
end
if ~exist('color','var')
   color = 'r';
end
switch elmtype
   case 'Q4'
      patch(ax,'Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:),...
         'Faces',reshape(1:nen*numel,nen,numel)',...
         'FaceColor','none', 'EdgeColor',color, 'Marker', 'o', 'MarkerFaceColor',color);
      
   case 'Q9'
      el2 = el(:,[1 5 2 6 3 7 4 8 9]);
      patch(ax,'Vertices',nde(reshape(el2(:,1:nen-1)',numel*(nen-1),1),:),...
         'Faces',reshape(1:(nen-1)*numel,(nen-1),numel)',...
         'FaceColor','none', 'EdgeColor', color, 'Marker','o', 'MarkerFaceColor',color);

   case 'T3'
      patch(ax,'Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:),...
         'Faces',reshape(1:nen*numel,nen,numel)',...
         'FaceColor','none', 'EdgeColor',color, 'Marker','o', 'MarkerFaceColor',color);
      
   case 'T6'
      el2 = el(:,[1 4 2 5 3 6]);
      patch(ax,'Vertices',nde(reshape(el2(:,1:nen)',numel*(nen),1),:),...
         'Faces',reshape(1:(nen)*numel,(nen),numel)',...
         'FaceColor','none', 'EdgeColor',color, 'Marker','o', 'MarkerFaceColor',color);

   case {'T4'}
      i = [...
         1 2 4, 1 4 8,...
         3 4 8, 4 2 8,...
         1 5 2, 5 2 8,...
         6 2 8, 6 5 8,...
         3 2 7, 3 2 8,...
         6 7 8, 7 2 8];
      fac = reshape(1:8*numel/6,8,numel/6);
      fac = reshape(fac(i,:),3,numel*2)';
      Q8el = zeros(numel/6,8);
      for i = 1:numel/6
        Q8el(i,:) = [...
           el(((i-1)*6+1),1), el(((i-1)*6+1),2), el(((i-1)*6+2),1), el(((i-1)*6+1),3),...
           el(((i-1)*6+3),2), el(((i-1)*6+4),1), el(((i-1)*6+5),3), el(((i-1)*6+1),4)];
      end
      nodes = nde(reshape(Q8el(:,1:8)',numel/6*8,1),:);
      patch(ax,'Vertices',nodes, 'Faces',fac,...
         'FaceColor','none','EdgeColor',color, 'Marker', 'o', 'MarkerFaceColor', color);
      view (ax,[0.2 -0.4 0.1])
      set(ax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])      

   case {'Q8'}
      i = [1 2 6 5, 2 3 7 6, 3 4 8 7, 4 1 5 8, 1 2 3 4, 5 6 7 8];
      fac = reshape(1:nen*numel,nen,numel);
      fac = reshape(fac(i,:),nen/2,6*numel)';
      patch(ax,'Vertices',nde(reshape(el(:,1:nen)',numel*nen,1),:), 'Faces',fac,...
         'FaceColor','none','EdgeColor',color, 'Marker', 'o', 'MarkerFaceColor', color);
      view (ax,[0.2 -0.4 0.1])
      set(ax,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
end
axis equal
end