function [] = PlotMesh(elmtype, nde, el, nen, numel, color, ax)
if ~exist('ax','var')
   figure('Position',[100,50,1200,700]);
   ax = gca;
end
if ~exist('color','var')
   color = 'r';
end
switch elmtype
   case {'Q4', 'T3'}
      conn = el(1:numel,1:nen)';
      X = reshape(nde(conn,1), nen, numel);
      Y = reshape(nde(conn,2), nen, numel);
      patch(ax,'XData',X, 'YData',Y, 'FaceColor','none', 'EdgeColor',color, 'Marker', 'o',...
         'MarkerFaceColor',color, 'MarkerSize',3);
      
   case 'Q9'
      conn = el(1:numel,[1 5 2 6 3 7 4 8])';
      X = reshape(nde(conn,1), nen-1, numel);
      Y = reshape(nde(conn,2), nen-1, numel);
      hold on
      patch(ax,'XData',X, 'YData',Y, 'FaceColor','none', 'EdgeColor',color, 'Marker', 'o',...
         'MarkerFaceColor',color, 'MarkerSize',3);
      scatter(nde(el(1:numel,9)',1), nde(el(1:numel,9)',2), 'Marker', 'o',...
         'MarkerFaceColor',color, 'MarkerEdgeColor',color, 'SizeData', 10);
      
   case 'T6'
      conn = el(:,[1,4,2,5,3,6])';
      X = reshape(nde(conn,1), nen, numel);
      Y = reshape(nde(conn,2), nen, numel);
      patch(ax,'XData',X, 'YData',Y, 'FaceColor','none', 'EdgeColor',color, 'Marker', 'o',...
         'MarkerFaceColor',color, 'MarkerSize',3);
      
   case 'T4'
      fac = zeros(3,numel*4);
      for i = 1:numel
         fac(:,(i-1)*4+(1:4)) = [el(i,[1,2,3]); el(i,[1,2,4]); el(i,[1,3,4]); el(i,[2,3,4])]';
      end
      X = reshape(nde(fac,1), size(fac,1), size(fac,2));
      Y = reshape(nde(fac,2), size(fac,1), size(fac,2));
      Z = reshape(nde(fac,3), size(fac,1), size(fac,2));
   
      patch(ax, 'XData',X, 'YData',Y, 'ZData',Z, 'FaceColor','none', 'EdgeColor',color,...
         'Marker','o', 'MarkerFaceColor',color, 'MarkerSize',3);
      view (gca,[0.2 -0.4 0.1])
      set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
   case 'T10'
     fac = zeros(6,numel*4);
      for i = 1:numel
         fac(:,(i-1)*4+(1:4)) = [...
            el(i,[1,7,3,6,2,5]);el(i,[1,7,3,10,4,8]);el(i,[2,9,4,8,1,5]);el(i,[2,6,3,10,4,9])]';
      end
      X = reshape(nde(fac,1), size(fac,1), size(fac,2));
      Y = reshape(nde(fac,2), size(fac,1), size(fac,2));
      Z = reshape(nde(fac,3), size(fac,1), size(fac,2));
   
      patch(ax, 'XData',X, 'YData',Y, 'ZData',Z, 'FaceColor','none', 'EdgeColor',color,...
         'Marker','o', 'MarkerFaceColor',color, 'MarkerSize',3);
      view (gca,[0.2 -0.4 0.1])
      set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
      
   case 'Q8'
      fac = zeros(4,numel*6);
      for i = 1:numel
         fac(:,(i-1)*6+(1:6)) = [...
            el(i,1:4);el(i,5:8);el(i,[1,2,6,5]);el(i,[4,3,7,8]);el(i,[1,5,8,4]);el(i,[2,6,7,3])]';
      end
      X = reshape(nde(fac,1), size(fac,1), size(fac,2));
      Y = reshape(nde(fac,2), size(fac,1), size(fac,2));
      Z = reshape(nde(fac,3), size(fac,1), size(fac,2));
   
      patch(ax, 'XData',X, 'YData',Y, 'ZData',Z, 'FaceColor','none', 'EdgeColor',color,...
         'Marker','o', 'MarkerFaceColor',color, 'MarkerSize',3);
      view (gca,[0.2 -0.4 0.1])
      set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
   case 'Q27'
      fac = zeros(8,numel*6);
      cen = zeros(7*numel,3);
      for i = 1:numel
         fac(:,(i-1)*6+(1:6)) = [...
            el(i,[1  2  3  6  9  8  7  4]); el(i,[19 20 21 24 27 26 25 22])
            el(i,[1  2  3 12 21 20 19 10]); el(i,[ 7  8  9 18 27 26 25 16])
            el(i,[3 12 21 24 27 18  9  6]); el(i,[ 1 10 19 22 25 16  7  4])]';
         cen((i-1)*7+(1:7),:) = [...
            nde(el(i,5) ,:); nde(el(i,23),:)
            nde(el(i,11),:); nde(el(i,17),:)
            nde(el(i,15),:); nde(el(i,13),:); nde(el(i,14),:)];
      end
      X = reshape(nde(fac,1), size(fac,1), size(fac,2));
      Y = reshape(nde(fac,2), size(fac,1), size(fac,2));
      Z = reshape(nde(fac,3), size(fac,1), size(fac,2));
      hold on
      patch(ax, 'XData',X, 'YData',Y, 'ZData',Z, 'FaceColor','none', 'EdgeColor',color,...
         'Marker','o', 'MarkerFaceColor',color, 'MarkerSize',3);
      scatter3(cen(:,1), cen(:,2), cen(:,3), 'Marker', 'o',...
         'MarkerFaceColor',color, 'MarkerEdgeColor',color, 'SizeData', 10);
      view (gca,[0.2 -0.4 0.1])
      set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
      
end
axis equal
end