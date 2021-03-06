function PostProcess(hist, num, gp)
hist.strss.voigt    = hist.stre;
hist.strss.mat      = zeros(3, 3, num.gp, num.el, num.steps); % matrix form
hist.strss.mat_n    = zeros(3, 3, num.np, num.steps); % nodal matrix form
hist.strss.P        = zeros(3,    num.np, num.steps);  % Principle stress
hist.strss.vm       = zeros(num.np, num.steps); % Von mises
hist.strss.avg      = zeros(6, num.el, num.steps); % Average stress (voigt)

hist        = Compute_mat(hist, num);
hist        = Compute_mat_nodal(hist, num, gp.Ninv);
hist        = Compute_principal(hist, num);
hist        = Compute_von(hist);
hist        = Compute_avg(hist, num);

mainH = figure;
mainH.Position = [25,50,1500,700];
subplot(2,2,1)
plotFvD(hist,num.steps)

subH = subplot(2,2,2);
plotcont(hist, num, mainH, subH)

ax = subplot(2,2,3);
plotDeformation(hist,num, ax)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hist = Compute_mat(hist,num)
% Converts voigt notation stress into 3x3 matrix:
ngp = num.gp; nel = num.el; nsteps = num.steps;

if num.ndm == 2
   ctan = reshape(hist.ctan(3,3,1:2,1:2,:,:,:), 4, ngp, nel, nsteps);
   
   hist.strss.mat(1:2,1:2,:,:,:) = reshape(hist.stre([1,3,3,2],:,:,:), 2,2,ngp,nel,nsteps);
   hist.strss.mat(3,3,:,:,:) = sum(ctan([1,4,2],:,:,:) .* hist.eps);
   
elseif num.ndm == 3
   i = [1,4,6,4,2,5,6,5,3];
   hist.strss.mat = reshape(hist.stre(i,:,:,:) ,3,3,ngp,nel,nsteps);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hist = Compute_mat_nodal(hist, num, Ninv)
temp_strss = zeros(3,3,num.nen,num.el,num.steps);
mat_perm = permute(hist.strss.mat,[3,1,2,4,5]);

% Extrapolate from gp to nodes
for i= 1:num.nen
   temp_strss(:,:,i,:,:) = squeeze(sum(Ninv(i,:)'.*mat_perm, 1));
end

% Average values at the nodes
for i=1:num.np
   sm = 0;
   for j=1:hist.nconn(i,end)
      el = hist.nconn(i,j);
      node = hist.conn(el,:)==i;
      
      sm = sm + temp_strss(:,:,node,el,:);
   end
   hist.strss.mat_n(:,:,i,:) = sm./hist.nconn(i,end);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hist = Compute_principal(hist,num)
for k=1:num.steps
   for i=1:num.np
      hist.strss.P(:,i,k)= eig(hist.strss.mat_n(:,:,i,k));
   end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hist = Compute_von(hist)
P1=squeeze(hist.strss.P(1,:,:));
P2=squeeze(hist.strss.P(2,:,:));
P3=squeeze(hist.strss.P(3,:,:));

hist.strss.vm = ( 0.5*((P1-P2).^2 + (P2-P3).^2 + (P1-P3).^2) ).^(1/2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hist = Compute_avg(hist,num)

S = reshape(hist.strss.mat ,9,num.gp,num.el,num.steps);

hist.strss.avg = squeeze(mean(S([1,5,9,2,6,3],:,:,:),2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFvD(hist,last)
D = zeros(last,1);

sumForce = squeeze(sum(hist.force))';
[~,f_indx]  = max(sumForce(last,:));
U = squeeze(hist.disp(:,f_indx,:));

for i=1:last
   u = U(:,i);
   if ( max(u)== max(abs(u)) )
      u_face = u>=0.98.*max(u);
   elseif (abs(min(u)) == max(abs(u)))
      u_face = u<=0.98.*min(u);
   end
   D(i) = mean(u(u_face));
end

plot(D, sumForce(:,f_indx),'r','LineWidth',2)
title('Force (total surface not nodal) vs displacement')
xlabel('displacement')
ylabel('force')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotDeformation(hist,num, ax)

conn = hist.conn(1:num.els,1:num.nen);
title('Undeformed vs Deformed')

PlotMesh(hist.eltype, squeeze(hist.coor(:,:,  1)), conn, num.nen, num.els, 'g', ax)
PlotMesh(hist.eltype, squeeze(hist.coor(:,:,end)), conn, num.nen, num.els, 'r', ax)

ax.XLim = [min(min(hist.coor(:,1,:))) max(max(hist.coor(:,1,:)))];
ax.YLim = [min(min(hist.coor(:,2,:))) max(max(hist.coor(:,2,:)))];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotcont(hist, num, contH,subH)
title('Stress contour GUI');
set(subH,'Position',[0.6,0.02,0.35,0.75])
axis equal
if num.nen == 9
   conn = hist.conn(:,[1 5 2 6 3 7 4 8])';
elseif num.nen == 6
   conn = hist.conn(:,[1,4,2,5,3,6])';
else
   conn = hist.conn(1:num.els,1:num.nen)';
end
if num.ndm == 2
   X = reshape(hist.coor(conn,1,1),   num.nen,num.els);
   Y = reshape(hist.coor(conn,2,1),   num.nen,num.els);
   C = reshape(hist.strss.P(3,conn,1),num.nen,num.els);
   patchH = patch(X,Y,C, 'FaceAlpha',.7,'EdgeColor','r',...
      'Marker', 'o','MarkerFaceColor','r','MarkerSize',3);
   fac = conn;
elseif num.ndm == 3
   if strcmp(hist.eltype,'Q8')
      fac = zeros(4,num.els*6);
      for i = 1:num.els
         fac(:,(i-1)*6+(1:6)) = [...
            conn(1:4,i), conn(5:8,i),...
            conn([1,2,6,5],i), conn([3,4,8,7],i),...
            conn([4,1,5,8],i), conn([5,6,7,8],i)];
      end
   elseif strcmp(hist.eltype,'T4')
      fac = zeros(3,num.els*4);
      for i = 1:num.els
         fac(:,(i-1)*4+(1:4)) = [conn([1,2,3],i), conn([1,2,4],i),conn([1,3,4],i), conn([2,3,4],i)];
      end
   elseif strcmp(hist.eltype,'T10')
      fac = zeros(6,num.els*4);
      for i = 1:num.els
         fac(:,(i-1)*4+(1:4)) = [...
            conn([1,7,3,6,2,5],i),conn([1,7,3,10,4,8],i),conn([2,9,4,8,1,5],i),conn([2,6,3,10,4,9],i)];
      end
   end
   X = reshape(hist.coor(fac,1,1),   size(fac,1), size(fac,2));
   Y = reshape(hist.coor(fac,2,1),   size(fac,1), size(fac,2));
   Z = reshape(hist.coor(fac,3,1),   size(fac,1), size(fac,2));
   C = reshape(hist.strss.P(3,fac,1),size(fac,1), size(fac,2));
   patchH = patch(X,Y,Z,C, 'FaceAlpha',.7,'EdgeColor','r',...
      'Marker', 'o','MarkerFaceColor','r','MarkerSize',3);
   view (gca,[0.2 -0.4 0.1])
   set(gca,'DataAspectRatio',[1 1 1],'PlotBoxAspectRatio',[570.5 570.5 570.5])
end
cbH = colorbar;
cbH.Limits = [min(min(hist.strss.P(3, :,:))) max(max(hist.strss.P(3, :,:)))];

colormap jet
subH.XLim = [min(min(hist.coor(:,1,:))) max(max(hist.coor(:,1,:)))];
subH.YLim = [min(min(hist.coor(:,2,:))) max(max(hist.coor(:,2,:)))];
uifunctions(contH, patchH, cbH, hist, num, fac)
end

function uifunctions(hFig, patchH, colorbarH, hist, num, conn)

sliderH=uicontrol('Style','slider','Units','normalized',...
   'Position',[0.75,0.8,0.15,0.03],'Value',1,...
   'Max',num.steps,'Min',1,'SliderStep',[1/(num.steps-1) 2/(num.steps-1)],...
   'Visible','on');

slidertext = text(0.65,1.3, '$step~~0$','FontSize',16,'Units','normalized','Interpreter','latex' );

legH = text(1.04,1.05,'$\sigma_{1}$','FontSize',25,'Units','normalized','Interpreter','latex' );

names = {'Principle stress 1', 'Principle stress 2', 'Principle stress 3',...
   'stress 11 ','stress 22 ','stress 33 ',...
   'stress 12 ','stress 23 ','stress 13 ',...
   'Von mises stress '};
listH=uicontrol('Parent',hFig,'Style','listbox','Units','normalized',...
   'String',names,'Position', [0.6 0.77 0.1 0.19], 'Visible','on');

sliderH.Callback =@(sliderH ,event) ChangeStep(sliderH, slidertext, patchH, hist, num, conn);
listH.Callback =@(listH,event) ChangeStress(listH, legH, colorbarH, patchH, hist, slidertext, num, conn);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChangeStep(sliderH, slidertext, h, hist, num, conn)
i = sliderH.Value;
step = int8(i);
if (step < 1.0001*i && step > 0.9999*i)
   
   h.XData = reshape(hist.coor(conn,1,step),size(conn,1),size(conn,2));
   h.YData = reshape(hist.coor(conn,2,step),size(conn,1),size(conn,2));
   h.CData = reshape(hist.strss.P(3,conn,step),size(conn,1),size(conn,2));
   if num.ndm == 3
      h.ZData = reshape(hist.coor(conn,3,step),size(conn,1),size(conn,2));
   end
   
   slidertext.String = ['$step~~',num2str(step-1),'$'];
   axis equal
else
   error('check step size of the slider');
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ChangeStress(listH, legH, colorbarH, h, hist, slidertext, conn)
i = listH.Value;
indx1 = strfind(slidertext.String,'~');
indx2 = strfind(slidertext.String,'$');
step = str2double(slidertext.String(indx1(end)+1:indx2(end)-1))+1;
switch i
   case 1
      h.CData = reshape(hist.strss.P(3,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{1}$';
      L_limit = min(min(hist.strss.P(1, :,:)));
      U_limit = max(max(hist.strss.P(1, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 2
      h.CData = reshape(hist.strss.P(2,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{2}$';
      L_limit = min(min(hist.strss.P(1, :,:)));
      U_limit = max(max(hist.strss.P(1, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 3
      h.CData = reshape(hist.strss.P(1,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{3}$';
      L_limit = min(min(hist.strss.P(1, :,:)));
      U_limit = max(max(hist.strss.P(1, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 4
      h.CData = reshape(hist.strss.mat_n(1,1,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{11}$';
      L_limit = min(min(hist.strss.mat_n(1,1, :,:)));
      U_limit = max(max(hist.strss.mat_n(1,1, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 5
      h.CData = reshape(hist.strss.mat_n(2,2,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{22}$';
      L_limit = min(min(hist.strss.mat_n(2,2, :,:)));
      U_limit = max(max(hist.strss.mat_n(2,2, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 6
      h.CData = reshape(hist.strss.mat_n(3,3,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{33}$';
      L_limit = min(min(hist.strss.mat_n(3,3, :,:)));
      U_limit = max(max(hist.strss.mat_n(3,3, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 7
      h.CData = reshape(hist.strss.mat_n(1,2,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{12}$';
      L_limit = min(min(hist.strss.mat_n(1,2, :,:)));
      U_limit = max(max(hist.strss.mat_n(1,2, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 8
      h.CData = reshape(hist.strss.mat_n(2,3,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{23}$';
      L_limit = min(min(hist.strss.mat_n(2,3, :,:)));
      U_limit = max(max(hist.strss.mat_n(2,3, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 9
      h.CData = reshape(hist.strss.mat_n(1,3,conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{13}$';
      L_limit = min(min(hist.strss.mat_n(1,3, :,:)));
      U_limit = max(max(hist.strss.mat_n(1,3, :,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit==L_limit)*0.001];
   case 10
      h.CData = reshape(hist.strss.vm(conn,step),size(conn,1),size(conn,2));
      legH.String = '$\sigma_{vm}$';
      L_limit = min(min(hist.strss.vm(:,:)));
      U_limit = max(max(hist.strss.vm(:,:)));
      colorbarH.Limits = [L_limit U_limit+(U_limit>L_limit)*0.001];
end
end