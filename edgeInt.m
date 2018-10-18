%% Integrate normal vectors
function [intedge, c1, nvect] = edgeInt(sGP, surfTan)
sGP.iel = 1;
ndm = size(surfTan,2)+1;
if ndm == 2
   J = norm(surfTan);
   n = cross([surfTan;0],[0;0;1] )/J;
   nvect = [...
      n(1) 0    n(2)
      0    n(2) n(1)];
elseif ndm == 3
   J = norm(cross(surfTan(:,2),surfTan(:,1)));
   n = cross(surfTan(:,2),surfTan(:,1))/J;
   nvect = [...
      n(1) 0    0    n(2) 0    n(3)
      0    n(2) 0    n(1) n(3) 0
      0    0    n(3) 0    n(2) n(1)];
end

c1 = J * sGP.weights;
intedge = sum(J * sGP.weights); % length of the edge
end