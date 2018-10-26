function [gRot, schmidt] = gRotSchmidt(angles, bi, ni, nslip)
cr = cosd(angles(1));	cs = cosd(angles(2));   ct = cosd(angles(3));
sr = sind(angles(1));   ss = sind(angles(2));   st = sind(angles(3));

gRot = [...
   -sr*st-cr*ct*cs,	+cr*st-sr*ct*cs, ct*ss
   +sr*ct-cr*st*cs,	-cr*ct-sr*st*cs, st*ss
   +cr*ss         ,	+sr*ss         , cs   ];

bs = bi * gRot;
ns = ni * gRot;
schmidt = zeros(3,3,nslip);
for s = 1:nslip
   schmidt(:,:,s) = bs(s,:)'*ns(s,:);
end
end

