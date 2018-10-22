%% integrate unmatching interfaces
function [xlintL, xlintR, drdrL, drdrR, xiL, xiR] = intBounds2(coorL, coorR, xiL, xiR, ndm)
if ndm == 2
   xlintL = zeros(size(coorL));
   xlintR = zeros(size(coorR));
   
   lengthL = norm(diff(coorL(1:2,:),1));
   lengthR = norm(diff(coorR(1:2,:),1));
   % Find the larger element to be shortend:
   if	lengthL>lengthR
      coorB = coorL;	xlintB = xlintL;
      coorS = coorR; xlintS = coorR;
   elseif lengthR>lengthL
      coorB = coorR; xlintB = xlintR;
      coorS = coorL; xlintS = coorL;
   end
   nelB = size(coorB,2);
   nelS = size(coorS,2);
   eS1 = -1;
   eS2 =  1;
   % There are two points in the interface, check if either or both of the larger don't match the shorter element, and force match them
   % Point 1:
   if norm(coorB(:,1) - coorS(:,2))<1e-8
      xlintB(:,[1,4]) = coorB(:,[1,4]);
      eB1 = -1;
   else
      POUxi = POU_Coord(coorS(1,2),coorS(2,2),coorB,1,nelB);
      eB1 = POUxi(1);
      xlintB(:,1) = coorS(:,2);
      N = shapeF([eB1,1], nelB);
      xlintB(:,4) = coorB*N;
   end
   % Point 2:
   if norm(coorB(:,2) - coorS(:,1))<1e-8
      xlintB(:,[2,3]) = coorB(:,[2,3]);
      eB2 = 1;
   else
      POUxi = POU_Coord(coorS(1,1),coorS(2,1),coorB,1,nelB);
      eB2 = POUxi(1);
      xlintB(:,2) = coorS(:,1);
      N = shapeF([eB2,1], nelB);
      xlintB(:,3) = coorB*N;
   end
   
   % shorten support for bubble, Big
   lside23 = norm(xlintB(:,2)-xlintB(:,3));
   lside14 = norm(xlintB(:,1)-xlintB(:,4));
   lside12 = norm(xlintB(:,2)-xlintB(:,1));
   shortyn = (lside12 < lside23) && (lside12 < lside14);
   if shortyn && lside14 < lside23
      lintratio = lside12/lside14;
      lintxy1 = (xlintB(:,4)-xlintB(:,1))*lintratio + xlintB(:,1);
      POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintB,0,nelB);
      
      xlintB(:,4) = lintxy1;
      shl = shapeF([+1,POUxi(2)], nelB);
      lintxy2 = xlintB*shl;
      xlintB(:,3) = lintxy2;
   elseif shortyn && lside14 >= lside23
      lintratio = lside12/lside23;
      lintxy1 = (xlintB(:,3)-xlintB(:,2))*lintratio + xlintB(:,2);
      POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintB,0,nelB);
      xlintB(:,3) = lintxy1;
      shl = shapeF([-1,POUxi(2)], nelB);
      lintxy2 = xlintB*shl;
      xlintB(:,4) = lintxy2;
   end
   
   % shorten support for bubble, Small
   lside23 = norm(xlintS(:,2)-xlintS(:,3),2);
   lside14 = norm(xlintS(:,1)-xlintS(:,4),2);
   lside12 = norm(xlintS(:,2)-xlintS(:,1),2);
   shortyn = (lside12 < lside23) && (lside12 < lside14);
   if shortyn && lside14 < lside23
      lintratio = lside12/lside14;
      lintxy1 = (xlintS(:,4)-xlintS(:,1))*lintratio + xlintS(:,1);
      POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintS,0,nelS);
      xlintS(:,4) = lintxy1;
      shl = shapeF([1,POUxi(2)],nelS);
      lintxy2 = xlintS*shl;
      xlintS(:,3) = lintxy2;
   elseif shortyn && lside14 >= lside23
      lintratio = lside12/lside23;
      lintxy1 = (xlintS(:,3)-xlintS(:,2))*lintratio + xlintS(:,2);
      POUxi = POU_Coord(lintxy1(1),lintxy1(2),xlintS,0,nelS);
      xlintS(:,3) = lintxy1;
      shl = shapeF([-1,POUxi(2)],nelS);
      lintxy2 = xlintS*shl;
      xlintS(:,4) = lintxy2;
   end
   
   if	lengthL>lengthR
      xlintL =	xlintB;
      xlintR = xlintS;
      eL = [eB1 eB2];
      eR = [eS1 eS2];
   elseif lengthR>lengthL
      xlintR = xlintB;
      xlintL = xlintS;
      eR = [eB1 eB2];
      eL = [eS1 eS2];
   end
   if nelB == 4
      drL = 2;
      drR = 2;
      ro = -1;
   elseif nelB == 3
      drL = 1;
      drR = 1;
      ro = 0;
   end
   drdrL = (eL(2) - eL(1))/drL;
   drdrR = (eR(2) - eR(1))/drR;
   m  = (eR(2) - eR(1))/(eL(1)-eL(2));
   rL = drdrL*(xiL(:,1)-ro) + eL(1);
   rR = m*(rL-eL(2)) + eR(1);
   
   xiL = [rL xiL(:,2)];
   xiR = [rR xiR(:,2)];
elseif ndm == 3
   xlintL = coorL;
   xlintR = coorR;
   drdrL = 1;
   drdrR = 1;
end
end