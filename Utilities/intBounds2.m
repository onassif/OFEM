%% integrate unmatching interfaces
function [xlintL,xlintR, drdrL,drdrR, eGPL,eGPR] = intBounds2(coorL,coorR, eGPL,eGPR)
xiL = eGPL.xi;
xiR = eGPR.xi;
ndm = size(coorL,1);
if ndm == 2
   xlintL = zeros(size(coorL));
   xlintR = zeros(size(coorR));
   
   lengthL = norm(diff(coorL(:,1:2),2));
   lengthR = norm(diff(coorR(:,1:2),2));
   % Find the larger element to be shortend:
   if	lengthL>lengthR
      coorB = coorL;	xlintB = xlintL;
      coorS = coorR; xlintS = coorR;
   elseif lengthR>=lengthL
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
   elseif lengthR>=lengthL
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
   rL = drdrL*(xiL(1,:)-ro) + eL(1);
   rR = m*(rL-eL(2)) + eR(1);
   
   eGPL.xi = [rL; xiL(2,:)];
   eGPR.xi = [rR; xiR(2,:)];
elseif ndm == 3
   xlintL = zeros(size(coorL));
   xlintR = zeros(size(coorR));
   
   surfL = coorL(2:3,1:4)';
   surfR = coorR(2:3,1:4)';
   areaL = polyarea(surfL(:,1), surfL(:,2));
   areaR = polyarea(surfR(:,1), surfR(:,2));
   % Find the larger element to be shortend:
   if	areaL>areaR
      coorB = coorL;	xlintB = xlintL;
      coorS = coorR; xlintS = coorR;
   elseif areaR>=areaL
      coorB = coorR; xlintB = xlintR;
      coorS = coorL; xlintS = coorL;
   end
   nelB = size(coorB,2);

   eS1 = [-1;-1];
   eS2 = [-1; 1];
   eS3 = [ 1; 1];
   eS4 = [ 1;-1];
   % There are two points in the interface, check if either or both of the larger don't match the shorter element, and force match them
   % Point 1:
   if norm(coorB(:,1) - coorS(:,1))<1e-8
      xlintB(:,[1,5]) = coorB(:,[1,5]);
      eB1 = [-1;-1];
   else
      POUxi = POU_Coord3(coorS(1,1),coorS(2,1),coorS(3,1),coorB,1,nelB);
      eB1 = POUxi(1:2);
      xlintB(:,1) = coorS(:,1);
      N = shapeF([eB1;1], nelB);
      xlintB(:,5) = coorB*N;
   end
   % Point 2:
   if norm(coorB(:,2) - coorS(:,4))<1e-8 
      xlintB(:,[2,6]) = coorB(:,[2,6]);
      eB2 = [1;-1];
   else
      POUxi = POU_Coord3(coorS(1,4),coorS(2,4),coorS(3,4),coorB,1,nelB);
      eB2 = POUxi(1:2);
      xlintB(:,2) = coorS(:,4);
      N = shapeF([eB2;1], nelB);
      xlintB(:,6) = coorB*N;
   end
   % Point 3:
   if norm(coorB(:,3) - coorS(:,3))<1e-8
      xlintB(:,[3,7]) = coorB(:,[3,7]);
      eB3 = [1;1];
   else
      POUxi = POU_Coord3(coorS(1,3),coorS(2,3),coorS(3,3),coorB,1,nelB);
      eB3 = POUxi(1:2);
      xlintB(:,3) = coorS(:,3);
      N = shapeF([eB3;1], nelB);
      xlintB(:,7) = coorB*N;
   end
   % Point 4:
   if norm(coorB(:,4) - coorS(:,2))<1e-8
      xlintB(:,[4,8]) = coorB(:,[4,8]);
      eB4 = [-1;1];
   else
      POUxi = POU_Coord3(coorS(1,2),coorS(2,2),coorS(3,2),coorB,1,nelB);
      eB4 = POUxi(1:2);
      xlintB(:,4) = coorS(:,2);
      N = shapeF([eB4;1], nelB);
      xlintB(:,8) = coorB*N;
   end   

   
   %step 1 - see CMAME Truster and Masud page 18
   lside12 = norm(xlintB(:,1)-xlintB(:,2),2);
   lside23 = norm(xlintB(:,2)-xlintB(:,3),2);
   lside34 = norm(xlintB(:,3)-xlintB(:,4),2);
   lside51 = norm(xlintB(:,5)-xlintB(:,1),2);
   lside62 = norm(xlintB(:,6)-xlintB(:,2),2);
   lside73 = norm(xlintB(:,7)-xlintB(:,3),2);
   lside84 = norm(xlintB(:,8)-xlintB(:,4),2);
   %step 2
   lsideb = max([lside12,lside23,lside34]);
   lsidedc = min([lside51,lside62,lside73]);
   shortyn = lsideb < lsidedc;
   
   % shorten support for bubble
   if shortyn
      if abs(lside51 - lsidedc)<1e-14
         xlintB(:,5) = xlintB(:,1) + [lsideb;0;0];
      end
      if abs(lside62 - lsidedc)<1e-14
         xlintB(:,6) = xlintB(:,2) + [lsideb;0;0];
      end
      if abs(lside73 - lsidedc)<1e-14
         xlintB(:,7) = xlintB(:,3) + [lsideb;0;0];
      end
      if abs(lside84 - lsidedc)<1e-14
         xlintB(:,8) = xlintB(:,4) + [lsideb;0;0];
      end
   end
   
   
   %step 1 - see CMAME Truster and Masud page 18
   lside12 = norm(xlintS(:,1)-xlintS(:,2),2);
   lside23 = norm(xlintS(:,2)-xlintS(:,3),2);
   lside34 = norm(xlintS(:,3)-xlintS(:,4),2);
   lside51 = norm(xlintS(:,5)-xlintS(:,1),2);
   lside62 = norm(xlintS(:,6)-xlintS(:,2),2);
   lside73 = norm(xlintS(:,7)-xlintS(:,3),2);
   lside84 = norm(xlintS(:,8)-xlintS(:,4),2);
   %step 2
   lsideb = max([lside12,lside23,lside34]);
   lsidedc = min([lside51,lside62,lside73]);
   shortyn = lsideb < lsidedc;
   
   % shorten support for bubble
   if shortyn
      if abs(lside51 - lsidedc)<1e-14
         xlintS(:,5) = xlintS(:,1) - [lsideb;0;0];
      end
      if abs(lside62 - lsidedc)<1e-14
         xlintS(:,6) = xlintS(:,2) - [lsideb;0;0];
      end
      if abs(lside73 - lsidedc)<1e-14
         xlintS(:,7) = xlintS(:,3) - [lsideb;0;0];
      end
      if abs(lside84 - lsidedc)<1e-14
         xlintS(:,8) = xlintS(:,4) - [lsideb;0;0];
      end
   end
   
   if	areaL>areaR
      xlintL =	xlintB;
      xlintR = xlintS;
      eL = [eB1 eB2 eB3 eB4];
      eR = [eS1 eS2 eS3 eS4];
   elseif areaR>=areaL
      xlintR = xlintB;
      xlintL = xlintS;
      eR = [eB1 eB2 eB3 eB4];
      eL = [eS1 eS2 eS3 eS4];
   end
   if nelB == 8
      drL = 4;
      drR = 4;
      ro = [...
         -1,-1, 1, 1
         -1, 1, 1,-1];
   else
      error('Shape not implemented');
   end
   areaL = polyarea(eL(1,:)',eL(2,:)');
   areaR = polyarea(eR(1,:)',eR(2,:)');
   drdrL = areaL/drL;
   drdrR = areaR/drR;
   m  = areaR/-areaL;
   rL = drdrL*(xiL(1:2,:)-ro) + eL;
   
   mx  = abs(eR(:,1)-eR(:,2))/2+ eR(:,1);
   mxl = norm(mx - eR(:,1));
   my  = abs(eR(:,1)-eR(:,4))/2+ eR(:,1);
   myl = norm(my - eR(:,1));
   sq3 = 1/sqrt(3);
   
   rR = [...
   mx(1)-sq3*mxl, mx(1)-sq3*mxl, mx(1)+sq3*mxl, mx(1)+sq3*mxl
   my(2)-sq3*myl, my(2)+sq3*myl, my(2)-sq3*myl, my(2)+sq3*myl];
   
%    rR = -m*(rL-eL) + eR;
   
   eGPL.xi = [rL; xiL(3,:)];
   eGPR.xi = [rR; xiR(3,:)];
   
end

eGPL = eGPL.shapeIso;
eGPR = eGPR.shapeIso;
end