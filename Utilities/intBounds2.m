%% integrate unmatching interfaces
function [xlintL,xlintR, drdrL,drdrR, eGPL,eGPR] = intBounds2(coorL,coorR, sGPL,sGPR, eGPL,eGPR)
xiL = eGPL.xi;
xiR = eGPR.xi;
ndm = size(coorL,1);
nen = size(coorL,2);
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
   
   [sGPL.det_dXdxi_list, sGPL.dNdX_list, sGPL.dXdxi_list] = shapeRef([eL(1);(eL(1)+eL(2))/2;eL(2)], 1:3, sGPL.dNdxi_list);
   [sGPR.det_dXdxi_list, sGPR.dNdX_list, sGPR.dXdxi_list] = shapeRef([eR(1);(eR(1)+eR(2))/2;eR(2)], 1:3, sGPR.dNdxi_list);
   
   drdrL = abs(sGPL.det_dXdxi_list(eGPL.i));
   drdrR = abs(sGPR.det_dXdxi_list(eGPR.i));
   %m  = (eR(2) - eR(1))/(eL(1)-eL(2));
   rL = drdrL*(xiL(1,:)-ro) + eL(1);
   rR = drdrR*(eL(2)   -rL) + eR(1);
   
   eGPL.xi = [rL; xiL(2,:)];
   eGPR.xi = [rR; xiR(2,:)];
elseif ndm == 3
   if nen == 8
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
      
      [sGPL.det_dXdxi_list, sGPL.dNdX_list, sGPL.dXdxi_list] = shapeRef(eL', 1:4, sGPL.dNdxi_list);
      [sGPR.det_dXdxi_list, sGPR.dNdX_list, sGPR.dXdxi_list] = shapeRef(eR', 1:4, sGPR.dNdxi_list);
      
%       areaL = polyarea(eL(1,:)',eL(2,:)');
%       areaR = polyarea(eR(1,:)',eR(2,:)');
%       drdrL = areaL/drL;
%       drdrR = areaR/drR;
      drdrL = abs(sGPL.det_dXdxi_list(eGPL.i));
      drdrR = abs(sGPR.det_dXdxi_list(eGPR.i));
      
      rL = drdrL*(xiL(1:2,:)-ro) + eL;
      
      gpR = Q4(0);
      gpR = gpR.shapeIso();
      [gpR.det_dXdxi_list, gpR.dNdX_list, gpR.dXdxi_list] = shapeRef(...
         [-1,-1; +1,-1; +1,+1; -1,+1], 1:4, gpR.dNdxi_list);

      rR = eR*gpR.Nmat;
      rR = rR(:,[1,4,2,3]);

%       mx  = abs(eR(:,1)-eR(:,2))/2+ eR(:,1);
%       mxl = norm(mx - eR(:,1));
%       my  = abs(eR(:,1)-eR(:,4))/2+ eR(:,1);
%       myl = norm(my - eR(:,1));
%       sq3 = 1/sqrt(3);
%       
%       rR = [...
%          mx(1)-sq3*mxl, mx(1)-sq3*mxl, mx(1)+sq3*mxl, mx(1)+sq3*mxl
%          my(2)-sq3*myl, my(2)+sq3*myl, my(2)-sq3*myl, my(2)+sq3*myl];
%       
% %       polyin = polyshape(eR(1,:),eR(2,:));
% %       [mx,my] = centroid(polyin);
% 
%       m12 = (eR(:,1)+eR(:,2))/2;
%       m23 = (eR(:,2)+eR(:,3))/2;
%       m34 = (eR(:,3)+eR(:,4))/2;
%       m41 = (eR(:,4)+eR(:,1))/2;
%       
%       pD1 = m12 + sq3*(eR(:,1)-m12);
%       pD2 = m12 + sq3*(eR(:,2)-m12);
%       pR2 = m23 + sq3*(eR(:,2)-m23);
%       pR3 = m23 + sq3*(eR(:,3)-m23);
%       pU3 = m34 + sq3*(eR(:,3)-m34);
%       pU4 = m34 + sq3*(eR(:,4)-m34);
%       pL4 = m41 + sq3*(eR(:,4)-m41);
%       pL1 = m41 + sq3*(eR(:,1)-m41);
%       
%       x1 = pD1(1); y1 = pD1(2);
%       x2 = pU4(1); y2 = pU4(2);
%       x3 = pR2(1); y3 = pR2(2);
%       x4 = pL1(1); y4 = pL1(2);
%       g1 =  [y2-y1, x1-x2; y4-y3, x3-x4] \ [x1*y2-x2*y1; x3*y4-x4*y3];
%       
%       x1 = pD2(1); y1 = pD2(2);
%       x2 = pU3(1); y2 = pU3(2);
%       x3 = pR2(1); y3 = pR2(2);
%       x4 = pL1(1); y4 = pL1(2);
%       g2 =  [y2-y1, x1-x2; y4-y3, x3-x4] \ [x1*y2-x2*y1; x3*y4-x4*y3];
%       
%       x1 = pD2(1); y1 = pD2(2);
%       x2 = pU3(1); y2 = pU3(2);
%       x3 = pR3(1); y3 = pR3(2);
%       x4 = pL4(1); y4 = pL4(2);
%       g3 =  [y2-y1, x1-x2; y4-y3, x3-x4] \ [x1*y2-x2*y1; x3*y4-x4*y3];
%       
%       x1 = pD1(1); y1 = pD1(2);
%       x2 = pU4(1); y2 = pU4(2);
%       x3 = pR3(1); y3 = pR3(2);
%       x4 = pL4(1); y4 = pL4(2);
%       g4 =  [y2-y1, x1-x2; y4-y3, x3-x4] \ [x1*y2-x2*y1; x3*y4-x4*y3];
%       
%       m = sum(eR,2)/4;
% %       m1_m = norm(m-eR(:,1));    m1_v = (m-eR(:,1))/m1_m;
% %       m2_m = norm(m-eR(:,2));    m2_v = (m-eR(:,2))/m2_m;
% %       m3_m = norm(m-eR(:,3));    m3_v = (m-eR(:,3))/m3_m;
% %       m4_m = norm(m-eR(:,4));    m4_v = (m-eR(:,4))/m4_m;
% %       
%       m1= eR(:,1) - m;
%       m3= eR(:,2) - m;
%       m4= eR(:,3) - m;
%       m2= eR(:,4) - m;
%       
%       rR = m + sq3*[m1, m2, m3, m4];
%       rR= [g1 g4 g2 g3];
      
      
      eGPL.xi = [rL; xiL(3,:)];
      eGPR.xi = [rR; xiR(3,:)];
   elseif nen == 4 || nen == 10
      xlintL = coorL;
      xlintR = coorR;
      drdrL = 1;
      drdrR = 1;
   end
end

eGPL = eGPL.shapeIso;
eGPR = eGPR.shapeIso;
end