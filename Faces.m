classdef Faces
   properties (SetAccess = public, GetAccess = public)
      nodes
      coor
      gp
      ifc
   end
   properties (SetAccess = private, GetAccess = public)
      numFaces
      numGP
      numEN
      faces
      indices
      faceList
   end
   properties (SetAccess = private, Hidden)
      elements
   end
   methods
      %% Construct
      function obj = Faces(elmtype, el, numel)
         switch elmtype
            case 'Q4'
               obj.numFaces = 4*numel;
               obj.numGP    = 2;
               obj.numEN    = 2;
               obj.faceList = zeros(obj.numFaces, 2);
               for i = 1:numel
                  obj.faceList(1+(4*(i-1): 4*i-1),:) =[...
                     el(i,[1,2]); el(i,[2,3]); el(i,[4,3]); el(i,[1,4])];
               end
               obj.gp = L2(0);
            case 'Q9'
               obj.numFaces = 4*numel;
               obj.numGP    = 3;
               obj.numEN    = 3;
               obj.faceList = zeros(obj.numFaces, 3);
               for i = 1:numel
                  obj.faceList(1+(4*(i-1): 4*i-1),:) = [...
                     el(i,[1,5,2]); el(i,[2,6,3]); el(i,[3,7,4]); el(i,[4,8,1])];
               end
               obj.gp = L3(0);
            case 'T3'
               obj.numFaces = 3*numel;
               obj.numGP    = 2;
               obj.numEN    = 2;
               obj.faceList = zeros(obj.numFaces, 2);
               for i = 1:numel
                  obj.faceList(1+(3*(i-1): 3*i-1),:) = [...
                     el(i,[1,2]); el(i,[2,3]); el(i,[3,1])];
               end
               obj.gp = L2(0);
            case 'T6'
               obj.numFaces = 3*numel;
               obj.numGP    = 3;
               obj.numEN    = 3;
               obj.faceList = zeros(obj.numFaces, 3);
               for i = 1:numel
                  obj.faceList(1+(3*(i-1): 3*i-1),:) = [...
                     el(i, [1,4,2]); el(i, [2,5,3]); el(i, [3,6,1])];
               end
               obj.gp = L3(0);
            case 'T4'
               obj.numFaces = 4*numel;
               obj.numGP    = 1;
               obj.numEN    = 3;
               obj.faceList = zeros(obj.numFaces, 3);
               for i = 1:numel
                  obj.faceList((i-1)*4+(1:4),:) = [...
                     el(i, [1,2,3]) 
                     el(i, [1,2,4]) 
                     el(i, [1,3,4])
                     el(i, [2,3,4])];
               end
               obj.gp = T3(0);
            case {'Q8'}
               obj.numFaces = 6*numel;
               obj.numGP    = 4;
               obj.numEN    = 4;
               obj.faceList = zeros(obj.numFaces, 4);
               for i = 1:numel
                  obj.faceList(1+(6*(i-1): 6*i-1),:) =[...
                     el(i,[1,2,3,4])
                     el(i,[2,3,7,6])
                     el(i,[5,6,7,8])
                     el(i,[1,4,8,5])
                     el(i,[2,1,5,6])
                     el(i,[3,4,8,7])];
               end
               obj.gp  = Q4(0);
         end
         obj.gp = obj.gp.shapeIso();
      end
      %% get functions
      function value = get.faces(obj)
         indc = false(size(obj.faceList,1),1);
         for i = 1:size(obj.faceList,1)
            if sum(sum(obj.faceList(i,:) == obj.nodes)) == obj.numEN
               indc(i) = true;
            end
         end
         value = obj.faceList(indc,:);
      end
      
      function value = get.indices(obj)
         value = zeros(obj.numEN,1);
         face = obj.faces(obj.ifc,:)';
         for i = 1:obj.numEN
            value(i)= find(obj.nodes == face(i));
         end
      end
   end
end

