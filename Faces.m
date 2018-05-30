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
               obj.faceList = zeros(obj.numFaces, 2);
               for i = 1:numel
                  obj.faceList(1+(4*(i-1): 4*i-1),:) =[...
                     el(i,[1,2]); el(i,[2,3]); el(i,[4,3]); el(i,[1,4])];
               end
               num = struct('el',obj.numFaces, 'nen', 2, 'gp', 2, 'ndm', 1);
               obj.gp = L2(num, 0);
            case 'Q9'
               obj.numFaces = 4*numel;
               obj.numGP    = 3;
               obj.faceList = zeros(obj.numFaces, 3);
               for i = 1:numel
                  obj.faceList(1+(4*(i-1): 4*i-1),:) = [...
                     el(i,[1,5,2]); el(i,[2,6,3]); el(i,[3,7,4]); el(i,[4,8,1])];
               end
               num = struct('el',obj.numFaces, 'nen', 3, 'gp', 3, 'ndm', 1);
               obj.gp = L3(num, 0);
            case 'T3'
               obj.numFaces = 3*numel;
               obj.numGP    = 2;
               obj.faceList = zeros(obj.numFaces, 2);
               for i = 1:numel
                  obj.faceList(1+(3*(i-1): 3*i-1),:) = [...
                     el(i,[1,2]); el(i,[2,3]); el(i,[3,1])];
               end
               num = struct('el',obj.numFaces, 'nen', 2, 'gp', 2, 'ndm', 1);
               obj.gp = L2(num, 0);
            case 'T6'
               obj.numFaces = 3*numel;
               obj.numGP    = 3;
               obj.faceList = zeros(obj.numFaces, 3);
               for i = 1:numel
                  obj.faceList(1+(3*(i-1): 3*i-1),:) = [...
                     el(i, [1,4,2]); el(i, [2,5,3]); el(i, [3,6,1])];
               end
               num = struct('el',obj.numFaces, 'nen', 3, 'gp', 3, 'ndm', 1);
               obj.gp = L3(num, 0);
            case 'Q8'
               obj.numFaces = 6*numel;
               obj.numGP    = 4;
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
               num = struct('el',obj.numFaces, 'nen', 4, 'gp', 4, 'ndm', 2); 
               obj.gp  = Q4(num, 0);
         end
      end
      %% get functions
      function value = get.faces(obj)
         comb  = [];
         value = [];
         c = nchoosek(obj.nodes, obj.numGP);
         for i = 1:size(c, 1)
           cc   = perms(c(i,:)); 
           comb = [comb; cc];
         end
         for i = 1:size(comb,1)
            t = obj.faceList( sum( obj.faceList==comb(i,:), 2)==obj.numGP, :);
            value = [value; t];
         end 
      end
      
      function value = get.indices(obj)
         value = zeros(obj.numGP,1);
         face = obj.faces(obj.ifc,:)';
         for i = 1:obj.numGP
           value(i)= find(obj.nodes == face(i));
         end
      end
   end
end

