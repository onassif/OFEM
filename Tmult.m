function [C] = Tmult(A,B)
%function to multiply tensors
dimA = ndims(A);
dimB = ndims(B);

if (dimA==4 && dimB==2 && size(B,2)~=1)
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm3A = size(A,3);
   ndm4A = size(A,4);
   ndm1B = size(B,1);
   ndm2B = size(B,2);
   C = zeros(ndm1A,ndm2A);
   if ~(ndm3A==ndm1B && ndm4A==ndm2B)
      error('wrong dimensions');
   end
   
   for i=1:ndm1A
      for j=1:ndm2A
         for k=1:ndm3A
            for l=1:ndm4A
               C(i,j) = C(i,j) + A(i,j,k,l)*B(k,l);
            end
         end
      end
   end
   
elseif (dimA==4 && dimB==3)
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm3A = size(A,3);
   ndm4A = size(A,4);
   ndm1B = size(B,1);
   ndm2B = size(B,2);
   ndm3B = size(B,3);
   C = zeros(ndm1A,ndm2A,ndm3B);
   if ~(ndm3A==ndm1B && ndm4A==ndm2B)
      error('wrong dimensions');
   end
   
   for i=1:ndm1A
      for j=1:ndm2A
         for k=1:ndm3B
            for l=1:ndm3A
               for m=1:ndm4A
                  C(i,j,k) = C(i,j,k) + A(i,j,l,m)*B(l,m,k);
               end
            end
         end
      end
   end
elseif(dimA==3 && dimB==2 && size(B,2)~=1)
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm3A = size(A,3);
   ndm1B = size(B,1);
   ndm2B = size(B,2);
   C = zeros(ndm1A,ndm2B,ndm3A);
   if ~(ndm2A==ndm1B)
      error('wrong dimensions');
   end
   
   for i =1:ndm1A
      for j=1:ndm2B
         for k=1:ndm3A
            for l=1:ndm2A
               C(i,j,k) = C(i,j,k) + A(i,l,k)*B(l,j);
            end
         end
      end
   end
   
elseif(dimA==2 && size(A,2)~=1 && dimB==3)
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm1B = size(B,1);
   ndm2B = size(B,2);
   ndm3B = size(B,3);
   C = zeros(ndm1A,ndm2B,ndm3B);
   if ~(ndm2A==ndm1B)
      error('wrong dimensions');
   end
   
   for i =1:ndm1A
      for j=1:ndm2B
         for k=1:ndm3B
            for l=1:ndm2A
               C(i,j,k) = C(i,j,k) + A(i,l)*B(l,j,k);
            end
         end
      end
   end
   
elseif(dimA==3 && dimB==2 && size(B,2) == 1 && size(A,3)==size(B,1))
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm3A = size(A,3);
   ndm1B = size(B,1);
   C = zeros(ndm1A,ndm2A,ndm3A);
   
   for s =1:ndm1B
      C(:,:,s) = C(:,:,s) + A(:,:,s)*B(s);
   end
   
 elseif(dimA==2 && dimB==2 && size(B,2) == 1 && size(A,2) ~= 1)
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm1B = size(B,1);
   C = zeros(ndm1A,ndm2A,ndm1B);
   
   for s =1:ndm1B
      C(:,:,s) = C(:,:,s) + A(:,:)*B(s);
   end  
   
elseif(dimA==3 && dimB==3 && size(A,3)==size(B,3))
   %%
   ndm1A = size(A,1);
   ndm2A = size(A,2);
   ndm3A = size(A,3);
   ndm1B = size(B,1);
   ndm2B = size(B,2);
   C = zeros(ndm1A,ndm2A,ndm1B,ndm2B);
   
   for i=1:ndm1A
      for j=1:ndm2A
         for k=1:ndm1B
            for l=1:ndm2B
               for s =1:ndm3A
                  C(i,j,k,l) = C(i,j,k,l) + A(i,j,s)*B(k,l,s);
               end
            end
         end
      end
   end
end
end