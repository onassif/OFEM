function ulres = getInitU(mult, inpt, el)
BC    = inpt.BC;
ulres = zeros(size(el.conn,1),size(el.conn,2)*el.ndm);
for i = 1:size(el.conn,1)
   for j = 1:size(el.conn,2)
      indc = find(el.conn(i,j)== BC(:,1));
      for k=1:size(indc)
         ulres(i, el.ndm*(j-1)+BC(indc(k),2)) = mult*BC(indc(k),3);
      end
   end
end
% ulres = ulres';
end

