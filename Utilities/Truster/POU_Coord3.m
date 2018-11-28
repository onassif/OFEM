function POUxi = POU_Coord3(xc,yc,zc,xl,~,nel)
% Function to evaluate POU coordinates for 3D brick element in
% terms of the physical coordinates (xc,yc,zc). If the element is a
% parallelogram, the coordinates (xi,eta) are determined explicitly;
% otherwise, the relationship between (xc,yc) and (xi,eta) is nonlinear.

%solve nonlinear system
        
xp = [xc; yc; zc];
xi = zeros(3,1);
err = 1;
i = 0;

while abs(err) > 10^-11 && i < 20

    i = i + 1;
    
    if nel == 4 || nel == 10
      [shl,shld,~,~] = shltt(xi,nel,nel,0,0);
    elseif nel == 6
      [shl,shld,~,~] = shlw(xi,nel,nel,0,0);
    else
      [shl,shld,~,~] = shlb(xi,nel,nel,0,0);
    end
    x = xl(:,1:nel)*shl;
    xs = xl(:,1:nel)*shld;%zeros(2,2);

    del_xi = xs\(xp - x);
    err = norm(del_xi,inf);
    xi = xi + del_xi;
end

POUxi = xi;  
end