ident.twoD.I = eye(2);
ident.twoD.second.I4_sym = [...
    1   0   0   0
    0   1   0   0
    0   0 1/2   0
    0   0   0 1/2];
ident.twoD.second.I4_bulk= [...
    1   1   0   0
    1   1   0   0
    0   0   0   0
    0   0   0   0];
ident.twoD.second.I4_dev = [...
    2./3 -1/3   0   0
    -1/3  2/3   0   0
    0.00    0 1/2   0
    0.00    0   0 1/2];
ident.twoD.fourth.I4_sym  = zeros(2,2,2,2);
ident.twoD.fourth.I4_sym([1 16]) = 1;
ident.twoD.fourth.I4_sym([6 7 10 11]) = 1/2;

ident.twoD.fourth.I4_bulk = zeros(2,2,2,2);
ident.twoD.fourth.I4_bulk([1 4 13 16]) = 1;

ident.twoD.fourth.I4_dev  = zeros(2,2,2,2);
ident.twoD.fourth.I4_dev([1 16]) = 2/3;
ident.twoD.fourth.I4_dev([4 13]) =-1/3;
ident.twoD.fourth.I4_dev([6 7 10 11]) = 1/2;

ident.threeD.I = eye(3);
ident.threeD.second.I4_sym = [...
    1   0   0   0   0   0
    0   1   0   0   0   0
    0   0   1   0   0   0
    0   0   0 1/2   0   0
    0   0   0   0 1/2   0
    0   0   0   0   0 1/2];
ident.threeD.second.I4_bulk = [...
    1   1   1   0   0   0
    1   1   1   0   0   0
    1   1   1   0   0   0
    0   0   0   0   0   0
    0   0   0   0   0   0
    0   0   0   0   0   0];
ident.threeD.second.I4_dev = [...
    2./3 -1/3 -1/3   0   0   0
    -1/3  2/3 -1/3   0   0   0
    -1/3 -1/3  2/3   0   0   0
    0.00    0    0 1/2   0   0
    0.00    0    0   0 1/2   0
    0.00    0    0   0   0 1/2];

ident.threeD.fourth.I4_sym  = zeros(3,3,3,3);
ident.threeD.fourth.I4_sym([1 41 81])  = 1;
ident.threeD.fourth.I4_sym([11 13 21 25 29 31 51 53 57 61 69 71])  = 1/2;

ident.threeD.fourth.I4_bulk = zeros(3,3,3,3);
ident.threeD.fourth.I4_bulk([1 5 9 37 41 45 73 77 81]) = 1;

ident.threeD.fourth.I4_dev  = zeros(3,3,3,3);
ident.threeD.fourth.I4_dev([1 41 81])  =  2/3;
ident.threeD.fourth.I4_dev([5 9 37 45 73 77])  = -1/3;
ident.threeD.fourth.I4_dev([11 13 21 25 29 31 51 53 57 61 69 71])  = 1/2;