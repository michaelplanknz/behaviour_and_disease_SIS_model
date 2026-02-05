function Ytrunc = trunc01(Y)

% Keep only rows of the matrix Y where all elements in the row are in [0,1]

inFlag = (Y >= 0 & Y <= 1);
keepFlag = all(inFlag, 2);
Ytrunc = Y(keepFlag, :);


