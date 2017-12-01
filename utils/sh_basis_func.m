function mat = sh_basis_func( k, X, Y, Z)
% Returns values of the non-normalized spherical harmonic polYnomials that
% describe the shim coil fields.
% |----------+-------------------------------------+---+----|
% | number k | formula                             | l |  m |
% |----------+-------------------------------------+---+----|
% |      0   | 1                                   |   |    |
% |      1   | z                                   | 1 |  0 |
% |      2   | x                                   | 1 |  1 |
% |      3   | y                                   | 1 | -1 |
% |      4   | z * z - 0.5 * (x * x + y * y)       | 2 |  0 |
% |      5   | z * x                               | 2 |  1 |
% |      6   | z * y                               | 2 | -1 |
% |      7   | x * x - y * y                       | 2 |  2 |
% |      8   | 2 * x * y                           | 2 | -2 |
% |      9   | z * (z * z - 1.5 * (x * x + y * y)) | 3 |  0 |
% |     10   | x * (4 * z * z - x * x - y * y)     | 3 |  1 |
% |     11   | y * (4 * z * z - x * x - y * y)     | 3 | -1 |
% |     12   | z * ( x * x - y * y)                | 3 |  2 |
% |     13   | 2 * x * y * z                       | 3 | -2 |
% |     14   | x * x * x - 3 * x * y * y           | 3 |  3 |
% |     15   | 3 * x * x * y - y * y * y           | 3 | -3 |
% |----------+-------------------------------------+---+----|
    
    mat = zeros(size(X));

    switch k
      case 1
	mat = Z;
      case 2
	mat = X;
      case 3
	mat = Y;
      case 4
	mat = Z .* Z - 0.5 .* (X .* X + Y .* Y);
      case 5
	mat = Z .* X;
      case 6
	mat = Z .* Y;
      case 7
	mat = X .* X - Y .* Y;
      case 8
	mat = 2 .* X .* Y;
      case 9
        mat = Z .* (Z .* Z - 1.5 .* (X .* X + Y .* Y));
      case 10
        mat = X .* (4 .* Z .* Z - X .* X - Y .* Y);
      case 11
        mat = Y .* (4 .* Z .* Z - X .* X - Y .* Y);
      case 12
        mat = Z .* ( X .* X - Y .* Y);
      case 13
        mat = 2 .* X .* Y .* Z;
      case 14
        mat = X .* X .* X - 3 .* X .* Y .* Y;
      case 15
        mat = 3 .* X .* X .* Y - Y .* Y .* Y;
    end

end