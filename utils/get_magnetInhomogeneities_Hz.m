function magnetInhomogeneities_Hz = get_magnetInhomogeneities_Hz(B0params, matrixSize, transform)
% magnetInhomogeneities_Hz = get_magnetInhomogeneities_Hz(B0params, matrixSize, transform)
% 
% calculate the inhomogeneities of main magnetic field in scanner xyz
% coordinates from the B0 parameters and
% transform it to the left-handed (flipped z-axis) scanner system
% 
% Input:
%        - B0params: struct to describe the magnet inhomogeneities.
%                    dB0 = sum C[n,m] * Yc[n,m]
%                          n,m               
% 
%                    with the coefficients C[n,m] in [uT] stored in 
%          B0params.B0coeffs_uT with the following structure.
%                    1: C[2,0]   2: C[4,0]   3: C[6,0]
%                    4: C[8,0]   5: C[10,0]  6: C[12,0]
%                    7: C[14,0]  8: C[16,0]  9: C[18,0]
%                    10: C[20,0] ...
%                    (For even n or m, these coefficients are negligible due to symmetry.
%                    The maximum degree n_max=30.)
%          
%                    and Yc[n,m] being the solid harmonic functions
%                                                     n
%                    Yc[n,m](r,phi,theta) = ( r / r0 )  * P[n,m]( cos( theta ) ) * cos( m * phi )
%                    with:       r<=r0 ;  m=0,1,2,... ;  n=m, m+1, m+2,...
%                    where P[n,m]( x ) is the Legendre-polynomial of degree n and order m.
%                    The reference radius r0 in [mm] is stored in
%          B0params.B0refRadius_mm: 
%
%          B0params.B0coeffs_uT: B0 field coefficients C[n,m] [uT]
%        - matrixSize: 1d array of 3 integers, size of the returned array
%        - transform: 4d matrix representation of an affine transformation from
%                     the image coordinate system (REC in MRecon) to 
%                     the left-handed scanner xyz system
% Output:
%        - magnetInhomogeneities_Hz: 3d array of size matrixSize
    
    B0rr_mm = B0params.B0refRadius_mm;
    B0coeffs_uT = B0params.B0coeffs_uT;

    % ijk coordinates
    [I, J, K] = ndgrid(1:matrixSize(1), 1:matrixSize(2), 1:matrixSize(3)); % [pix]

    % transfrom to xyz coordinates
    CSijk = [I(:), J(:), K(:), ones(size(I(:)))]';
    CSxyz = transform * CSijk;

    X = reshape(CSxyz(1, :), matrixSize); % [mm]
    Y = reshape(CSxyz(2, :), matrixSize); % [mm]
    Z = reshape(CSxyz(3, :), matrixSize); % [mm]
    
    R = sqrt(X.^2 + Y.^2 + Z.^2);       % [mm]
    CosTheta = Z ./ R;                  % []

    magnetInhomogeneities_uT = zeros(size(X));
    
    for i = 1:numel(B0coeffs_uT)
        
        if B0coeffs_uT(i) == 0
            continue;
        end

        % spherical harmonics
        n = 2 * i;
        m = 0;
        Yc = legendrePnm(n, m, CosTheta);
        magnetInhomogeneities_uT = magnetInhomogeneities_uT + ...
            B0coeffs_uT(i) * (R/B0rr_mm).^n .* Yc;

    end
    
    gammabar = 42.58e6;                 % gyromagnetic ratio / 2 pi [Hz/T]
    magnetInhomogeneities_Hz = gammabar * 1e-6 * magnetInhomogeneities_uT;

end