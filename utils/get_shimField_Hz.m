function shimField_Hz = get_shimField_Hz(shimValues, matrixSize, transform)
% shimField_Hz = get_shimField_Hz(shimValues, matrixSize, transform)
%
% from shim values calculate shimfield in scanner xyz coordinates and
% transform it to the left-handed (flipped z-axis) scanner system
% 
% Input: 
%        - shimValues:  1d array of 15 doubles
%                       order of shimValues: 
%                       Z, X, Y, Z2, ZX, ZY, X2-Y2, 2XY, Z3, Z2X, Z2Y, ZX2Y2, ZXY, X3, Y3
%                       shimValues have units of [mT/m] or [mT/m^2]
%        - matrixSize: 1d array of 3 integers, size of the returned array
%        - transform:   4d matrix representation of an affine transformation from
%                       the image coordinate system (REC in MRecon) to the 
%                       left-handed scanner xyz system
% 
% Output:
%        - shimfield_T: 3d array of size matrixSize
    
    if sum(shimValues(:)) == 0
        shimfield_T = 0;
        return
    end

    % ijk coordinates
    [I, J, K] = ndgrid(1:matrixSize(1), 1:matrixSize(2), 1:matrixSize(3)); % [pix]

    % transfrom to xyz coordinates
    CSijk = [I(:), J(:), K(:), ones(size(I(:)))]';
    CSxyz = transform * CSijk;

    X = reshape(CSxyz(1, :), matrixSize); % [mm]
    Y = reshape(CSxyz(2, :), matrixSize); % [mm]
    Z = reshape(CSxyz(3, :), matrixSize); % [mm]
    
    % allocate shimfield_T matrix
    shimfield_T = zeros(size(X));
    
    for i = 1:numel(shimValues)
        
        % conversio_factor to convert to [T]
        if i <= 3
            conversio_factor = 1e-6;            % mT/m * mm = 1e-3 T * 1e-3 m / m = 1e-6 T
        elseif i <= 8
            conversio_factor = 1e-9;    % mT/m^2 * mm^2 = 1e-3 T * 1e-6 m^2 / m^2 = 1e-9 T
        elseif i <= 15
            conversio_factor = 1e-12;  % mT/m^3 * mm^3 = 1e-3 T * 1e-9 m^3 / m^3 = 1e-12 T
        end
        
        % spherical harmonics
        sh = sh_basis_func(i, X, Y, Z); % [mm^l]
        
        shimfield_T = shimfield_T + conversio_factor * shimValues(i) * sh;
    end
    
    gammabar = 42.58e6;                 % gyromagnetic ratio / 2 pi [Hz/T]    
    shimField_Hz = gammabar * shimfield_T;

end
