function residualLinearField_Hz = get_residualLinearField_Hz(signal, TE_s, cropDim)
% residualLinearField_Hz = get_residualLinearField_Hz(signal, TE_s, cropDim)
% 
% measure shifts of the k-space maximum along dimension cropDim and
% translate to linear field in image space

    sz = size(signal);
    matrixSize = sz(1:3);
    nTE = sz(4);

    padsize = ceil((matrixSize + 1) / 2);
    padfactor = 10;
    padsize(cropDim) = padfactor * padsize(cropDim);

    delta_k = zeros(3, 3);
    for iTE = 1:nTE
        % zero-pad signal
        s = squeeze(signal(:, :, :, iTE));
        S = padarray(s, padsize);
        
        % Fourier transform
        Sk = ifftshift(fftn(S));
        
        % find k-space maximum
        [kmax, index_kmax] = max(abs(Sk(:)));
        [ix, iy, iz] = ind2sub(size(Sk), index_kmax);
        
        % calculate k-space shifts delta_k in [Hz/pix]
        center = ceil((size(Sk) + 1) / 2);
        delta_k(:, iTE) = ([ix, iy, iz] - center) ./ size(Sk);
    end
    
    % replace first echo with origin of the (TE-delta_k)-plane
    TE_s = [0, TE_s(2:end)];
    delta_k0 = [0, delta_k(cropDim, 2:end)];
    
    % fit delta_k0 vs. TE_s
    p = polyfit(TE_s, delta_k0, 1);
    slope_a = p(1);
    
    % assemble residual linear field
    [I, J, K] = ndgrid(1:matrixSize(1), 1:matrixSize(2), 1:matrixSize(3));
    IJK = cat(4, I, J, K);
    x = IJK(:, :, :, cropDim);
    residualLinearField_Hz = slope_a .* x;

end