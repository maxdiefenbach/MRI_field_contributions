function RDF_ppm = forwardSimulate_RDF_ppm(DataParams)
    
    voxelSize_mm = DataParams.voxelSize_mm(:)';
    B0dir = DataParams.B0dir;
    chimap_ppm = DataParams.chimap_ppm;
    
    matrixSize = size(chimap_ppm);
    padsize = ceil((matrixSize + 1) / 2);
    kernelSize = matrixSize + 2 * padsize;
    FOV_mm = kernelSize .* voxelSize_mm;
    DCoffset = 0;
    
    D = get_dipoleKernel_kspace(kernelSize, FOV_mm, B0dir, DCoffset);

    chiBig = padarray(chimap_ppm, padsize, 'replicate');
    Chi = fftshift(fftn(ifftshift(chiBig)));
    
    Psi = D .* Chi;
    
    psiBig = fftshift(ifftn(ifftshift(Psi)));
    psi = depad_array3d(psiBig, padsize);
    
    RDF_ppm = real(psi);

end
