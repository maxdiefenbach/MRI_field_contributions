function objectBasedFieldmap_Hz = get_objectBasedFieldmap_Hz(signal, voxelSize_mm, B0dir, centerFreq_Hz)
% objectBasedFieldmap_Hz = get_objectBasedFieldmap_Hz(signal, voxelSize_mm, B0dir, centerFreq_Hz)
% 
% implementation of the object-based fast field map estimation described in
% Sharma, S. D., Artz, N. S., Hernando, D., Horng, D. E., & Reeder, S. B., 
% Improving chemical shift encoded water-fat separation using object-based 
% information of the magnetic field inhomogeneity, 
% Magnetic Resonance in Medicine, 73(2), 597–604 (2014).  
% http://dx.doi.org/10.1002/mrm.25163


    CHI_H20_ppm = -9.05;
    CHI_FAT_ppm = -7.79;
    CHI_TISSUE_ppm = (CHI_H20_ppm + CHI_FAT_ppm)/2.;  % -8.42ppm
    CHI_AIR_ppm = 0.36;

    magnitude = abs(signal);
    echoMIP = squeeze(sqrt(sum(magnitude.^2, ndims(signal))));
    threshold = 0.05;
    tissueMask = echoMIP >= threshold .* max(echoMIP(:));
    chimap_ppm = ~tissueMask * CHI_AIR_ppm + tissueMask * CHI_TISSUE_ppm;
    
    DataParams.voxelSize_mm = voxelSize_mm;
    DataParams.B0dir = B0dir;
    DataParams.chimap_ppm = chimap_ppm;
    RDF_ppm = forwardSimulate_RDF_ppm(DataParams);
    
    objectBasedFieldmap_Hz = centerFreq_Hz * 1e-6 * (RDF_ppm - mean(RDF_ppm(:)));
end
