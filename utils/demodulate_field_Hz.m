function signal = demodulate_field_Hz(signal, TE_s, field_Hz)
    
    for iTE = 1:length(TE_s)
        signal(:, :, :, iTE) = signal(:, :, :, iTE) .* ...
            exp(-1j * 2 * pi * field_Hz * TE_s(iTE));
    end

end