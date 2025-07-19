    A_fft = fftshift(fft(ifftshift(A_in)))
    A_fft *= dispersion_operator
    A_lin = fftshift(ifft(ifftshift(A_fft)))

    A_nl = A_lin * np.exp(-1j * gammaconstant * getPower(A_lin) * dz - (alpha_dB_per_m / 2) * dz)

    A_fft = fftshift(fft(ifftshift(A_nl)))
    A_fft *= dispersion_operator
    A_out = fftshift(ifft(ifftshift(A_fft)))
    return A_out