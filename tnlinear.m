function [tnlx, tnly] = tnlinear(ut, vt, utfis, vtfis, nx, ny, kx, ky)
    %----------------------------------------------------------------------
    uu = (utfis .* utfis);
    vv = (vtfis .* vtfis);
    uv = (utfis .* vtfis);

    uu = ifft2(uu);
    uv = ifft2(uv);
    vv = ifft2(vv);
    %----------------------------------------------------------------------
    tnlxa = zeros(nx, ny);
    tnlya = zeros(nx, ny);
    for i = 1:nx
        for j = 1:ny
            tnlxa(i, j) = 1i * kx(i) * uu(i, j) + 1i * ky(j) * uv(i, j);
            tnlya(i, j) = 1i * kx(i) * uv(i, j) + 1i * ky(j) * vv(i, j);
        end
    end
    % ---------------------------------------------------------------------
    dudx = ifft2(1i * kx' .* ut);
    dudy = ifft2(1i * ky .* ut);
    dvdx = ifft2(1i * kx' .* vt);
    dvdy = ifft2(1i * ky .* vt);
    % ---------------------------------------------------------------------
    dudx = real(dudx);
    dudy = real(dudy);
    dvdx = real(dvdx);
    dvdy = real(dvdy);
    % ---------------------------------------------------------------------
    tnlxc = fft2(utfis .* dudx + vtfis .* dudy);
    tnlyc = fft2(utfis .* dvdx + vtfis .* dvdy);
    % ---------------------------------------------------------------------
    tnlx = 0.5 * (tnlxa + tnlxc);
    tnly = 0.5 * (tnlya + tnlyc);
end
