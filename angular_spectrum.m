function img= angular_spectrum( dx, r, obj, z )
%angular spectrum propagation by circular convolution
%dx:pixel size r:wavelength obj:field before propagation; 
%z:propagation distance; img:field after propagation
dy = dx;
du = dx;
[nn,mm] = size(abs(obj));

dfx = 1/(dx*nn);
dfy = 1/(dy*mm);

% % Q operator with dfx, dfy at frequency plane
pha = zeros(nn,mm);
for ii = 1:nn
    for jj = 1:mm
        pha(ii,jj) = dfx^2*(ii-nn/2-0.5)^2 + dfy^2*(jj-mm/2-0.5)^2; % fx^2 + fy^2
    end
end

% pha = e_pha_dfx;
e_pha = exp(1i*2*pi*z/r.*sqrt(1-r^2.*pha));

%shibu=real(e_pha);
%figure; imshow(abs(shibu));

tmp = fftshift(fft2(fftshift(obj)));
tmp = tmp.*e_pha;
img = fftshift(ifft2(fftshift(tmp)));
%img=img/(4*mm*nn);

return