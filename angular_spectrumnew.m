function img= angular_spectrumnew( dx, r, obj, z )
%angular spectrum propagation by linear convolution
%dx:pixel size r:wavelength obj:field before propagation; 
%z:propagation distance; img:field after propagation

dy = dx;
du = dx;
[nn,mm] = size(abs(obj));

dfx = 1/(2*dx*nn);
dfy = 1/(2*dy*mm);


% % Q operator with dfx, dfy at frequency plane
pha = zeros(2*nn,2*mm);
for ii = 1:2*nn
    for jj = 1:2*mm
        pha(ii,jj) = dfx^2*(ii-nn-0.5)^2 + dfy^2*(jj-mm-0.5)^2; % fx^2 + fy^2
    end
end

% pha = e_pha_dfx;
e_pha = exp(1i*2*pi*z/r.*sqrt(1-r^2.*pha));

objnew=zeros(2*nn,2*mm);
objnew(round(nn/2.0+1):round(3.0*nn/2),round(mm/2.0+1):round(3.0*mm/2))=obj; %zero padding

tmp = fftshift(fft2(fftshift(objnew)));
tmp = tmp.*e_pha;
imgnew = fftshift(ifft2(fftshift(tmp)));
img=imgnew(round(nn/2.0+1):round(3.0*nn/2),round(mm/2.0+1):round(3.0*mm/2));%crop the center part
return