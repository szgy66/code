%This code is a simualtion for A multiplexed optical image security system
%with cascaded phase-only masks. Four pairs of input host images and output
%hidden images are tested.
clear;
N=6; %N-1: number of cascaded phase-only masks (L=N-1); the last mask is conjugate(P)
size1=512;
size2=512;
size3=4; %Number of input-output image pairs for image hiding

inputall=zeros(size1,size2,size3);
%Four input host images
inputall(:,:,1)=im2double(imread('p1.bmp'));
inputall(:,:,2)=im2double(imread('p2.bmp'));
inputall(:,:,3)=im2double(imread('p3.bmp'));
inputall(:,:,4)=im2double(imread('p4.bmp'));

%paramaters in simulating the Fresnel field propagation
dist=0.05;%distance between neighboring phase-only masks
lamda=632e-9;%wavelength: 632nm
psize=8e-6;%pixel size: 8 micrometer

phasemask=exp(1i*2*pi*rand(size1,size2,N)); %initial values for the cascaded phase-only masks

targetall=zeros(size1,size2,size3);
%target output results (hidden image displayed in a sub-window in the output imaging plane)
targetall(1:256,1:256,1)=imresize(im2double(imread('p5.bmp')),0.5);
targetall(1:256,257:512,2)=imresize(im2double(imread('p6.bmp')),0.5);
targetall(257:512,1:256,3)=imresize(im2double(imread('p7.bmp')),0.5);
targetall(257:512,257:512,4)=imresize(im2double(imread('p8.bmp')),0.5);

%display the target output results
imwrite(targetall(:,:,1),'target1.bmp','bmp');
imwrite(targetall(:,:,2),'target2.bmp','bmp');
imwrite(targetall(:,:,3),'target3.bmp','bmp');
imwrite(targetall(:,:,4),'target4.bmp','bmp');
temp=targetall(:,:,1)+targetall(:,:,2)+targetall(:,:,3)+targetall(:,:,4);
imwrite(temp,'targetall.bmp','bmp');

%Each input host image is multiplied with a random phase mask
for ii=1:size3    
    inputall(:,:,ii)=inputall(:,:,ii).*exp(1i*2*pi*rand(size1,size2));
end

%Wavefront-matching algorithm for designing the phase-only masks
for iter=1:30 %number of iterations in the optimization
    iter
    for ii=1:N %1:N
        summation=0;        
        for mm=1:size3
            inputpat=inputall(:,:,mm);
            temp1=inputpat;
            temp1=angular_spectrum(psize,lamda,temp1,dist);
            %diffractive field propagation with angular spectrum method (dist can be positive (forward) or negative (backward)) 
            if ii>1
            for kk=1:(ii-1)                
                temp1=temp1.*phasemask(:,:,kk);
                temp1=angular_spectrum(psize,lamda,temp1,dist);
            end
            end
                                   
            outputpat=targetall(:,:,mm);
            temp2=outputpat;
            
            if ii<N
            for kk1=(ii+1):N
                kk=(N+ii+1)-kk1;
                temp2=temp2.*conj(phasemask(:,:,kk));
                temp2=angular_spectrum(psize,lamda,temp2,-dist);                                
            end
            end
            maskcom=temp2.*conj(temp1);
            summation=summation+maskcom;
        end        
        phasemask(:,:,ii)=exp(1i*angle(summation));
    end
end

%save the phase-only masks
save phasemasknew.mat phasemask

%Output result when Host image 1 is input to the system individually
temp1=inputall(:,:,1);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'result1.bmp','bmp');

%Output result when Host image 2 is input to the system individually
temp1=inputall(:,:,2);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'result2.bmp','bmp');

%Output result when Host image 3 is input to the system individually
temp1=inputall(:,:,3);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'result3.bmp','bmp');

%Output result when Host image 4 is input to the system individually
temp1=inputall(:,:,4);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'result4.bmp','bmp');

%Output result when all the four host images are input to the system
temp1=inputall(:,:,1)+inputall(:,:,2)+inputall(:,:,3)+inputall(:,:,4);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'resultall.bmp','bmp');

%Output result when Host image 1 and Host image 3 are input to the system
temp1=inputall(:,:,1)+inputall(:,:,3);
for kk=1:N
    temp1=angular_spectrum(psize,lamda,temp1,dist);
    temp1=temp1.*phasemask(:,:,kk);
end
finalmag=abs(temp1);
vmax=max(max(finalmag));
vmin=min(min(finalmag));
vnorm=(finalmag-vmin)/(vmax-vmin);
imwrite(vnorm,'resultpart.bmp','bmp');



