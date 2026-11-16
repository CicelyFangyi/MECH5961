function [frequency,binPower,CSM]=block_spectra(data,fs,L,frequencies,windowName)
% data: samples by channels, in documented input units. No calibration assumed.
% One-sided power per FFT bin, energy-normalized window, mean of block powers.
% Rectangular window exactly matches current supplied CBF.m for even L.
assert(isreal(data) && all(isfinite(data),'all'));
assert(fs>0 && L>=2 && mod(L,1)==0 && mod(size(data,1),L)==0);
if strcmp(windowName,'rectangular')
    w=ones(L,1);
elseif strcmp(windowName,'hann-periodic')
    w=0.5-0.5*cos(2*pi*(0:L-1)'/L);
else
    error('Unknown window');
end
B=size(data,1)/L; M=size(data,2); K=floor(L/2)+1;
frequency=(0:K-1)'*fs/L;
idx=round(frequencies*L/fs)+1;
assert(all(idx>=1 & idx<=K));
assert(all(abs(frequency(idx)-frequencies(:))<1e-8*fs),'Requested frequency off FFT bin');
scale=sqrt(2)*ones(K,1); scale(1)=1;
if mod(L,2)==0, scale(end)=1; end
pf=complex(zeros(M,B,numel(idx))); binPower=zeros(K,M);
for m=1:M
    F=fft(reshape(double(data(:,m)),L,B).*w,[],1)/sqrt(L*sum(w.^2));
    F=F(1:K,:).*scale;
    binPower(:,m)=mean(abs(F).^2,2);
    for k=1:numel(idx), pf(m,:,k)=F(idx(k),:); end
end
CSM=complex(zeros(M,M,numel(idx)));
for k=1:numel(idx), CSM(:,:,k)=pf(:,:,k)*pf(:,:,k)'/B; end
end
