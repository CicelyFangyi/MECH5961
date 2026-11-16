function [power,level]=cbf_map(CSM,xy,x,y,h,f,c0)
% Supplied CBF.m Green function and unit-response weights, unchanged model.
% level is normalized beamforming power, NOT microphone SPL.
assert(isequal(size(CSM),[size(xy,1),size(xy,1)]));
assert(h>0 && f>0 && c0>0 && all(isfinite(CSM),'all'));
[X,Y]=meshgrid(x,y); mapSize=size(X); X=X(:); Y=Y(:); power=zeros(numel(X),1);
for first=1:8192:numel(X)
    q=first:min(first+8191,numel(X));
    D=sqrt((xy(:,1)-X(q)').^2+(xy(:,2)-Y(q)').^2+h^2);
    T=exp(-1i*2*pi*f*D/c0)./(4*pi*D);
    w=T./sum(abs(T).^2,1);
    power(q)=real(sum(conj(w).*(CSM*w),1)).';
end
assert(min(power)>=-1e-10*max(1,max(abs(power))),'Negative beamforming power');
power=reshape(max(power,0),mapSize);
if max(power,[],'all')==0
    level=nan(size(power));
else
    level=10*log10(max(power,realmin)/max(power,[],'all'));
end
end
