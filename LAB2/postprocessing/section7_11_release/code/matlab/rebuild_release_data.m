function rebuild_release_data(labRoot,reuseValidated)
% Optional raw-data rebuild. No figures; never writes external inputs.
% labRoot contains Beamforming/Coor_56Array.txt and Group2-.../MicData_*.mat.
% Default rebuild requires only project-owned helpers distributed beside this.
if nargin<1,error('Supply authorized local LAB2 input directory.');end
if nargin<2,reuseValidated=false;end
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
out=fullfile(root,'data_frozen'); if ~isfolder(out),mkdir(out);end
xy=readmatrix(fullfile(labRoot,'Beamforming','Coor_56Array.txt'));
assert(isequal(size(xy),[56 2]));
x=linspace(-.5,.5,301); xb=linspace(-.5,.5,151); xf=linspace(-.5,.5,601);
frequencies=1000:250:6000; fs=50000; L=10000; c0=343; h=1;
save(fullfile(out,'geometry.mat'),'xy','x','xb','xf','frequencies','fs','L','c0','h','-v7');
for run=0:3
    target=fullfile(out,sprintf('run%d.mat',run));
    if isfile(target),fprintf('Existing frozen run %d retained; remove explicitly for rebuild.\n',run);continue;end
    D=load(fullfile(labRoot,'Group2-20261116T100220Z-1-001',sprintf('MicData_%d.mat',run)),'ConvertedData');
    md=D.ConvertedData.Data.MeasuredData;
    if run<3, f=frequencies;else,f=[2000 3000 6185];end
    pf=complex(zeros(56,50,numel(f)));
    for ch=1:56
        assert(strcmp(md(ch+2).Name,sprintf('Untitled/Sound Pressure_%d',ch-1)));
        d=double(md(ch+2).Data(:));assert(numel(d)==500000 && all(isfinite(d)));
        F=fft(reshape(d,L,50),[],1)/L;
        pf(ch,:,:)=permute(sqrt(2)*F(round(f/5)+1,:),[3 2 1]);
        if ch==1 && run==3
            [spectrumFrequency,spectrumPower]=block_spectra(d,fs,L,[],'hann-periodic');
        end
    end
    clear D md d F;
    C=zeros(56,56,numel(f)); macroCSM=complex(zeros(56,56,numel(f),5));
    meanMaps=zeros(301,301,numel(f)); macroMaps=zeros(151,151,numel(f),5);
    fineF=f(ismember(f,[1000 3000 6000]));if run==3,fineF=f;end
    fineMaps=zeros(601,601,numel(fineF)); referenceDeviation=zeros(size(fineF));
    for k=1:numel(f)
        A=pf(:,:,k);C(:,:,k)=A*A'/50;
        for b=1:5
            Q=A(:,(b-1)*10+(1:10));macroCSM(:,:,k,b)=Q*Q'/10;
            macroMaps(:,:,k,b)=cbf_map(macroCSM(:,:,k,b),xy,xb,xb,h,f(k),c0);
        end
        meanMaps(:,:,k)=cbf_map(C(:,:,k),xy,x,x,h,f(k),c0);
        j=find(fineF==f(k));
        if ~isempty(j)
            fineMaps(:,:,j)=cbf_map(C(:,:,k),xy,xf,xf,h,f(k),c0);
            old=fullfile(labRoot,'postprocessing','output','cache',sprintf('map_run%d_f%d_n601.mat',run,f(k)));
            if reuseValidated && isfile(old)
                V=load(old); referenceDeviation(j)=max(abs(V.power-fineMaps(:,:,j)),[],'all')/max(V.power,[],'all');
                assert(referenceDeviation(j)<1e-12);
            end
        end
        fprintf('run%d f%d completed\n',run,f(k));
    end
    if run==3
        save(target,'run','f','meanMaps','macroMaps','fineMaps','fineF','C','macroCSM','referenceDeviation','spectrumFrequency','spectrumPower','-v7');
    else
        save(target,'run','f','meanMaps','macroMaps','fineMaps','fineF','C','macroCSM','referenceDeviation','-v7');
    end
end
end
