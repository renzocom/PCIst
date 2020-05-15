%PCI State Transitions (Comolatti et al, Brain Stimulation 2019)
%
%function [pci,dNST,parameters]=PCIst(signal_evk,times,parameters)
%
%Please cite this paper if you use this code:
%Comolatti R et al., "A fast and general method to empirically estimate the
%complexity of brain responses to transcranial and intracranial 
%stimulations" Brain Stimulation
%https://doi.org/10.1016/j.brs.2019.05.013
%
% INPUTS: 
% signal_evk = channels x samples
% times = 1 x samples (in milliseconds)
% parameters = struct (optional) - for more details see nested function
%               "checkparameters" 
% OUTPUTS:
% pci = PCIst value
% dNST = PCIst decomposition
%
% [Authors: Renzo Comolatti and Adenauer G. Casali
% Tested on MatlabR2015a, R2016b, R2017a and R2018a
% Last update: 12-may-2020]

function [pci,dNST,parameters]=PCIst(signal_evk,times,parameters)

if nargin<3
    parameters=[];
end
parameters=checkparameters(parameters);
[signal_svd]=dimensionality_reduction(signal_evk,times,parameters);
if size(signal_svd,1)==0
     warning('No components --> PCIst=0');
     pci=0;
     dNST=[];
     return
else
     dNST=statetransitions(signal_svd,times,parameters);
     pci=sum(dNST);
end

 
function parameters=checkparameters(parameters)

stdparameters.max_var=99; %Percentage of variance accounted for by the selected principal components.
stdparameters.min_snr=1.1; %Selects principal components with a signal-to-noise ratio (SNR) > min_snr.
stdparameters.k=1.2;%Noise control parameter.
stdparameters.baseline=[-400 -50];%Signal's baseline time interval [ini,end] in milliseconds.
stdparameters.response=[0 300];%Signal's response time interval [ini,end] in milliseconds.
stdparameters.nsteps=100;% Number of steps used to search for the threshold that maximizes ∆NST.
stdparameters.l=1;%Number of embedding dimensions (1 = no embedding).
stdparameters.tau=2;%Number of timesamples of embedding delay

flds=fieldnames(stdparameters);
for i=1:numel(flds)
    if ~isfield(parameters,flds{i})
        parameters.(flds{i})=stdparameters.(flds{i});
    end
end
flds=fieldnames(parameters);
for i=1:numel(flds)
    if ~isfield(stdparameters,flds{i})
        warning(['Parameters: field ''',flds{i},''' ignored!'])
        parameters=rmfield(parameters,flds{i});
    end
end

 

function [signal_svd,eigenvalues]=dimensionality_reduction(signal,times,parameters)

inds=(times>=parameters.response(1) & times<=parameters.response(2));
indsbase=(times>=parameters.baseline(1) & times<=parameters.baseline(2));
[U,S]=svd(signal(:,inds));
eigenvalues=diag(S);
PCs=U'*signal;
vars=cumsum(eigenvalues.^2);
vars=100.*vars./vars(end);
max_dim=find(vars>=parameters.max_var,1,'first');
signal_svd=PCs(1:max_dim,:);
eigenvalues=eigenvalues(1:max_dim);
if parameters.min_snr>0
    snr=sqrt(mean(signal_svd(:,inds).^2,2)./mean(signal_svd(:,indsbase).^2,2));
    x=find(snr>parameters.min_snr);
    signal_svd=signal_svd(x,:);
    eigenvalues=eigenvalues(x);
end


function [dNST,max_thresholds,D_base,D_resp]=statetransitions(signal,times,parameters)
   
    [vec,td]=embedding(signal,parameters.l,parameters.tau,times);
    inds=find(td>=parameters.response(1) & td<=parameters.response(2));
    indsbase=find(td>=parameters.baseline(1) & td<=parameters.baseline(2));
    D_base=zeros(size(signal,1),numel(indsbase),numel(indsbase));
    D_resp=zeros(size(signal,1),numel(inds),numel(inds));
    for i=1:size(signal,1)
        D_base(i,:,:)=distancia(vec(:,indsbase,i));
        D_resp(i,:,:)=distancia(vec(:,inds,i));
    end
    thresholds=zeros(parameters.nsteps,size(signal,1));
    N_base=zeros(parameters.nsteps,size(signal,1));
    N_resp=zeros(parameters.nsteps,size(signal,1));
    for i=1:size(signal,1)
        minthr=median(reshape(D_base(i,:,:),1,size(D_base,2)*size(D_base,2)));
        maxthr=max(reshape(D_resp(i,:,:),1,size(D_resp,2)*size(D_resp,2)));
        thresholds(:,i)=linspace(minthr,maxthr,parameters.nsteps);
        for j=1:size(thresholds,1)
            RP=zeros(size(D_base,2),size(D_base,2));
            RP(D_base(i,:,:)<=thresholds(j,i))=1;
            T_base=abs(diff(RP,1,2));
            RP=zeros(size(D_resp,2),size(D_resp,2));
            RP(D_resp(i,:,:)<=thresholds(j,i))=1;
            T_resp=abs(diff(RP,1,2));
            N_base(j,i)=sum(sum(T_base))/(size(D_base,2)^2);
            N_resp(j,i)=sum(sum(T_resp))/(size(D_resp,2)^2);
        end
    end
    NST_diff=N_resp-parameters.k*N_base;
    [~,b]=max(NST_diff);
    max_thresholds=diag(thresholds(b,:));
    dNST=diag(NST_diff(b,:))'.*size(D_resp,2);
    


function [vec,td]=embedding(Y,L,tau,t)
    N=size(Y,2);
    M=(L-1)*tau+1;
    td=t(M:end);
    vec=zeros(L,N-M+1,size(Y,1));
    if L==1
        vec(1,:,:)=Y';
    else
        for j=1:size(Y,1)
            for k=1:L
                vec(k,:,j)=Y(j,M-(k-1)*tau:1:N-(k-1)*tau);
            end
        end
    end

    
function d=distancia(vec)
    d=zeros(size(vec,2),size(vec,2));
    for j=1:size(vec,2)
        d(j,:)=sqrt(sum((vec(:,:)-repmat(vec(:,j),1,size(vec,2))).^2,1));
    end
        
