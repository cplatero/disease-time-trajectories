function [mask_SMC,Converts,convertTime,onset]= convertTime_SMC_v1_1(ADNIMERGE)


 
mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');
mask_SMC=(ADNIMERGE.DX_bl=='SMC') & mask_DX; %NC
ADNIMERGE=ADNIMERGE(mask_SMC,:); %SMC: 456v 105s

% mask_SMC= ADNIMERGE.DX_bl=='SMC';
% ADNIMERGE=ADNIMERGE(mask_SMC,:); 

%% Clinical converts
mask_pSMC= (ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI'); % people who convert to MCI/AD
Converts=false(size(ADNIMERGE,1),1);
pSMC_subjs=unique(ADNIMERGE.RID(mask_pSMC));
for i=1:length(pSMC_subjs)
    index=find(ADNIMERGE.RID==pSMC_subjs(i));
    DX=ADNIMERGE.DX(index);
    if((DX(end)=='Dementia') || (DX(end)=='MCI'))
        Converts(index)=true;
    end 
end




%% pSMC (clinical) 22s 91 v
RID=unique(ADNIMERGE.RID(Converts==1));
convertTime=zeros(size(ADNIMERGE,1),1);
onset=convertTime;
for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
    viscode=ADNIMERGE.M(index);
    [viscode,ord]=sort(viscode);
    if(length(ord)>1)
        if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
            fprintf(2,'No sorting visits from RID %d.\n',RID(i));
        end
    end
    DX=ADNIMERGE.DX(index(ord));
    ct=find((DX=='MCI') | (DX=='Dementia'),1);
    convertTime(index)=viscode(ct);
    time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
    Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
        ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;
    onset(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30.4167;
    
end


%% sSMC (clinical)
RID=unique(ADNIMERGE.RID(Converts==0));
for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
    viscode=ADNIMERGE.M(index);
    [viscode,ord]=sort(viscode);
    if(length(ord)>1)
        if(sum(ord(2:end)-ord(1:end-1)-1)~=0)
            fprintf(2,'No sorting visits from RID %d.\n',RID(i));
        end
    end    
    convertTime(index)=viscode(end);
    time_zero=ADNIMERGE.EXAMDATE(index(ord(end)));
    onset(index)=days(ADNIMERGE.EXAMDATE(index)-time_zero)/30.4167;    
end





end


