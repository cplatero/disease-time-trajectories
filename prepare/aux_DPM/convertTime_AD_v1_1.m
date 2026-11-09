function [mask_AD,Converts,convertTime,onset]= convertTime_AD_v1_1(ADNIMERGE)



mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');

mask_AD=(ADNIMERGE.DX_bl=='AD') & mask_DX;

%AD: 3474v 542s: 1920v 261s_MCI 1554v 281s_pMCI
ADNI_AD=ADNIMERGE(mask_AD,:); 
%RID=unique(ADNIMERGE.RID); %395s 1192v

% mask_AD=(ADNIMERGE.DX_bl=='AD') ;
% ADNI_AD=ADNIMERGE(mask_AD,:); 

%% clinical regresion
Converts=false(size(ADNIMERGE,1),1);

mask_regresion= ADNI_AD.DX=='MCI' | ADNI_AD.DX=='CN';
rSubjs=unique(ADNI_AD.RID(mask_regresion));
for i=1:length(rSubjs)
    index=find(ADNI_AD.RID==rSubjs(i));
    DX=ADNI_AD.DX(index);
    if(DX(end)=='MCI' || DX(end)=='CN')
        index=find(ADNIMERGE.RID==rSubjs(i));
        Converts(index)=true;
    end
end


mask_AD(Converts)= 0;
convertTime=zeros(sum(mask_AD),1);
Converts=convertTime;
onset=convertTime;

%% sDementia
ADNIMERGE=ADNIMERGE(mask_AD,:); 
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
    time_zero=ADNIMERGE.EXAMDATE(index(ord(1)));
    onset(index)=days(ADNIMERGE.EXAMDATE(index)-time_zero)/30.4167;    


end

% ADNI_rAD= ADNIMERGE(Converts>0,:); %AD: 8v 2s
% ADNI_AD= ADNIMERGE(Converts==0,:); %AD: 752v 231s





end


