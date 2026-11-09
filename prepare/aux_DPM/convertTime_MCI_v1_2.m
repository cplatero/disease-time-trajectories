function [mask_MCI,Converts,convertTime,onset]= convertTime_MCI_v1_2(ADNIMERGE)

mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');

% mask_MCI= ADNIMERGE.DX_bl=='LMCI' & mask_DX;
% mask_MCI= (ADNIMERGE.DX_bl=='EMCI') & mask_DX;
mask_MCI=((ADNIMERGE.DX_bl=='LMCI') | (ADNIMERGE.DX_bl=='EMCI')) & mask_DX;

%LMCI: 3474v 542s: 1920v 261s_MCI 1554v 281s_pMCI
%eMCI: v s:  v 262s_MCI v 46s_pMCI
%elMCI: v s: v 523s_MCI v 327s_pMCI
ADNI_all=ADNIMERGE;
ADNIMERGE=ADNIMERGE(mask_MCI,:); 
% RID=unique(ADNIMERGE.RID); % 1015s 5694v 

% ADNI_all=ADNIMERGE;
% mask_MCI=(ADNIMERGE.DX_bl=='LMCI') | (ADNIMERGE.DX_bl=='EMCI');
% ADNIMERGE=ADNIMERGE(mask_MCI,:); 

%% Delete MCI regresion to CN
[~, ~, ridGroups] = unique(ADNIMERGE.RID); % Agrupar por RID
maxMonthIdx = accumarray(ridGroups, ADNIMERGE.M, [], @max); % Obtener mes máximo por grupo
isLastVisit = ismember([ADNIMERGE.RID, ADNIMERGE.M], [ADNIMERGE.RID, maxMonthIdx(ridGroups)], 'rows'); % Índices donde M es máximo
lastVisitIdx = find(isLastVisit); % Índices finales
last_MCI=ADNIMERGE(lastVisitIdx,:);
mask_del=last_MCI.DX=='CN';
RID_del=last_MCI.RID(mask_del);
for i=1:length(RID_del)
    ADNIMERGE(ADNIMERGE.RID==RID_del(i),:)=[];
    mask_MCI(ADNI_all.RID==RID_del(i))=0;
end



%% Clinical converts
mask_pMCI= (ADNIMERGE.DX=='Dementia'); % people who convert to AD
Converts=false(size(ADNIMERGE,1),1);
pMCI_subjs=unique(ADNIMERGE.RID(mask_pMCI));
for i=1:length(pMCI_subjs)
    index=find(ADNIMERGE.RID==pMCI_subjs(i));
    DX=ADNIMERGE.DX(index);
    if(DX(end)=='Dementia')
        Converts(index)=true;
    end
end

% ADNI_pMCI= ADNIMERGE(Converts>0,:); %pMCI: 1768v 256s
% ADNI_sMCI= ADNIMERGE(Converts==0,:); %sMCI: 2731v 444s



%% pMCI 345s 2304v
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
    ct=find(DX=='Dementia',1);
    convertTime(index)=viscode(ct);
    time_zero=ADNIMERGE.EXAMDATE(index(ord(ct)));
    if(ct>1)
        Delta_days=days(ADNIMERGE.EXAMDATE(index(ord(ct)))-...
            ADNIMERGE.EXAMDATE(index(ord(ct-1))))/2;
    else
        Delta_days=0;      
    end
    onset(index)=(days(ADNIMERGE.EXAMDATE(index)-time_zero)+Delta_days)/30.4167;
    

end



%% sMCI
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


