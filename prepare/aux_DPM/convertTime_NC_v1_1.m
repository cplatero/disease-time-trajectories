function [mask_NC,Converts,convertTime,onset]=convertTime_NC_v1_1(ADNIMERGE)

 
mask_DX=(ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI') | (ADNIMERGE.DX=='CN');
mask_NC=(ADNIMERGE.DX_bl=='CN') & mask_DX; %NC
% mask_NC=(ADNIMERGE.DX_bl=='CN'); %NC
ADNIMERGE=ADNIMERGE(mask_NC,:); %CU:2952v 515s, NC: 2771v 409s, (old)4233v 512s

%% Clinical converts
mask_pNC= (ADNIMERGE.DX=='Dementia') | (ADNIMERGE.DX=='MCI'); % people who convert to MCI/AD
Converts=false(size(ADNIMERGE,1),1);
pNC_subjs=unique(ADNIMERGE.RID(mask_pNC));
for i=1:length(pNC_subjs)
    index=find(ADNIMERGE.RID==pNC_subjs(i));
    DX=ADNIMERGE.DX(index);
    if((DX(end)=='Dementia') || (DX(end)=='MCI'))
        Converts(index)=true;
    end 
end

% ADNIMERGE.Convert=Converts;
% ADNI_pNC= ADNIMERGE(Converts,:); %pNC: 721v 81s
% ADNI_sNC= ADNIMERGE(Converts==0,:); %sNC: 1638v 255s



%% pCU (clinical) 810v 96s
RID=unique(ADNIMERGE.RID(Converts==1));
convertTime=zeros(size(ADNIMERGE,1),1);
onset=convertTime;
for i=1:length(RID)
    index=find(ADNIMERGE.RID==RID(i));
%     if(RID(i)==926)
%         fprintf('Error\n');       
%     end        
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



%% sCU (clinical) 419s 2142v
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


