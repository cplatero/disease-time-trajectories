function ADNI = prepare_ADNI_Amyloid(group, biomarkerType, bEstimatedBaseline, bAlternative,bLong)
%% ADNI data
bDelete_M6_M18=false;
[ADNI,~,~,~,~] = getDataADNIMERGE_v1_3b('CAM',bDelete_M6_M18);
if(bLong)
    ADNI.convertTime_bl = putConvertTime_bl(ADNI);
    ADNI=ADNI(ADNI.convertTime_bl>0,:);
end

% load('./data/ADNIMERGE_250321.mat','ADNI');
% [ADNI,~]=sortrows(ADNI,{'RID','M'});

ADNI=sortrows(ADNI,{'RID','M'});
%% WU
WU_cognition=readtable('../data/uwnpsychsum.xlsx');
WU_cognition=WU_cognition(:,3:end);
WU_cognition.M = convertVISCODE_v1(WU_cognition.VISCODE);
WU_cognition=sortrows(WU_cognition,{'RID','M'});
% bShow=true;
% WU_cognition = deleteDuplicates_VISCODE(WU_cognition,bShow);

ADNI = outerjoin(ADNI, WU_cognition, 'Keys', {'RID','M'}, ...
    'MergeKeys', true, 'Type', 'left');

%% Carga de Biomarcadores
if strcmpi(biomarkerType, 'CSF')
    if(bAlternative)
        tlb_CSF = readtable('../data/UPENNBIOMK_ROCHE_ELECSYS_Jun2025.csv');
        tlb_CSF.Properties.VariableNames{'VISCODE2'} = 'VISCODE';
        tlb_CSF.VISCODE = categorical(tlb_CSF.VISCODE);
        tlb_CSF.M = convertVISCODE_v1(tlb_CSF.VISCODE);
        tlb_CSF = sortrows(tlb_CSF,{'RID','M'});
        ADNI = outerjoin(ADNI, tlb_CSF, 'Keys', {'RID','M'}, 'MergeKeys', true, 'Type', 'left');
        ADNI.ptauAB = ADNI.PTAU_tlb_CSF ./ ADNI.ABETA42;
        if(bEstimatedBaseline)
            ADNI.ptauAB_bl=getMarker_bl(ADNI.ptauAB,ADNI.RID,ADNI.M);        
        end
        
    else
        ADNI.ptauAB = ADNI.PTAU_bl./ ADNI.ABETA_bl;
        ADNI.ptauAB_bl=getMarker_bl(ADNI.ptauAB,ADNI.RID,ADNI.M);
    end    
    
elseif strcmpi(biomarkerType, 'PET')
    %load('../data/PET_Amyloid_ADNI.mat', 'tbl_PET');
    load('../data/PET_Amyloid_ADNI_260224.mat', 'tbl_PET');
    ADNI = outerjoin(ADNI, tbl_PET, 'Keys', {'RID','M'}, 'MergeKeys', true, 'Type', 'left');
    if(bEstimatedBaseline)
        ADNI.AMYLOID_STATUS_bl=getMarker_bl(ADNI.AMYLOID_STATUS,ADNI.RID,ADNI.M);
        ADNI.CENTILOIDS_bl=getMarker_bl(ADNI.CENTILOIDS,ADNI.RID,ADNI.M);
    end
end

%% WMH
WMH=readtable('../data/ucd_wmh_data.xlsx');
WMH.VISCODE=categorical(WMH.VISCODE);

% new version
WMH=sortrows(WMH,{'RID','EXAMDATE'});
WMH.VISCODE=categorical(WMH.VISCODE);
WMH.M = convertVISCODE_v1(WMH.VISCODE);
WMH=WMH(~isnan(WMH.M),:);
WMH=sortrows(WMH,{'RID','M'});
bShow=false;
WMH = deleteDuplicates_VISCODE(WMH,bShow);

ADNI = outerjoin(ADNI, WMH, 'Keys', {'RID','M'}, ...
'MergeKeys', true, 'Type', 'left');


ADNI.WMH_bl=getMarker_bl(ADNI.TOTAL_WMH,ADNI.RID,ADNI.M);

%% Selected group
if(group==1)  
    mask_sMCI = (ADNI.DX_bl == 'EMCI' | ADNI.DX_bl == 'LMCI') & ADNI.Convert == 0;
    mask_cog = ADNI.DX_bl == 'CN' | ADNI.DX_bl == 'SMC' | mask_sMCI;
elseif(group==2)
    mask_cog = ADNI.DX_bl == 'EMCI' | ADNI.DX_bl == 'LMCI' | ADNI.DX_bl == 'AD';
else

end

ADNI = ADNI(mask_cog,:);
end



function convertTime_bl = putConvertTime_bl(ADNI)
convertTime_bl=nan(size(ADNI,1),1);
converTime_p=-ADNI.onset/12;

RID=unique(ADNI.RID);
for i=1:length(RID)
    mask_RID= ADNI.RID==RID(i);
    converTime_bl_subj=converTime_p(mask_RID);
    convertTime_bl(mask_RID)=converTime_bl_subj(1);
end
end
