function [Tresults,features,features_names,ni,features_SNan] = getDataADNIMERGE_v1_3b(feat_sel,bDelete_M6_M18)    
    % This function is charged of importing and cleaning the Tresults
    % table, taking only candidate features. However, Tresults is equally
    % returned so cleaned variables are still accessible. Extracted
    % features can be Rois, neuropsycological tests, clusters
    % or all of them.
    %
    %
    % possible values of feat_sel:
    %       'ROI'
    %       'CAM'
    %       'MASS'
    %       'all'
    
    %% Data input
    %load('./data/ADNIMERGE_220706','ADNI');
    load('../data/ADNIMERGE_250321','ADNI');
    %ADNI(ADNI.RID==6505,:)=[];%for NaN in AGE
    ADNI=sortrows(ADNI,{'RID','M'});
    %ADNI(ADNI.M==3,:)=[];
    if(bDelete_M6_M18)
        ADNI(ADNI.M==6,:)=[];
        ADNI(ADNI.M==18,:)=[];
    end
    
    Tresults=ADNI;

    ni=get_ni(Tresults.M,length(unique(Tresults.RID)));
   
    %% MRI ROI markers
    %Combination of both hemispheres (volumes)
    NHV=Tresults.Hippocampus./Tresults.ICV*1e3;
    NVV=Tresults.Ventricles./Tresults.ICV*1e3;
    NFV=Tresults.Fusiform./Tresults.ICV*1e3;
    NMV=Tresults.MidTemp./Tresults.ICV*1e3;
    NEV=Tresults.Entorhinal./Tresults.ICV*1e3;
    
    feat_Roi=[NHV, NVV, NFV, NMV, NEV];
    
    feat_Roi_names={'NHV','NVV','NFV','NMV','NEV'};
    feat_Roi_SNan = clean_invalid_CAM(feat_Roi,ni,Tresults.M);
    
    
    %clear variables for debugging
    clear NHV NVV NFV NMV NEV;
    
    
    %% CAM MARKERS
%     feat_CAM = [Tresults.RAVLT_forgetting,...%cell2mat(Tresults.RAVLT_forgetting),...
%         Tresults.RAVLT_immediate,...
%         Tresults.RAVLT_learning,...
%         Tresults.RAVLT_perc_forgetting,...%? k paso?
%         %Tresults.ADAS11,...
%         Tresults.ADAS13,...
%         Tresults.FAQ,...
%         Tresults.MMSE,...
%         Tresults.CDRSB,...
%         %Tresults.ADASQ4,...
%         %Tresults.MOCA,...
%         %Tresults.EcogPtTotal,...
%         %Tresults.EcogSPTotal,...
%         Tresults.LDELTOTAL];%,...
%         %Tresults.DIGITSCOR];

    feat_CAM = [Tresults.RAVLT_forgetting,...%cell2mat(Tresults.RAVLT_forgetting),...
        Tresults.RAVLT_immediate,...
        Tresults.RAVLT_learning,...
        Tresults.RAVLT_perc_forgetting,...%? k paso?
        Tresults.ADAS13,...
        Tresults.FAQ,...
        Tresults.MMSE,...
        Tresults.CDRSB,...
        Tresults.LDELTOTAL];
    
%     feat_CAM_names = {'RAVLT_forgetting',...
%         'RAVLT_immediate',...
%         'RAVLT_learning',...
%         'RAVLT_perc_forgetting',...
%         %'ADAS11',...
%         'ADAS13',...
%         'FAQ',...
%         'MMSE',...
%         'CDRSB',...
%         %'ADASQ4',...
%         %'MOCA',....
%         %'EcogPtTotal',...
%         %'EcogSPTotal',...
%         'LDELTOTAL'};%,...
%         %'DIGITSCOR'};
        
    feat_CAM_names = {'RAVLT_forgetting',...
        'RAVLT_immediate',...
        'RAVLT_learning',...
        'RAVLT_perc_forgetting',...
        'ADAS13',...
        'FAQ',...
        'MMSE',...
        'CDRSB',...
        'LDELTOTAL'};
        
    feat_CAM_SNan = clean_invalid_CAM(feat_CAM,ni,Tresults.M);
    
    %% CSF
    if(strcmp(feat_sel,'all') || strcmp(feat_sel,'CAMCSF'))
        feat_CSF = [Tresults.ABETA,...
            Tresults.TAU,...
            Tresults.PTAU,...
            Tresults.TAU./Tresults.ABETA,...
            Tresults.PTAU./Tresults.ABETA];

        feat_CSF_names = {'ABETA',...
            'TAU',...
            'PTAU',...
            'TAU_AB',...
            'PTAU_AB'};
        
        feat_CSF_SNan = clean_invalid_CAM(feat_CSF,ni,Tresults.M);
        
        
    end
    %% glucose metabolism: average normalized FDG uptake 
    feat_FDG=Tresults.FDG;
    feat_FDG_names={'FDG'};
    feat_FDG_SNan = clean_invalid_CAM(feat_FDG,ni,Tresults.M);
    
    %% brain amyloidosys: average normalized AV45 uptake
    feat_AV45=Tresults.AV45;
    feat_AV45_names={'AV45'};
    feat_AV45_SNan = clean_invalid_CAM(feat_AV45,ni,Tresults.M);
    %% PACC
    feat_PACC = [Tresults.mPACCdigit,Tresults.mPACCtrailsB,Tresults.TRABSCOR];
    feat_PACC_names={'PACCdigit','PACCtrailsB','TRABSCOR'};
    feat_PACC_SNan = clean_invalid_CAM(feat_PACC,ni,Tresults.M);
     
    
    %% Cortical thickness clusters
    
    % As number of clusters is undefined, columns are selected thanks to
    % clusters name: ?h_clxx, where ? is r or l and xx is the cluster
    % number.
    
    %feat_Mass_names = original_names(contains(original_names,'h_cl'));
    %feat_Mass = table2array(Tresults(:,feat_Mass_names));
    
    
    %% Appropriate combination
    
    if strcmp(feat_sel,'ROI')
        features = feat_Roi;
        features_SNan = feat_Roi_SNan;        
        features_names = feat_Roi_names';
    elseif strcmp(feat_sel,'CAM')
        features = [feat_CAM, feat_PACC];
        features_SNan = [feat_CAM_SNan, feat_PACC_SNan];        
        features_names = [feat_CAM_names, feat_PACC_names]';
    elseif strcmp(feat_sel,'all')
        features = [feat_Roi, feat_CAM, feat_PACC, feat_CSF, feat_FDG, feat_AV45];
        features_SNan = [feat_Roi_SNan, feat_CAM_SNan,feat_PACC_SNan, feat_CSF_SNan, feat_FDG_SNan, feat_AV45_SNan, ];
        features_names = [feat_Roi_names, feat_CAM_names,feat_PACC_names,feat_CSF_names,feat_FDG_names,feat_AV45_names]';        
     elseif strcmp(feat_sel,'ROICAM')
        features = [feat_Roi, feat_CAM, feat_PACC];
        features_SNan = [feat_Roi_SNan, feat_CAM_SNan, feat_PACC_SNan];
        features_names = [feat_Roi_names, feat_CAM_names, feat_PACC_names]';
    elseif strcmp(feat_sel,'CAMCSF')
        features = [feat_CAM, feat_PACC, feat_CSF];
        features_SNan = [feat_CAM_SNan, feat_PACC_SNan, feat_CSF_SNan];        
        features_names = [feat_CAM_names, feat_PACC_names, feat_CSF_names]';
        
    else
        error('Invalid input statement. feat_sel must be ROI, CAM, MASS or all');
    end

%% CU
[mask_NC,Convert,convertTime,onset]=convertTime_NC_v1_1(Tresults);

Tresults_NC=Tresults(mask_NC,:);
Tresults_NC.Convert=Convert;
Tresults_NC.ConvertTime=convertTime;
Tresults_NC.onset=onset;


features_NC=features(mask_NC,:);
features_SNan_NC=features_SNan(mask_NC,:);

%% MCI
[mask_MCI,Convert,convertTime,onset]= convertTime_MCI_v1_2(Tresults);

Tresults_MCI=Tresults(mask_MCI,:);
Tresults_MCI.Convert=Convert;
Tresults_MCI.ConvertTime=convertTime;
Tresults_MCI.onset=onset;


features_MCI=features(mask_MCI,:);
features_SNan_MCI=features_SNan(mask_MCI,:);

%% AD
[mask_AD,Convert,convertTime,onset]= convertTime_AD_v1_1(Tresults);

Tresults_AD=Tresults(mask_AD,:);
Tresults_AD.Convert=Convert;
Tresults_AD.ConvertTime=convertTime;
Tresults_AD.onset=onset;


features_AD=features(mask_AD,:);
features_SNan_AD=features_SNan(mask_AD,:);

%% SMC
[mask_SMC,Convert,convertTime,onset]=convertTime_SMC_v1_1(Tresults);

Tresults_SMC=Tresults(mask_SMC,:);
Tresults_SMC.Convert=Convert;
Tresults_SMC.ConvertTime=convertTime;
Tresults_SMC.onset=onset;


features_SMC=features(mask_SMC,:);
features_SNan_SMC=features_SNan(mask_SMC,:);

%% Join
Tresults=[Tresults_NC;Tresults_MCI;Tresults_AD;Tresults_SMC];
MCI=[zeros(size(Tresults_NC,1),1);ones(size(Tresults_MCI,1),1);...
    ones(size(Tresults_AD,1),1)*2;ones(size(Tresults_SMC,1),1)*(-1)];
Tresults.stMCI=MCI;
features=[features_NC;features_MCI;features_AD;features_SMC];
features_SNan=[features_SNan_NC;features_SNan_MCI;features_SNan_AD;features_SNan_SMC];
ni=get_ni(Tresults.M,length(unique(Tresults.RID)));

end

