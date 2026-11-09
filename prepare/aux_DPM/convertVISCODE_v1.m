function VISNUM = convertVISCODE_v1(VISCODE)
VISCODE = upper(string(VISCODE));
viscode2num = @(v) arrayfun(@(x) viscodeTOnum(x), v);
VISNUM = viscode2num(VISCODE);
end

function num = viscodeTOnum(v)
    if ismissing(v) || strlength(v) == 0
        num = NaN;
        return;
    end
    if strcmp(v, "BL")
        num = 0;
    elseif startsWith(v, "M")
        % Extrae los números después de la M
        num = str2double(extractAfter(v, "M"));
        if isnan(num), num = NaN; end
    elseif strcmp(v, "SCMRI")
        num = 0;
    elseif strcmp(v, "SC")
        num = 0;        
    elseif startsWith(v, "Y")
        % Extrae los números después de la M
        num = str2double(extractAfter(v, "Y"))*12;
        if isnan(num), num = NaN; end
        
    else
        num = NaN; % Por si hay otros códigos
    end
end