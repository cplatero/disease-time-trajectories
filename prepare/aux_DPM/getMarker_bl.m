function marker_bl=getMarker_bl(marker,RID,M)
% Inicializar la nueva columna
T=table;
T.RID=RID;
T.M=M;
T.marker=marker;
marker_bl = nan(height(T),1);

% Obtener lista de pacientes
rids = unique(T.RID);

for i = 1:numel(rids)
    idx = T.RID == rids(i);        % Filas del paciente
    Tp = T(idx,:);                 % Subtabla del paciente
    
    % Ordenar por visita
    Tp = sortrows(Tp, 'M');
    
    % Encontrar primera visita con PLASMA no NaN
    first_valid = find(~isnan(Tp.marker), 1, 'first');
    
    if ~isempty(first_valid)
        plasma_bl = Tp.marker(first_valid);
        marker_bl(idx) = plasma_bl;
    end
end
