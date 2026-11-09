function T_sinDuplicados = deleteDuplicates_TimeCoherence(T,bShow)
% Elimina duplicados de visitas (RID, M) en la tabla T, priorizando la coherencia temporal.
% La visita elegida es la que tiene la EXAMDATE que mejor se ajusta a un patrón 
% de separación de 6 o 12 meses respecto a las otras visitas del mismo paciente (RID).
% 
% T debe tener las columnas: RID (ID paciente), M (Momento/Visita), EXAMDATE (Fecha de la visita).
% Se asume que T ya está ordenada por RID y EXAMDATE.

% 1. Definir los intervalos de coherencia en días
INTERVALS_MONTHS = [6, 12];
% Usamos una aproximación robusta para la duración de los meses
INTERVALS_DAYS = INTERVALS_MONTHS * (365.25 / 12); 

% 2. Identificar Grupos Duplicados (RID, M)
G_RID_M = findgroups(T(:, {'RID', 'M'}));
counts = splitapply(@numel, T.RID, G_RID_M);
duplicados_G = counts > 1;

% Filas de la tabla original T que NO son duplicados (se conservan automáticamente)
T_no_dups = T(ismember(G_RID_M, find(~duplicados_G)), :);
T_duplicadas = T(ismember(G_RID_M, find(duplicados_G)), :);

if isempty(T_duplicadas)
    T_sinDuplicados = T;
    % disp('No se encontraron visitas duplicadas (RID, M). La tabla se devuelve sin cambios.')
    return;
end

if(bShow)
    % Muestra los duplicados únicos
    disp('Combinaciones RID-M duplicadas:')
    unique(T_duplicadas(:, {'RID','M'}))
end


% 3. Identificar pacientes con duplicados y preparar para el análisis
uniqueRIDs = unique(T_duplicadas.RID);
T_dups_chosen = table(); % Tabla para almacenar las filas elegidas de los duplicados


% Iterar sobre cada paciente con duplicados
for i = 1:numel(uniqueRIDs)
    currentRID = uniqueRIDs(i);
    
    % Sub-tabla con TODAS las visitas del paciente en T_full, ordenadas por fecha (Contexto completo)
    T_patient_full = sortrows(T(T.RID == currentRID, :), 'EXAMDATE');
    
    % Sub-tabla con solo las visitas duplicadas del paciente actual
    T_patient_dups = T_duplicadas(T_duplicadas.RID == currentRID, :);
    
    % Agrupar por momento (M) dentro de las duplicadas del paciente
    uniqueMs = unique(T_patient_dups.M);
    
    % Iterar sobre cada momento (M) que está duplicado para este paciente
    for j = 1:numel(uniqueMs)
        currentM = uniqueMs(j);
        
        % Filas del grupo duplicado actual (RID, M)
        T_group = T_patient_dups(T_patient_dups.M == currentM, :);

        code_viscode = extractBetween(string(T_group.VISCODE),1,1)=="m" |...
                       extractBetween(string(T_group.VISCODE),1,1)=="b";

        if(isscalar(unique(code_viscode)))        
            % 4. Aplicar la Lógica de Coherencia Temporal
            %    Calcula la puntuación de coherencia para cada opción en T_group
            coherenceScores = calculateTemporalCoherence(T_group, T_patient_full, INTERVALS_DAYS);
            
            % Encontrar la visita con la MÁS BAJA desviación (mejor coherencia)
            [~, best_idx_in_group] = min(coherenceScores);
        else
            best_idx_in_group=find(code_viscode==1);
        end

            
            

        
        % La fila elegida es la que mejor se ajusta al patrón
        T_chosen_row = T_group(best_idx_in_group, :);
        
        % Agregar la fila elegida a la tabla de resultados
        T_dups_chosen = [T_dups_chosen; T_chosen_row];
    end
end

% 5. Combinar las filas no duplicadas y las filas duplicadas elegidas
T_sinDuplicados = [T_no_dups; T_dups_chosen];
T_sinDuplicados = sortrows(T_sinDuplicados, {'RID', 'EXAMDATE'});

%fprintf('Proceso completado. Filas finales: %d.\n', height(T_sinDuplicados));

end

% --- Función de Lógica de Coherencia Temporal (Clave) ---
function coherenceScores = calculateTemporalCoherence(T_group, T_patient_full, INTERVALS_DAYS)
% T_group: Sub-tabla de las visitas duplicadas (mismo RID, mismo M)
% T_patient_full: La tabla con TODAS las visitas del paciente (necesaria para el contexto)
% INTERVALS_DAYS: [Dias_6meses, Dias_12meses]

num_options = height(T_group);
coherenceScores = zeros(num_options, 1);
Dates_patient = T_patient_full.EXAMDATE;

% Iterar sobre cada visita duplicada para evaluar su coherencia
for k = 1:num_options
    current_date = T_group.EXAMDATE(k);
    
    % Encontrar la posición de esta fecha en la secuencia completa del paciente
    idx_in_full = find(Dates_patient == current_date, 1, 'first');
    
    % --- 1. Contexto Anterior (Prev_Date) ---
    if idx_in_full > 1
        prev_date = Dates_patient(idx_in_full - 1);
        diff_prev_days = days(current_date - prev_date);
        % Penalización: Distancia mínima a un múltiplo
        min_dev_prev = min(abs(diff_prev_days - INTERVALS_DAYS));
    else
        min_dev_prev = 0; % No penalizar si es la primera visita del paciente
    end
    
    % --- 2. Contexto Posterior (Next_Date) ---
    if idx_in_full < height(T_patient_full)
        next_date = Dates_patient(idx_in_full + 1);
        diff_next_days = days(next_date - current_date);
        % Penalización: Distancia mínima a un múltiplo
        min_dev_next = min(abs(diff_next_days - INTERVALS_DAYS));
    else
        min_dev_next = 0; % No penalizar si es la última visita del paciente
    end

    % Puntuación de coherencia: suma de las desviaciones mínimas
    % (La menor puntuación indica mejor coherencia con el patrón 6/12 meses)
    coherenceScores(k) = min_dev_prev + min_dev_next;
end
end