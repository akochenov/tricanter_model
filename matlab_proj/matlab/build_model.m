function build_model(mdl, ddPath)
%BUILD_MODEL Собирает tricanter_top.slx программно: блоки, провода, решатель.
%
%   build_data_dictionary        % сначала словарь
%   build_model                  % потом модель
%
% Пересоздаёт модель с нуля. Ручные правки существующей .slx теряются.
%
% Структура:
%   Inputs (Constant)  ──┐
%   Params (Constant)  ──┼──> RHS ──> Cake mass ──┐
%                        │      ├──> phi_s state ─┤
%                        │      ├──> phi_so state┤ обратные связи
%                        │      ├──> phi_w state ─┤
%                        └──────┴──────────────────┘
%                               └──> log y (To Workspace)
%
% Владелец: A.

if nargin < 1 || isempty(mdl),    mdl    = 'tricanter_top';       end
if nargin < 2 || isempty(ddPath), ddPath = 'tricanter_data.sldd'; end

if ~exist(ddPath, 'file')
    error('Словаря %s нет. Сначала build_data_dictionary.', ddPath);
end

%% --- пересоздание ---
if bdIsLoaded(mdl), close_system(mdl, 0); end
if exist([mdl '.slx'], 'file'), delete([mdl '.slx']); end

new_system(mdl);
open_system(mdl);
set_param(mdl, 'DataDictionary', ddPath);

%% --- блоки ---
add_block('simulink/Sources/Constant', [mdl '/Inputs'], ...
    'Position', [40 60 140 100], ...
    'Value', 'TricInputs0', ...
    'OutDataTypeStr', 'Bus: BusTricInputs', ...
    'VectorParams1D', 'off');

add_block('simulink/Sources/Constant', [mdl '/Params'], ...
    'Position', [40 130 140 170], ...
    'Value', 'TricParams', ...
    'OutDataTypeStr', 'Bus: BusTricParams', ...
    'VectorParams1D', 'off');

add_block('simulink/User-Defined Functions/MATLAB Function', [mdl '/RHS'], ...
    'Position', [260 40 400 200]);

% Discrete-Time Integrator: масса кека. Нижняя отсечка обязательна —
% без неё при большой dn масса уходит в минус и Rtr превышает Rd.
add_block('simulink/Discrete/Discrete-Time Integrator', [mdl '/Cake mass'], ...
    'Position', [500 45 570 85], ...
    'IntegratorMethod', 'Integration: Forward Euler', ...
    'InitialCondition', 'm0', ...
    'SampleTime', 'Ts', ...
    'LimitOutput', 'on', ...
    'UpperSaturationLimit', 'inf', ...
    'LowerSaturationLimit', '0');

% Концентрации уже посчитаны точной экспонентой внутри функции.
% Накапливать их не надо — только вернуть на следующий шаг.
add_block('simulink/Discrete/Unit Delay', [mdl '/phi_s state'], ...
    'Position', [500 105 570 145], ...
    'InitialCondition', 'phi0', ...
    'SampleTime', 'Ts');

% Твёрдое, пришедшее внутри нефти — третья дисперсная популяция.
add_block('simulink/Discrete/Unit Delay', [mdl '/phi_so state'], ...
    'Position', [500 165 570 205], ...
    'InitialCondition', 'phi0', ...
    'SampleTime', 'Ts');

add_block('simulink/Discrete/Unit Delay', [mdl '/phi_w state'], ...
    'Position', [500 225 570 265], ...
    'InitialCondition', 'phi0', ...
    'SampleTime', 'Ts');

add_block('simulink/Sinks/To Workspace', [mdl '/log y'], ...
    'Position', [500 285 570 325], ...
    'VariableName', 'y', ...
    'SaveFormat', 'Structure With Time', ...
    'MaxDataPoints', 'inf');

%% --- тело функции и порты ---
set_rhs_body(mdl);

%% --- провода ---
% Порядок входов RHS: m, phi_s, phi_so, phi_w, u, pc
% Порядок выходов:    dm, phi_s_next, phi_so_next, phi_w_next, y
add_line(mdl, 'Inputs/1', 'RHS/5', 'autorouting', 'smart');
add_line(mdl, 'Params/1', 'RHS/6', 'autorouting', 'smart');

add_line(mdl, 'RHS/1', 'Cake mass/1',    'autorouting', 'smart');
add_line(mdl, 'RHS/2', 'phi_s state/1',  'autorouting', 'smart');
add_line(mdl, 'RHS/3', 'phi_so state/1', 'autorouting', 'smart');
add_line(mdl, 'RHS/4', 'phi_w state/1',  'autorouting', 'smart');
add_line(mdl, 'RHS/5', 'log y/1',        'autorouting', 'smart');

add_line(mdl, 'Cake mass/1',    'RHS/1', 'autorouting', 'smart');
add_line(mdl, 'phi_s state/1',  'RHS/2', 'autorouting', 'smart');
add_line(mdl, 'phi_so state/1', 'RHS/3', 'autorouting', 'smart');
add_line(mdl, 'phi_w state/1',  'RHS/4', 'autorouting', 'smart');

%% --- решатель ---
% Шаг задан буквами. Если разойдётся с Ts из словаря, экспонента
% посчитается для одного интервала, а шаг сделается на другой:
% ошибки не будет, результат тихо испортится.
set_param(mdl, ...
    'SolverType',   'Fixed-step', ...
    'Solver',       'FixedStepDiscrete', ...
    'FixedStep',    'Ts', ...
    'StartTime',    '0', ...
    'StopTime',     '1500');

save_system(mdl);
fprintf('Модель собрана: %s.slx\n', mdl);
fprintf('Проверка: Ctrl+D, затем Run. Ожидается 1500 шагов, всё по нулям.\n');

end


% =====================================================================
function set_rhs_body(mdl)
%SET_RHS_BODY Пишет код в MATLAB Function и настраивает порты.
%
% Ts объявлен как Parameter, а не Input: значение берётся из словаря,
% отдельный провод не нужен.

[n, J] = tric_dims();

root  = sfroot;
chart = root.find('-isa', 'Stateflow.EMChart', '-and', 'Path', [mdl '/RHS']);

chart.Script = sprintf([ ...
    'function [dm, phi_s_next, phi_so_next, phi_w_next, y] = ' ...
    'fcn(m, phi_s, phi_so, phi_w, u, pc, Ts)\n' ...
    '%%#codegen\n' ...
    '[dm, phi_s_next, phi_so_next, phi_w_next, y] = ' ...
    'tric_rhs(m, phi_s, phi_so, phi_w, u, pc, Ts);\n']);

% Порты создаются разбором кода; ждём и настраиваем.
setPort(chart, 'm',          'Input',  'double',              sprintf('[%d 1]', n));
setPort(chart, 'phi_s',      'Input',  'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'phi_so',     'Input',  'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'phi_w',      'Input',  'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'u',          'Input',  'Bus: BusTricInputs',  '-1');
setPort(chart, 'pc',         'Input',  'Bus: BusTricParams',  '-1');
setPort(chart, 'Ts',         'Parameter', 'double',           '1');

setPort(chart, 'dm',         'Output', 'double',              sprintf('[%d 1]', n));
setPort(chart, 'phi_s_next', 'Output', 'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'phi_so_next','Output', 'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'phi_w_next', 'Output', 'double',              sprintf('[%d %d]', n, J));
setPort(chart, 'y',          'Output', 'Bus: BusTricOut',     '-1');

end


% =====================================================================
function setPort(chart, name, scope, dtype, dims)
%SETPORT Настраивает одну переменную MATLAB Function блока.

d = chart.find('-isa', 'Stateflow.Data', '-and', 'Name', name);
if isempty(d)
    warning('Порт %s не найден — проверьте вручную в Edit Data.', name);
    return
end
d = d(1);
d.Scope = scope;
d.DataType = dtype;
if ~strcmp(dims, '-1')
    d.Props.Array.Size = dims;
end
end
