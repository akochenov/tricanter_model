function build_data_dictionary(ddPath)
%BUILD_DATA_DICTIONARY Собирает tricanter_data.sldd целиком: шины + значения.
%
%   build_data_dictionary                       % создаёт ./tricanter_data.sldd
%   build_data_dictionary('data/tricanter_data.sldd')
%
% ВНИМАНИЕ, СКРИПТ ОТСТАЛ ОТ РЕАЛЬНОГО СЛОВАРЯ. В tricanter_data.sldd
% сейчас есть то, чего здесь нет:
%   * шина BusDrives (M1, P1, M2, P2) и ветка drives в BusTricOut;
%   * флаги f_M1_limit, f_M2_limit в BusFlags (там 8 элементов, не 6);
%   * заполненные Description и DocUnits у всех элементов.
% Запуск пересоздаёт словарь С НУЛЯ и всё перечисленное сотрёт.
% Пока скрипт не догнал словарь — правьте .sldd руками по WIRING_OIL.md,
% а этот файл держите как справку о составе шин.
%
% Порядок элементов внутри шин критичен: он должен совпадать с порядком
% полей в структурах, которые собирает код. Не переставлять.
%
% Владелец: A.

if nargin < 1 || isempty(ddPath)
    ddPath = 'tricanter_data.sldd';
end

[n, J] = tric_dims();

%% --- пересоздание файла ---
if exist(ddPath, 'file')
    Simulink.data.dictionary.closeAll(ddPath, '-discard');
    delete(ddPath);
end
dd  = Simulink.data.dictionary.create(ddPath);
sec = dd.getSection('Design Data');

%% ---------- шины ----------

% Входы: 14 элементов
mkbus(sec, 'BusTricInputs', {
    'eta_o'   'double'  1
    'd50'     'double'  1
    'x50'     'double'  1
    'eps_s'   'double'  1
    'eps_wd'  'double'  1
    'eps_wf'  'double'  1
    'eps_o'   'double'  1
    'Rw'      'double'  1
    'Ro'      'double'  1
    'Q'       'double'  1
    'dn'      'double'  1
    'C'       'double'  1
    'x50_so'        'double'  1
    'frac_s_in_oil' 'double'  1 });

% Параметры: 20 элементов
mkbus(sec, 'BusTricParams', {
    'Rd'       'double'   1
    'Wsc'      'double'   1
    'Lpond'    'double'   1
    'alpha'    'double'   1
    'Lsep'     'double'   1
    'rho_s'    'double'   1
    'rho_w'    'double'   1
    'rho_o'    'double'   1
    'eta_w'    'double'   1
    'mu_s'     'double'   1
    'phi_sed'  'double'   1
    'phi_max'  'double'   1
    'phi_ref'  'double'   1
    'u_conv'   'double'   1
    'bs'       'double'   1
    'bw'       'double'   1
    'rrsb'     'boolean'  1
    'b_so'     'double'   1
    'rrsb_so'  'double'   1
    'inject_at_interface' 'boolean' 1 });

% Геометрия: 15 элементов
mkbus(sec, 'BusTricGeo', {
    'omega'    'double'   1
    'beta'     'double'   1
    'Lax'      'double'   1
    'L_cone'   'double'   1
    'f_clar'   'double'   1
    'u_ax'     'double'   1
    'r_i'      'double'   1
    'Qo'       'double'   1
    'Qw'       'double'   1
    'has_oil'  'boolean'  1
    'oil'      'boolean'  1
    'f_so'     'double'   1
    'phi_s0'   'double'   1
    'phi_so0'  'double'   1
    'phi_w0'   'double'   1 });

% Вложенные шины выхода
mkbus(sec, 'BusCake', {
    'Rtr'         'double'  n
    'dH'          'double'  n
    'm_cake'      'double'  n
    'm_cake_tot'  'double'  1
    'tau_s'       'double'  n
    'q_cake'      'double'  1 });

mkbus(sec, 'BusStreams', {
    'mdot_oil'      'double'  1
    'mdot_water'    'double'  1
    'mdot_cake'     'double'  1
    'x_w_in_oil'    'double'  1
    'x_o_in_water'  'double'  1
    'x_s_in_cake'   'double'  1 });

mkbus(sec, 'BusQuality', {
    'U'           'double'  1
    'E_s'         'double'  1
    'E_so'        'double'  1
    'E_w'         'double'  1
    'E_s_tot'     'double'  1
    'phi_s_out'   'double'  1
    'phi_so_out'  'double'  1
    'phi_w_out'   'double'  1
    'bsw'         'double'  1
    'phi_s_prof'  'double'  n
    'phi_so_prof' 'double'  n
    'phi_w_prof'  'double'  n
    'E_s_i'       'double'  n
    'tau_c'       'double'  1 });

mkbus(sec, 'BusFlags', {
    'bal_resid'         'double'   1
    'M_scroll'          'double'   1
    'f_oil_collapse'    'boolean'  1
    'f_water_collapse'  'boolean'  1
    'f_overload'        'boolean'  1
    'f_overfill'        'boolean'  1 });

% Корневая выходная шина
mkbus(sec, 'BusTricOut', {
    'cake'     'Bus: BusCake'     1
    'streams'  'Bus: BusStreams'  1
    'quality'  'Bus: BusQuality'  1
    'flags'    'Bus: BusFlags'    1 });

%% ---------- значения ----------

% Параметры. Порядок полей обязан совпасть с BusTricParams.
p = struct( ...
    'Rd',      0.040,  ...
    'Wsc',     0.025,  ...
    'Lpond',   0.17,   ...
    'alpha',   7.0,    ...
    'Lsep',    0.0,    ...
    'rho_s',   2720.0, ...
    'rho_w',   998.0,  ...
    'rho_o',   870.0,  ...
    'eta_w',   1.0e-3, ...
    'mu_s',    0.25,   ...
    'phi_sed', 0.45,   ...
    'phi_max', 0.64,   ...
    'phi_ref', 0.13,   ...
    'u_conv',  0.57,   ...
    'bs',      3.0,    ...
    'bw',      2.0,    ...
    'rrsb',    false,  ...
    'b_so',    0.0,    ...
    'rrsb_so', -1.0,   ...
    'inject_at_interface', false);

% Начальные входы. Порядок обязан совпасть с BusTricInputs.
% eta_o лежит здесь, а не в параметрах: сценарий 13.5 меняет её
% при охлаждении сырья.
inp = struct( ...
    'eta_o',  0.03,     ...
    'd50',    3.0e-6,   ...
    'x50',    1.5e-6,   ...
    'eps_s',  0.02,     ...
    'eps_wd', 0.08,     ...
    'eps_wf', 0.30,     ...
    'eps_o',  0.60,     ...
    'Rw',     0.0287,   ...
    'Ro',     0.028,    ...
    'Q',      30/3.6e6, ...
    'dn',     10.0,     ...
    'C',      1500.0,   ...
    'x50_so', 0.0,      ...
    'frac_s_in_oil', 0.25);

pPar = Simulink.Parameter(p);
pPar.DataType = 'Bus: BusTricParams';
pPar.CoderInfo.StorageClass = 'Auto';

iPar = Simulink.Parameter(inp);
iPar.DataType = 'Bus: BusTricInputs';
iPar.CoderInfo.StorageClass = 'Auto';

addEntry(sec, 'TricParams',  pPar);
addEntry(sec, 'TricInputs0', iPar);
addEntry(sec, 'm0',   zeros(n, 1));    % пустая машина
addEntry(sec, 'phi0', zeros(n, J));    % чистая жидкость в ячейках
addEntry(sec, 'Ts',   1.0);            % шаг; тот же в Solver

dd.saveChanges();
dd.close();

fprintf('Словарь собран: %s\n', ddPath);
fprintf('  шин: 8, значений: 5, ячеек n=%d, классов J=%d\n', n, J);
fprintf('  ВНИМАНИЕ: frac_s_in_oil = %.2f. Эталоны в ref_dyn/ выгружены\n', inp.frac_s_in_oil);
fprintf('  при frac_s_in_oil = 0 и действительны только на этом значении.\n');

end


% =====================================================================
function mkbus(sec, name, spec)
%MKBUS Создаёт Simulink.Bus и кладёт его в секцию словаря.
%   spec — cell {имя, тип, размер} по строкам, порядок значим.

b = Simulink.Bus;
b.Description = sprintf('Собрано build_data_dictionary, %s', datestr(now, 'yyyy-mm-dd'));

els = Simulink.BusElement.empty(0, 1);
for k = 1:size(spec, 1)
    e = Simulink.BusElement;
    e.Name       = spec{k, 1};
    e.DataType   = spec{k, 2};
    e.Dimensions = spec{k, 3};
    e.Complexity = 'real';
    els(k, 1) = e;
end
b.Elements = els;

addEntry(sec, name, b);
end


% =====================================================================
function addEntry(sec, name, value)
%ADDENTRY Добавляет или перезаписывает запись словаря.
try
    e = sec.getEntry(name);
    e.setValue(value);
catch
    sec.addEntry(name, value);
end
end
