function [u, pc, dt] = ref_case(tag)
%REF_CASE Входы и параметры для эталонного режима. Владелец: B.
%
%   [u, pc, dt] = ref_case('3ph')   % трёхфазный, основной
%   [u, pc, dt] = ref_case('2ph')   % двухфазный, контрольные значения §11
%
% Источник истины — словарь tricanter_data.sldd. Если он не собран или
% Simulink недоступен, берутся литералы, совпадающие с build_data_dictionary.m
% (правя один файл, править и второй).
%
% Двухфазный режим получается наложением правок two_phase() из питона
% на дефолтный набор — ровно так же, как в dump_reference_dyn.py.
%
% frac_s_in_oil ЗДЕСЬ ВСЕГДА НОЛЬ. Эталоны в ref_dyn/ выгружены до
% появления третьей популяции, и модель совпадает с ними только при
% frac_s_in_oil = 0 — тогда твёрдого в нефти нет, Qo, Qw и все входные
% доли вырождаются в прежние формулы побитово. В словаре рабочее
% значение 0.25, как в питоне; для сверки с эталонами его надо занулять.

[u, pc] = base_case();
dt = 1.0;

% Сверка с эталонами возможна только на вырожденной модели.
u.frac_s_in_oil = 0.0;

switch tag
    case '3ph'
        % как есть
    case '2ph'
        % two_phase(): нефти нет, r_i = Rw, другое твёрдое
        u.eps_o  = 0.0;
        u.eps_wd = 0.0;
        u.eps_wf = 0.98;
        u.eps_s  = 0.02;
        u.C      = 250.0;
        u.dn     = 5.0;
        u.x50    = 2.5e-6;
        u.Rw     = 0.034;
        pc.rho_s   = 1410.0;
        pc.mu_s    = 0.3;
        pc.phi_sed = 0.5;
        pc.phi_ref = 0.5;
        pc.bs      = 2.78;
        pc.rrsb    = true;
    otherwise
        error('ref_case: неизвестный режим "%s", ожидается 3ph или 2ph', tag);
end
end


% =====================================================================
function [u, pc] = base_case()
%BASE_CASE Дефолтный набор: из словаря, иначе литералами.

ddPath = 'tricanter_data.sldd';
try
    dd  = Simulink.data.dictionary.open(ddPath);
    sec = dd.getSection('Design Data');
    pc  = unwrap(sec.getEntry('TricParams').getValue());
    u   = unwrap(sec.getEntry('TricInputs0').getValue());
    dd.close();
    return
catch
    warning('ref_case:noDict', ...
            'Словарь %s недоступен, беру литералы из ref_case.m', ddPath);
end

pc = struct('Rd', 0.040, 'Wsc', 0.025, 'Lpond', 0.17, 'alpha', 7.0, ...
            'Lsep', 0.0, 'rho_s', 2720.0, 'rho_w', 998.0, 'rho_o', 870.0, ...
            'eta_w', 1.0e-3, 'mu_s', 0.25, 'phi_sed', 0.45, 'phi_max', 0.64, ...
            'phi_ref', 0.13, 'u_conv', 0.57, 'bs', 3.0, 'bw', 2.0, ...
            'rrsb', false, 'b_so', 0.0, 'rrsb_so', -1.0, ...
            'inject_at_interface', false);

% Порядок полей обязан совпадать с BusTricInputs: новые поля в КОНЦЕ.
u  = struct('eta_o', 0.03, 'd50', 3.0e-6, 'x50', 1.5e-6, 'eps_s', 0.02, ...
            'eps_wd', 0.08, 'eps_wf', 0.30, 'eps_o', 0.60, 'Rw', 0.0287, ...
            'Ro', 0.028, 'Q', 30/3.6e6, 'dn', 10.0, 'C', 1500.0, ...
            'x50_so', 0.0, 'frac_s_in_oil', 0.25);
end

function v = unwrap(v)
if isa(v, 'Simulink.Parameter'), v = v.Value; end
end
