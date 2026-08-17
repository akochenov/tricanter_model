function project_status()
%PROJECT_STATUS  Готовность функций модели трикантера.
%
%   Файл считается заглушкой, если где-то в нём есть маркер ПОКА_ЗАГЛУШКА.
%   Убрали формулу из заглушки — уберите и маркер, иначе отчёт соврёт.
%
%   Список функций задан здесь литералами, а не выборкой по каталогу:
%   иначе отсутствующий файл просто не попадёт в отчёт, и это будет
%   выглядеть как будто его и не должно быть.

MARKER = 'ПОКА_ЗАГЛУШКА';

files = { ...
    'tric_dims',            'A'; ...
    'tric_geometry',        'A'; ...
    'tric_rhs',             'A'; ...
    'tric_psd',             'B'; ...
    'tric_grade',           'B'; ...
    'tric_settle_dyn',      'B'; ...
    'tric_cascade_solids',  'B'; ...
    'tric_cascade_drops',   'B'; ...
    'tric_cake_geom',       'C'; ...
    'tric_cake_balance',    'C'; ...
    'tric_outputs',         'D'; ...
    'tric_limits',          'D'};

n     = size(files, 1);
state = zeros(n, 1);            % 0 нет файла, 1 заглушка, 2 готово

fprintf('\n  функция                зона   состояние\n');
fprintf('  %s\n', repmat('-', 1, 48));

for i = 1:n
    name = files{i, 1};
    p    = which(name);

    if isempty(p)
        state(i) = 0;
        st = 'НЕТ ФАЙЛА';
    elseif contains(fileread(p), MARKER)
        state(i) = 1;
        st = 'заглушка';
    else
        state(i) = 2;
        st = 'готово';
    end

    fprintf('  %-22s %-6s %s\n', name, files{i, 2}, st);
end

fprintf('  %s\n', repmat('-', 1, 48));
fprintf('  готово %d из %d, заглушек %d, нет файла %d\n\n', ...
        sum(state == 2), n, sum(state == 1), sum(state == 0));
end