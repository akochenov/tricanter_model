function r = ref_load(name)
%REF_LOAD Читает эталон из ref/<name>.csv.
%
%   Для файлов вида "name,value" возвращает struct с полями по именам.
%   Для файлов с колонками возвращает struct с полями-векторами.
%
%   r = ref_load('geometry_3ph');   r.r_i, r.Qw, ...
%   r = ref_load('cake_steady_3ph'); r.Rtr, r.dH, r.tr, r.dm, ...

here = fileparts(mfilename('fullpath'));
path = fullfile(here, '..', 'ref_dyn', [name '.csv']);
t    = readtable(path, 'TextType', 'string');

if width(t) == 2 && all(strcmp(t.Properties.VariableNames, {'name','value'}))
    r = struct();
    for i = 1:height(t)
        r.(char(t.name(i))) = t.value(i);
    end
else
    r = table2struct(t, 'ToScalar', true);
end
end
