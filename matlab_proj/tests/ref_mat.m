function M = ref_mat(name)
%REF_MAT Читает матричный эталон из ref_dyn/<name>.csv (без заголовка).
%
%   M = ref_mat('step1_Ts_3ph');   % 25x40

here = fileparts(mfilename('fullpath'));
M = readmatrix(fullfile(here, '..', 'ref_dyn', [name '.csv']));
end
