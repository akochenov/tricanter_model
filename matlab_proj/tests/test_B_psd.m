% tests/test_B_psd.m
[n, J] = tric_dims();

% твёрдое: x50, bs, rrsb из параметров
p = ref_load('psd_solids_3ph');
[x, w] = tric_psd(1.5e-6, 3.0, false);
fprintf('solids x:'); chk(x, p.x);
fprintf('solids w:'); chk(w, p.w);

% капли: d50, bw, rrsb всегда false
p = ref_load('psd_drops_3ph');
[x, w] = tric_psd(3.0e-6, 2.0, false);
fprintf('drops  x:'); chk(x, p.x);
fprintf('drops  w:'); chk(w, p.w);

% инварианты
assert(abs(sum(w) - 1) < 1e-15, 'доли не нормированы');
assert(all(diff(x) > 0),        'сетка не монотонна');
% --- общая сетка обеих твёрдых популяций ---
% При x50_so = 0 (наследование) сетка обязана совпасть с tric_psd
% побитово, а доли ws и wso — друг с другом.
[u, pc] = ref_case('3ph');
p = ref_load('psd_solids_3ph');
[xg, ws, wso] = tric_psd_solids(u, pc);
fprintf('общая x:'); chk(xg, p.x);
fprintf('общая w:'); chk(ws, p.w);
assert(isequal(ws, wso), 'при наследовании доли обеих популяций совпадают');
