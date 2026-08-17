function chk(got, want, tol)
%CHK Сравнение с эталоном. Порог по умолчанию 1e-12.

if nargin < 3, tol = 1e-12; end
if numel(got) ~= numel(want)
    fprintf('  РАЗМЕР  got %s, want %s\n', mat2str(size(got)), mat2str(size(want)));
    return
end
e  = max(abs(got(:) - want(:)));
st = 'ПРОВАЛ'; if e < tol, st = 'ок'; end
fprintf('  %-6s  max|d| = %.3e\n', st, e);
end
