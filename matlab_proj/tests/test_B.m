function test_B()
%TEST_B Все проверки участника B разом. Запуск: test_B

test_B_psd
test_B_grade
test_B_settle
test_B_cascade

fprintf('\n========================================\n');
fprintf('Модуль B пройден.\n');
fprintf('Оба уровня, step1 и steady, сверены с ref_dyn/ напрямую.\n');
fprintf('========================================\n');
end
