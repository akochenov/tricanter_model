function test_B()
%TEST_B Все проверки участника B разом. Запуск: test_B

test_B_psd
test_B_grade
test_B_settle
test_B_cascade
test_B_oil

fprintf('\n========================================\n');
fprintf('Модуль B пройден.\n');
fprintf('Оба уровня, step1 и steady, сверены с ref_dyn/ напрямую\n');
fprintf('при frac_s_in_oil = 0. Третья популяция (твёрдое в нефти)\n');
fprintf('эталона не имеет и проверяется свойствами — см. test_B_oil.\n');
fprintf('========================================\n');
end
