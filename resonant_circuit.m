function resonant_circuit()

%% 実験パラメータ
L = 6.00e-3;
C = 0.150e-6;
R_list = [20.0, 200, 400];
ratio = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95, 1.0, ...
         1.05, 1.1, 1.2, 1.5, 2.0, 4.0, 10.0];

%% 固有量
omega0 = 1/sqrt(L*C);
f0     = omega0/(2*pi);

fprintf('■ Rによらず一定な固有量\n');
fprintf('  L = %.4e H,  C = %.4e F\n', L, C);
fprintf('  ω0 = %.6g rad/s,  f0 = %.6g Hz\n\n', omega0, f0);

%% 出力先ディレクトリ
% filename('fullpath')
% → 実行中の この .m ファイルのフルパスを文字列として返す
% fileparts(...)
% → パスから「フォルダ部分」だけを取り出す
baseDir = fileparts(mfilename('fullpath'));
if isempty(baseDir), baseDir = pwd; end % もしbaseDirが空なら代わりにpwdを使う
outDir = fullfile(baseDir,'results'); % baseDirの中にresultsフォルダのパスを作る
if ~exist(outDir,'dir'), mkdir(outDir); end % outDirが存在しなければ作る

%% 計算＆CSV出力
% repmat(A,m,n)は配列Aをmxn行列に複製する関数
plotData = repmat(struct('R',0,'ratio',[],'gain_dB',[],'phase_deg',[]), numel(R_list),1);
all_gain_dB = []; % 縦軸範囲計算用

for k = 1:numel(R_list) %numel()は引数の要素数を返す
    R = R_list(k);
    zeta = (R/2)*sqrt(C/L);                          % 減衰係数
    omega = ratio * omega0;                          % 角周波数
    f = omega/(2*pi);                                % 周波数
    
    den_real = 1 - ratio.^2;                         % 周波数応答関数の分母の実部
    den_imag = 2*zeta.*ratio;                        % 周波数応答関数の分母の虚部

    gain_lin = 1 ./ sqrt(den_real.^2 + den_imag.^2); % ゲインの振幅比
    gain_dB  = 20*log10(gain_lin);                   % dBに変換
    phase_deg = -atan2(den_imag, den_real) * 180/pi; % 位相計算

    all_gain_dB = [all_gain_dB, gain_dB]; %#ok<AGROW>

    
    T = table(ratio(:), omega(:), f(:), gain_lin(:), gain_dB(:), phase_deg(:), ...
        'VariableNames', {'ω/ω_0','角周波数[rad/s]','周波数[Hz]','gain_linear','ゲイン[dB]','位相[deg]'});

    fprintf('■ R = %g Ω , ζ = %.6f のとき\n', R, zeta);
    disp(T); % コマンドウィンドウにテーブルを表示

    csvName = sprintf('共振回路_ゲイン・位相特性_%g.csv', R);
    writetable(T, fullfile(outDir, csvName));

    plotData(k).R = R;
    plotData(k).ratio = ratio;
    plotData(k).gain_dB = gain_dB;
    plotData(k).phase_deg = phase_deg;
end

%% ゲイン特性
fig1 = figure('Name','ゲイン(利得)－角周波数特性線図','Color','w');
hold on; grid on; set(gca,'XScale','log','XMinorGrid','on','YMinorGrid','on'); % x軸を対数スケールに変換
for k = 1:numel(plotData)
    semilogx(plotData(k).ratio, plotData(k).gain_dB, '-', ...
        'LineWidth', 1.6, 'DisplayName', sprintf('R=%gΩ',plotData(k).R));
end
xlabel('\omega / \omega_0'); ylabel('Gain [dB]');
title('ゲイン(利得)－角周波数特性線図');
y_min = min(all_gain_dB) - 3;
y_max = max(all_gain_dB) + 3;
ylim([y_min, y_max]);
legend('Location','best');
exportgraphics(fig1, fullfile(outDir,'ゲイン(利得)－角周波数特性線図.png'), 'Resolution', 200);

%% 位相特性
fig2 = figure('Name','位相－角周波数特性線図','Color','w');
hold on; grid on; set(gca,'XScale','log','XMinorGrid','on','YMinorGrid','on');
for k = 1:numel(plotData)
    semilogx(plotData(k).ratio, plotData(k).phase_deg, '-', ...
        'LineWidth', 1.6, 'DisplayName', sprintf('R=%gΩ',plotData(k).R));
end
xlabel('\omega / \omega_0'); ylabel('位相\phi [deg]');
title('位相－角周波数特性線図');
legend('Location','best');
exportgraphics(fig2, fullfile(outDir,'位相－角周波数特性線図.png'), 'Resolution', 200);

fprintf('\n■ 出力先: %s\n', outDir);
fprintf('  - ゲイン(利得)－角周波数特性線図.png\n');
fprintf('  - 位相－角周波数特性線図.png\n');
fprintf('  - 共振回路_ゲイン・位相特性_R.csv\n');
end
