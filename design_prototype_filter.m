%% ========================================================================
%  项目实验1 — 子频带信号分析与合成系统
%  原型低通滤波器系数搜索程序 (MATLAB 版)
%  参考：项目实验指导书 第3.2.5节
%  功能：
%    1. 搜索满足指标的最优原型低通滤波器系数 h(n)
%    2. 生成 D=8 路子带滤波器系数 hk_real / hk_imag
%    3. 绘制滤波器幅频响应图
%    4. 导出可直接用于 CCS C 工程的系数文件
% =========================================================================

clear; clc; close all;

%% ======================== 参数定义 ======================================
D = 8;                              % 子频带数
OMG_PASS = (15*pi)/(16*D);          % 通带截止角频率 (rad)
OMG_STOP  = pi/D;                   % 阻带起始角频率 (rad)
Fs = 8000;                          % 采样频率 (Hz)，参考示例

% 通带允许最大波动 0.4 dB，阻带要求最小衰减 60 dB
PASS_RIPPLE_DB = -0.4;              % 通带纹波上限 (dB)
STOP_ATTEN_DB  = -60;               % 阻带衰减下限 (dB)

fprintf('====================================================\n');
fprintf('  子频带信号分析/合成系统 — 原型滤波器设计\n');
fprintf('====================================================\n');
fprintf('  子频带数 D             = %d\n', D);
fprintf('  通带角频率 ωp          = %.10f rad\n', OMG_PASS);
fprintf('  阻带角频率 ωs          = %.10f rad\n', OMG_STOP);
fprintf('  通带最大波动           = %.1f dB\n', PASS_RIPPLE_DB);
fprintf('  阻带最小衰减           = %.1f dB\n', STOP_ATTEN_DB);
fprintf('  采样频率 Fs            = %d Hz\n', Fs);
fprintf('  滤波器阶数 N           = q × D (q=10~399)\n');
fprintf('====================================================\n\n');

%% ======================== 搜索最优滤波器参数 ============================
% 对 q=10~399（即 N=80~3192），以及对每个 q 调整理想低通截止频率 omiga，
% 找到一组满足通带/阻带指标的原型滤波器系数。

flag = 0;       % 找到标志
best_q = 0;
best_order = 0;
best_omega = 0;
best_h = [];

fprintf('正在搜索最优滤波器系数...\n');
fprintf('  q 范围      : 10 ~ 399\n');
fprintf('  omega 搜索  : %.10f → %.10f, 步长 1e-5\n\n', OMG_STOP, OMG_PASS);

tic;
for q = 10:399
    N = q * D;                      % 滤波器阶数
    
    % 对 omiga 从 OMG_STOP 向下搜索到 OMG_PASS
    for omega = OMG_STOP:-1e-5:OMG_PASS
        
        % --- 用 Hanning 窗设计原型低通滤波器 ---
        h = zeros(1, N);
        sum_h    = 0;    % H(0)        — 直流增益
        re_pass  = 0;    % H(ωp) 实部
        im_pass  = 0;    % H(ωp) 虚部
        re_stop  = 0;    % H(ωs) 实部
        im_stop  = 0;    % H(ωs) 虚部
        
        for n = 0:N-1
            % 理想低通冲激响应 (sinc 函数)
            t = n - (N-1)/2;
            if abs(t) < 1e-12
                vv = omega / pi;    % sin(0)/0 → 1, 乘以 ω/π
            else
                vv = sin(t * omega) / (pi * t);
            end
            
            % Hanning 窗加权 (周期型)
            h(n+1) = vv * (0.5 - 0.5*cos((2*n+1)*pi/N));
            
            % 累加频率响应
            sum_h = sum_h + h(n+1);                         % H(0)
            
            phi = n * OMG_PASS;
            re_pass = re_pass + h(n+1) * cos(phi);          % H(ωp)
            im_pass = im_pass + h(n+1) * sin(phi);
            
            phi = n * OMG_STOP;
            re_stop = re_stop + h(n+1) * cos(phi);          % H(ωs)
            im_stop = im_stop + h(n+1) * sin(phi);
        end
        
        % 归一化幅频响应 (dB)
        H0_sq     = sum_h^2;
        H_pass_sq = re_pass^2 + im_pass^2;
        H_stop_sq = re_stop^2 + im_stop^2;
        
        ripple_pass = 10*log10(H_pass_sq / H0_sq);   % 通带纹波 (dB)
        atten_stop  = 10*log10(H_stop_sq / H0_sq);   % 阻带衰减 (dB)
        
        % 检查是否满足指标
        if (ripple_pass > PASS_RIPPLE_DB) && (atten_stop < STOP_ATTEN_DB)
            flag = 1;
            best_q = q;
            best_order = N;
            best_omega = omega;
            best_h = h;
            break;  % 跳出 omega 循环
        end
    end
    
    if flag == 1
        break;  % 跳出 q 循环
    end
    
    % 进度显示
    if mod(q, 20) == 0
        fprintf('  已搜索到 q=%d, 当前最优 N=%d ...\n', q, q*D);
    end
end
elapsed = toc;

%% ======================== 输出搜索结果 ==================================
if flag == 1
    fprintf('\n✅ 找到满足要求的滤波器!\n');
    fprintf('  q               = %d\n', best_q);
    fprintf('  滤波器阶数 N    = %d\n', best_order);
    fprintf('  实际截止频率    = %.10f rad\n', best_omega);
    fprintf('  搜索耗时        = %.2f 秒\n\n', elapsed);
    
    % === 用 freqz 精确验证 ===
    H0 = sum(best_h);
    w_axis = linspace(0, pi, 8000);
    H_w = freqz(best_h, 1, w_axis);
    H_dB = 20*log10(abs(H_w) / abs(H0));
    
    H_all = freqz(best_h, 1, [OMG_PASS, OMG_STOP]);
    H_pass = abs(H_all(1));
    H_stop = abs(H_all(2));
    actual_pass_ripple = 20*log10(H_pass / abs(H0));
    actual_stop_atten  = 20*log10(H_stop / abs(H0));
    
    fprintf('  验证结果:\n');
    fprintf('    通带纹波 (ωp=%.6f): %.4f dB  (要求 > %.1f dB)\n', ...
            OMG_PASS, actual_pass_ripple, PASS_RIPPLE_DB);
    fprintf('    阻带衰减 (ωs=%.6f): %.4f dB  (要求 < %.1f dB)\n', ...
            OMG_STOP, actual_stop_atten, STOP_ATTEN_DB);
    
    %% ==================== 导出滤波器系数 ================================
    fprintf('\n📄 导出滤波器系数到 h_coeff.txt ...\n');
    fid = fopen('h_coeff.txt', 'w');
    fprintf(fid, '/* ========================================================\n');
    fprintf(fid, ' * 子频带信号分析/合成系统 — 原型低通滤波器系数\n');
    fprintf(fid, ' * 生成方式: Hanning 窗函数法\n');
    fprintf(fid, ' * 参数: q=%d, N=%d, D=%d, omega_s=%.10f\n', ...
            best_q, best_order, D, best_omega);
    fprintf(fid, ' * 通带纹波: %.2f dB @ ωp=%.6f\n', actual_pass_ripple, OMG_PASS);
    fprintf(fid, ' * 阻带衰减: %.2f dB @ ωs=%.6f\n', actual_stop_atten, OMG_STOP);
    fprintf(fid, ' * ======================================================== */\n\n');
    fprintf(fid, '#ifndef H_COEFF_H\n');
    fprintf(fid, '#define H_COEFF_H\n\n');
    fprintf(fid, '#define ORD2 %d\n\n', best_order);
    fprintf(fid, 'double h[ORD2] = {\n');
    for n = 0:length(best_h)-1
        if mod(n, 4) == 0
            fprintf(fid, '    ');
        end
        fprintf(fid, '%.15e, ', best_h(n+1));
        if mod(n, 4) == 3
            fprintf(fid, '\n');
        end
    end
    if mod(length(best_h), 4) ~= 0
        fprintf(fid, '\n');
    end
    fprintf(fid, '};\n\n');
    fprintf(fid, '#endif /* H_COEFF_H */\n');
    fclose(fid);
    fprintf('  ✅ 已保存到 h_coeff.txt (可直接用于 CCS C 工程)\n');
    
    %% ==================== 生成子带滤波器系数 ============================
    fprintf('\n📄 生成子带滤波器系数 (D=%d)...\n', D);
    hk_real = zeros(D, best_order);
    hk_imag = zeros(D, best_order);
    for k = 0:D-1
        for n = 0:best_order-1
            theta = 2*pi*n*k/D;
            hk_real(k+1, n+1) = best_h(n+1) * cos(theta);
            hk_imag(k+1, n+1) = best_h(n+1) * sin(theta);
        end
    end
    
    % 保存子带系数到 .mat 文件
    save('subband_filters.mat', 'hk_real', 'hk_imag', 'best_h', ...
         'best_order', 'D', 'Fs', 'OMG_PASS', 'OMG_STOP');
    fprintf('  ✅ 子带滤波器系数已保存到 subband_filters.mat\n');
    
    %% ==================== 绘图 ==========================================
    figure('Name', '项目实验1 — 滤波器设计', 'Position', [50, 50, 1400, 900]);
    
    % (1) 原型滤波器时域系数
    subplot(2, 4, 1);
    stem(0:best_order-1, best_h, 'b.', 'MarkerSize', 2);
    title(sprintf('原型滤波器 h(n) (N=%d)', best_order));
    xlabel('n'); ylabel('h(n)'); grid on; xlim([-5, best_order+5]);
    
    % (2) 幅频响应 (dB) — 全频段
    subplot(2, 4, 2);
    plot(w_axis/pi, H_dB, 'b-', 'LineWidth', 1.2); hold on;
    xline(OMG_PASS/pi, 'r--', 'ωp', 'LineWidth', 1);
    xline(OMG_STOP/pi, 'g--', 'ωs', 'LineWidth', 1);
    yline(PASS_RIPPLE_DB, 'r:', sprintf('%.1f dB', PASS_RIPPLE_DB), 'LineWidth', 1);
    yline(STOP_ATTEN_DB,  'g:', sprintf('%.1f dB', STOP_ATTEN_DB),  'LineWidth', 1);
    title('原型滤波器幅频响应'); grid on;
    xlabel('归一化频率 (×π rad)'); ylabel('幅度 (dB)');
    xlim([0, 1]); ylim([-80, 5]);
    legend('|H(e^{jω})|', 'ωp', 'ωs', 'Location', 'best');
    
    % (3) 通带细节
    subplot(2, 4, 3);
    idx_pass = w_axis <= OMG_PASS*1.2;
    plot(w_axis(idx_pass)/pi, H_dB(idx_pass), 'b-', 'LineWidth', 1.2); hold on;
    yline(PASS_RIPPLE_DB, 'r--', sprintf('%.1f dB', PASS_RIPPLE_DB));
    title('通带细节'); grid on;
    xlabel('归一化频率 (×π rad)'); ylabel('幅度 (dB)');
    ylim([-1, 0.5]);
    
    % (4) 阻带细节
    subplot(2, 4, 4);
    idx_stop = w_axis >= OMG_STOP*0.8;
    plot(w_axis(idx_stop)/pi, H_dB(idx_stop), 'b-', 'LineWidth', 1.2); hold on;
    yline(STOP_ATTEN_DB, 'r--', sprintf('%.1f dB', STOP_ATTEN_DB));
    title('阻带细节'); grid on;
    xlabel('归一化频率 (×π rad)'); ylabel('幅度 (dB)');
    xlim([OMG_STOP*0.8/pi, 1]); ylim([-80, 0]);
    
    % (5) 子带滤波器组 — 幅度响应
    subplot(2, 4, [5, 6]);
    colors = lines(D);
    w_fft = linspace(0, pi, 2048);
    hold on;
    for k = 0:D-1
        Hk = freqz(hk_real(k+1,:) + 1j*hk_imag(k+1,:), 1, w_fft);
        plot(w_fft/pi, 20*log10(abs(Hk)/abs(H0)), ...
            'Color', colors(k+1,:), 'LineWidth', 1.2);
    end
    title(sprintf('子带滤波器组 |H_k(e^{jω})|  (D=%d)', D));
    xlabel('归一化频率 (×π rad)'); ylabel('幅度 (dB)');
    xlim([0, 1]); ylim([-80, 5]); grid on;
    
    % (6) 子带滤波器组 — 相位响应 (前4路)
    subplot(2, 4, [7, 8]);
    hold on;
    for k = 0:min(3, D-1)
        Hk = freqz(hk_real(k+1,:) + 1j*hk_imag(k+1,:), 1, w_fft);
        plot(w_fft/pi, unwrap(angle(Hk)), ...
            'Color', colors(k+1,:), 'LineWidth', 1.2);
    end
    title('子带滤波器相位响应 (前4路)');
    xlabel('归一化频率 (×π rad)'); ylabel('相位 (rad)');
    xlim([0, 1]); grid on;
    legend(arrayfun(@(k) sprintf('k=%d', k), 0:min(3,D-1), 'UniformOutput', false));
    
    sgtitle(sprintf('项目实验1: 子频带信号分析/合成系统  q=%d, N=%d', best_q, best_order));
    
    %% ==================== 打印综合信息 ==================================
    fprintf('\n');
    fprintf('====================================================\n');
    fprintf('  设计完成! 关键参数总结\n');
    fprintf('====================================================\n');
    fprintf('  子频带数 D          = %d\n', D);
    fprintf('  原型滤波器阶数 N    = %d\n', best_order);
    fprintf('  q                   = %d\n', best_q);
    fprintf('  每组数据帧长度      = 120×D = %d\n', 120*D);
    fprintf('  通带纹波 @ ωp       = %.4f dB\n', actual_pass_ripple);
    fprintf('  阻带衰减 @ ωs       = %.4f dB\n', actual_stop_atten);
    fprintf('  输出文件:\n');
    fprintf('    - h_coeff.txt        (CCS C 头文件)\n');
    fprintf('    - subband_filters.mat (MATLAB 数据)\n');
    fprintf('====================================================\n');
    
else
    fprintf('\n❌ 未找到满足要求的滤波器!\n');
    fprintf('   建议: 增大 q 的搜索范围或放宽阻带衰减指标。\n');
end
