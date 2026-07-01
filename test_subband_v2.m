%% ========================================================================
%  项目实验1 — 子频带信号分析/合成系统 测试程序 v2
%  使用 MATLAB 的 filter() 函数，避免手动实现卷积的 bug
% =========================================================================

clear; clc; close all;

%% 加载系数
load('subband_filters.mat', 'hk_real', 'hk_imag', 'best_h', ...
     'best_order', 'D', 'Fs');
ORD2 = best_order;
FRM_LEN = 120 * D;
OMG_PASS = (15*pi)/(16*D);
OMG_STOP = pi/D;

fprintf('D=%d, N=%d, Fs=%dHz, FRM_LEN=%d\n', D, ORD2, Fs, FRM_LEN);

% 构建复数子带滤波器系数
hk = zeros(D, ORD2);
for k = 0:D-1
    hk(k+1,:) = hk_real(k+1,:) + 1j*hk_imag(k+1,:);
end

%% ====================================================================
%  测试①：子带幅频响应
% ====================================================================
fprintf('\n=== 测试①：子带滤波器幅频响应 ===\n');

w = linspace(0, pi, 10000);
figure('Name', '测试①');
subplot(1,2,1); hold on;
for k = 0:D-1
    Hk = freqz(hk(k+1,:), 1, w);
    plot(w/pi, 20*log10(abs(Hk)), 'LineWidth', 1);
end
xline(OMG_PASS/pi, 'k--'); xline(OMG_STOP/pi, 'k:');
xlabel('×π rad'); ylabel('dB'); title('子带幅频响应');
xlim([0,1]); ylim([-80,5]); grid on; legend(arrayfun(@(k)sprintf('k=%d',k),0:D-1,'un',0));

subplot(1,2,2); hold on;
for k = 0:D-1
    Hk = freqz(hk(k+1,:), 1, w);
    plot(w/pi, 20*log10(abs(Hk)), 'LineWidth', 1);
end
yline(-6, 'r--', '-6dB');
xlabel('×π rad'); ylabel('dB'); title('子带交叠（交点应在-6dB）');
xlim([0,1]); ylim([-10,3]); grid on;

% 检查 k=0 在 ωp 处的衰减 (手动计算频响，避免 freqz 标量歧义)
H0 = sum(best_h);
n_idx = (0:ORD2-1)';
Hk_wp = abs(sum(hk(1,:)' .* exp(-1j*OMG_PASS*n_idx)));
ripple_k0 = 20*log10(Hk_wp / abs(H0));
fprintf('  子带 k=0 @ ωp: %.2f dB ', ripple_k0);
if ripple_k0 > -0.4, fprintf('✅\n'); else, fprintf('⚠️\n'); end

%% ====================================================================
%  测试②：算法正确性 (用 filter 实现卷积)
% ====================================================================
fprintf('\n=== 测试②：算法正确性 ===\n');

% 生成测试信号
rng(0);
N_test = FRM_LEN * 20;  % 20帧
t = (0:N_test-1)' / Fs;
x = 0.5*sin(2*pi*440*t) + 0.3*sin(2*pi*880*t) + 0.2*sin(2*pi*1760*t) + 0.05*randn(N_test,1);

% 用 filter + 重叠保留法实现分帧处理
y = zeros(N_test, 1);
ana_zi = zeros(D, ORD2-1);   % 分析滤波器初始条件
syn_zi = zeros(D, ORD2-1);   % 合成滤波器初始条件

for frame = 0:floor(N_test/FRM_LEN)-1
    idx = (1:FRM_LEN) + frame*FRM_LEN;
    x_frame = x(idx);
    
    % --- 分析滤波：每路子带 ---
    ana_out = zeros(D, FRM_LEN);
    for k = 0:D-1
        [ana_out(k+1,:), ana_zi(k+1,:)] = filter(hk(k+1,:), 1, x_frame, ana_zi(k+1,:));
    end
    
    % --- 下采样：每D点保留一个，其余置零（上采样） ---
    syn_in = zeros(D, FRM_LEN);
    for k = 0:D-1
        syn_in(k+1, 1:D:end) = ana_out(k+1, 1:D:end);
    end
    
    % --- 合成滤波 ---
    syn_acc = zeros(1, FRM_LEN);
    for k = 0:D-1
        [syn_k, syn_zi(k+1,:)] = filter(hk(k+1,:), 1, syn_in(k+1,:), syn_zi(k+1,:));
        syn_acc = syn_acc + syn_k;
    end
    
    % 输出 (乘 D 归一化)
    y(idx) = real(syn_acc' * D);
end

% 计算 SNR (跳过前几帧的暂态)
skip = FRM_LEN * 3;
x_valid = x(skip+1:end);
y_valid = y(skip+1:end);
e = x_valid - y_valid;
SNR = 10*log10(var(x_valid) / var(e));
fprintf('  重建 SNR = %.2f dB\n', SNR);
if SNR > 20
    fprintf('  判定: ✅\n');
else
    fprintf('  判定: ⚠️ SNR偏低，检查滤波器实现\n');
end

figure('Name', '测试②');
subplot(3,1,1); plot(t, x); ylabel('幅度'); title('原始信号'); xlim([0, N_test/Fs]);
subplot(3,1,2); plot(t, y); ylabel('幅度'); title(sprintf('重建信号 (SNR=%.2fdB)', SNR)); xlim([0, N_test/Fs]);
subplot(3,1,3); plot(t, y-x); ylabel('幅度'); title('误差'); xlim([0, N_test/Fs]); xlabel('时间(s)');

%% ====================================================================
%  测试③：群延迟
% ====================================================================
fprintf('\n=== 测试③：群延迟 ===\n');
w_gd = linspace(0.01, pi-0.01, 2000);
for k = 0:D-1
    gd = grpdelay(hk(k+1,:), 1, w_gd);
    fprintf('  k=%d: 平均 %.2f 样点 (理论 %.1f)\n', k, mean(gd), (ORD2-1)/2);
end

%% ====================================================================
%  测试④：系统幅频响应 (扫频)
% ====================================================================
fprintf('\n=== 测试④：系统幅频响应 ===\n');
N_sweep = FRM_LEN * 30;
t_sweep = (0:N_sweep-1)'/Fs;
x_chirp = chirp(t_sweep, 50, t_sweep(end), Fs/2-100);

y_chirp = zeros(N_sweep, 1);
ana_zi = zeros(D, ORD2-1);
syn_zi = zeros(D, ORD2-1);
for frame = 0:floor(N_sweep/FRM_LEN)-1
    idx = (1:FRM_LEN) + frame*FRM_LEN;
    xf = x_chirp(idx);
    ao = zeros(D, FRM_LEN);
    for k = 0:D-1
        [ao(k+1,:), ana_zi(k+1,:)] = filter(hk(k+1,:), 1, xf, ana_zi(k+1,:));
    end
    si = zeros(D, FRM_LEN);
    for k = 0:D-1
        si(k+1, 1:D:end) = ao(k+1, 1:D:end);
    end
    sa = zeros(1, FRM_LEN);
    for k = 0:D-1
        [sk, syn_zi(k+1,:)] = filter(hk(k+1,:), 1, si(k+1,:), syn_zi(k+1,:));
        sa = sa + sk;
    end
    y_chirp(idx) = sa' * D;
end

nfft = 8192;
Xf = fft(x_chirp(1:nfft), nfft);
Yf = fft(y_chirp(1:nfft), nfft);
H_sys = Yf(1:nfft/2) ./ Xf(1:nfft/2);
f_ax = (0:nfft/2-1)'/nfft*Fs;

figure('Name', '测试④');
semilogx(f_ax, 20*log10(abs(H_sys)), 'b'); grid on;
xlabel('Hz'); ylabel('dB'); title('系统幅频响应');
xlim([20, Fs/2]); ylim([-6, 3]); hold on;
yline(0,'k'); yline(-3,'r--','-3dB');

%% ====================================================================
%  测试⑤：THD
% ====================================================================
fprintf('\n=== 测试⑤：THD ===\n');
test_freqs = [200, 500, 1000, 2000];
for fi = 1:length(test_freqs)
    ft = test_freqs(fi);
    Nt = round(Fs/ft)*50; Nt = floor(Nt/FRM_LEN)*FRM_LEN;
    tt = (0:Nt-1)'/Fs;
    xs = sin(2*pi*ft*tt);
    ys = zeros(Nt,1);
    ana_zi = zeros(D, ORD2-1); syn_zi = zeros(D, ORD2-1);
    for frame = 0:Nt/FRM_LEN-1
        idx = (1:FRM_LEN)+frame*FRM_LEN;
        xf = xs(idx);
        ao = zeros(D, FRM_LEN);
        for k = 0:D-1
            [ao(k+1,:), ana_zi(k+1,:)] = filter(hk(k+1,:), 1, xf, ana_zi(k+1,:));
        end
        si = zeros(D, FRM_LEN);
        for k = 0:D-1
            si(k+1,1:D:end) = ao(k+1,1:D:end);
        end
        sa = zeros(1, FRM_LEN);
        for k = 0:D-1
            [sk, syn_zi(k+1,:)] = filter(hk(k+1,:), 1, si(k+1,:), syn_zi(k+1,:));
            sa = sa + sk;
        end
        ys(idx) = real(sa' * D);
    end
    Yfft = fft(ys(round(end*0.3):end), 16384);
    Ym = abs(Yfft(1:end/2));
    [~, pi_] = max(Ym);
    fund = Ym(pi_)^2;
    harm = sum(Ym(2*pi_:10*pi_).^2);
    thd = sqrt(harm/fund)*100;
    fprintf('  %d Hz: THD = %.4f%%\n', ft, thd);
end

%% ====================================================================
%  存储与复杂度
% ====================================================================
fprintf('\n=== 测试⑥：存储 ===\n');
total = D*ORD2*2*8 + (ORD2-1+FRM_LEN)*8 + D*(ORD2-1+FRM_LEN)*2*8 + 3*FRM_LEN*2 + D*FRM_LEN*2*8;
fprintf('  总内存: %.1f KB (C6748 SRAM 256KB ', total/1024);
if total < 262144, fprintf('✅)\n'); else, fprintf('⚠️)\n'); end

fprintf('\n=== 测试⑦：复杂度 ===\n');
mac = 2*D*FRM_LEN*ORD2 * Fs/1024;
cpu = mac*2/456e6*100;
fprintf('  %.1f M MAC/s, CPU: %.1f%% ', mac/1e6, cpu);
if cpu < 80, fprintf('✅\n'); else, fprintf('⚠️\n'); end

fprintf('\n===== 测试完成 =====\n');
