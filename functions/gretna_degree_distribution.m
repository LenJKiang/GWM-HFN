function [Para, R2]= gretna_degree_distribution(Deg, Nbins)

%==========================================================================
% 此函数用于拟合网络G的度分布（也适用于介数和其他指标），包括三种类型：
%    类型1: 指数截断幂律拟合
%                p(k) = a * k^(alpha-1) * exp(-k/kc);
%    类型2: 幂律拟合
%                p(k) = a * k^(1-alpha);
%    类型3: 指数拟合
%                p(k) = a * exp(-alpha/k).
%
% 语法: [Para, R2]= gretna_degree_distribution_fit(Deg,Nbins)
%
% 输入:
%          Deg:
%             节点度（n*1数组）。
%          Nbins:
%             'Deg'中元素被分组的等间距容器数量（默认Nbins = 10）。
% 输出:
%          Para.exponentialtruncated:
%             指数截断幂律拟合的参数。
%          Para.powerlaw:
%             幂律拟合的参数。
%          Para.exponential:
%             指数拟合的参数。
%          R2.exponentialtruncated:
%             指数截断幂律拟合的拟合优度（值越接近1表示拟合越好）。
%          R2.powerlaw:
%             幂律拟合的拟合优度（值越接近1表示拟合越好）。
%          R2.exponential:
%             指数拟合的拟合优度（值越接近1表示拟合越好）。
%
% 参考文献:
% 1. Achard et al. (2006): A Resilient, Low-Frequency, Small-World Human
%    Brain Functional Network with Highly Connected Association Cortical Hubs.
% 2. Gong et al. (2009): Mapping Anatomical Connectivity Patterns of Human
%    Cerebral Cortex Using In Vivo Diffusion Tensor Imaging Tractography.
%
% Yong HE,     BIC,    MNI, McGill,  2006/08/01
% Jinhui WANG, NKLCNL, BNU, BeiJing, 2012/08/16, Jinhui.Wang.1982@gmail.com
%==========================================================================

% 如果只输入了一个参数，设置默认的Nbins值
if nargin == 1
    Nbins = 10; end

% 计算度的直方图
[n,xout] = hist(Deg, Nbins);

xdata = xout;

% 计算累积概率
ydata = fliplr(cumsum(n));
ydata = ydata./ydata(1);

% 设置非线性最小二乘法拟合选项
s = fitoptions('Method','NonlinearLeastSquares');

% 指数截断幂律拟合
f = fittype('a*x.^(b-1).*exp(-x./c)','independent','x','options',s);
[para1, gof1] = fit(xdata',ydata',f);

Para.exponentialtruncated.a     = para1.a;
Para.exponentialtruncated.alpha = para1.b;
Para.exponentialtruncated.kc    = para1.c;
R2.exponentialtruncated         = gof1.rsquare;
fity_truncated = para1.a*(xdata.^(-1+para1.b)).*exp(-xdata./para1.c);

% 幂律拟合
f = fittype('a*x.^(1-b)','independent','x','options',s);
[para2, gof2] = fit(xdata',ydata',f);

Para.powerlaw.a     = para2.a;
Para.powerlaw.alpha = para2.b;
R2.powerlaw         = gof2.rsquare;
fity_power = para2.a*(xdata.^(1-para2.b));

% 指数拟合
f = fittype('a.*exp(-b.*x)','independent','x','options',s);
[para3, gof3] = fit(xdata',ydata',f);

Para.exponential.a     = para3.a;
Para.exponential.alpha = para3.b;
R2.exponential         = gof3.rsquare;
fity_exponential       = para3.a.*exp(-para3.b.*xdata);

% 绘图
plot(log(xdata),log(ydata),'k*','MarkerSize',4)
xlabel('log(Degree)', 'FontSize', 23,'FontWeight','bold');
ylabel('log(Cumulative Distribution)', 'FontSize', 23,'FontWeight','bold');
hold on

plot(log(xdata),log(fity_truncated),'-','color','#D95319','Linewidth',2)
plot(log(xdata),log(fity_power),':k','Linewidth',2)
plot(log(xdata),log(fity_exponential),'--','color','#0072BD','Linewidth',2)

% 创建包含R²值的图例标签
truncated_label = sprintf('exponentially truncated power law (R² = %.4f)', R2.exponentialtruncated);
power_label = sprintf('power law (R² = %.4f)', R2.powerlaw);
exp_label = sprintf('exponential (R² = %.4f)', R2.exponential);

legend('real data',  truncated_label, power_label, exp_label, ...
    'Location', 'southwest', ...
    'FontSize', 16, ...
    'FontWeight','bold')

return