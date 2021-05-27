% clear all
clear I_HDmax
clear I_HD
% load('C:\Users\temp\Downloads\ber_vis_d32.mat')
load('C:\Users\Owner\Desktop\qkd\clicks_up_to_63.mat')
% ber_d =    [1.7104    2.3332    3.5656    5.8914    9.8209 3.5656/7];
ber_d =    [1.9080    2.5803    3.9062    6.7852   11.8197 11.819/31 3*11.819/31 7*11.819/31 15*11.819/31];
vmin = 1.0e+03*[ 1.6207    0.3982    0.1072    0.0305    0.0074 ];
vmax = 1.0e+04*[8.0037    3.8732    1.5522    0.5210    0.1582];
a = vmax(5);
vmax_t = [a*256 a*64 a*16 a*4 a];
vis_d =   100*(vmax_t-vmin)./(vmax_t+vmin);
vis_d_max =   100*((vmax_t-(vmin+sqrt(vmin)))./(vmax_t+(vmin+sqrt(vmin))));

vis_d = [vis_d vis_d(5) vis_d(5) vis_d(5) vis_d(5)];
% load('C:\Users\temp\Downloads\clicks_up_to_63.mat')
d_vector = [2,4,8,16,32,2,4,8,16];
clicks = clicks_tot1/2.1475;
%ֵֵֵֵֵ0.196
%cal = 10;
cal = 35;
%0.00056
mu_vector = cal*0.0056*10^(-40.5/10)*ones(1,64);
for i =1:64
    mu_vector(i) = mu_vector(i)*(10^(1/10))^(i-1);
end
%2.326
clicks1 = 800+(2e-9*cal*kron([2,4,8,16]',(1./mu_vector))+(10/2.32)*1e-6).^-1;
clicks = [clicks(1:5,:);clicks1];
% clicks = [clicks(1:5,:); clicks(1:4,:)];
% Q_vector = 0.01*[2,3.2,5.1,8.5,13,22,38]./(d_vector-1);
% V_vector = 0.01*[96,96,96,96,96,96,96];
Q_vector = 0.01*ber_d./(d_vector-1);
V_vector = 0.01*vis_d;
%0.0196 @ 31dB


% mu_vector = 0.005*[1,1.2589,1.2589^2,1.2589^3,1.2589^4,1.2589^5,1.2589^6];

for i = 1:size(clicks,1)
    for j = 1:size(mu_vector,2)
d = d_vector(i);
Q = Q_vector(i);
V = V_vector(i);
mu = mu_vector(j);
p = abs(exp(-mu/2).*sqrt(V)-sqrt(1-exp(-mu)).*sqrt(1-V));
a = (Q./d).*((d-2)*exp(-mu)+1);
b = (Q./d).*(1-exp(-mu));
a1 = ((1-(d-1).*Q)./d) .* ((d-1).*p.^2+1);
b1 = ((1-(d-1).*Q)./d) .* (1-p.^2);
holevo_hd = d.*(S(a)+(d-2).*S(b)-((d-1)./d).*S(Q))+S(a1)+(d-1).*S(b1)-S(1-(d-1).*Q);
I_HD(i,j) = log2(d)+(d-1).*Q.*log2(Q)+(1-Q.*(d-1)).*log2(1-Q.*(d-1))-holevo_hd;
    end
end
%%
figure(8)
for i = 1:4
    subplot(1,3,1)
    semilogx(mu_vector,I_HD(i,:),'LineWidth',2)
    hold on
    grid on
    xlabel('Occupation')
    ylabel('Secure bits per photon')
    legend('D=2','D=4','D=8','D=16')

    subplot(1,3,2)
    loglog(mu_vector,clicks(i,:),'LineWidth',2)
    hold on
    grid on
    xlabel('Occupation')
    ylabel('Detected photons per second')
    legend('D=2','D=4','D=8','D=16')

    subplot(1,3,3)
    semilogx(mu_vector,I_HD(i,:).*clicks(i,:),'LineWidth',2)
    hold on
    grid on
    xlabel('Occupation')
    ylabel('Secure bits per second')
    legend('D=2','D=4','D=8','D=16')

end
% colros = [['#0072BD'],['#D95319'],['#EDB120'],['#7E2F8E']];
colors = [[0 0.4470 0.7410];[0.8500 0.3250 0.0980];[0.9290 0.6940 0.1250];[0.4940 0.1840 0.5560]];
for i = 6:9
    subplot(1,3,1)
    loglog(mu_vector,I_HD(i,:),'--','color', colors(i-5,:),'LineWidth',2)
    set(gca,'FontSize',14)
    hold on
    grid on
    xlabel('Occupation')
    ylabel('Secure bits per photon')
    legend('D=2','D=4','D=8','D=16','D=2 model','D=4 model','D=8 model','D=16 model')
    xlim([10^-3 10^0])
    ylim([0 4])

    subplot(1,3,2)
    loglog(mu_vector,clicks(i,:),'--','color', colors(i-5,:),'LineWidth',2)
    set(gca,'FontSize',14)
    hold on
    grid on
    xlabel('Occupation')
    legend('D=2','D=4','D=8','D=16','D=2 model','D=4 model','D=8 model','D=16 model')   
    xlim([10^-3 10^0])
    ylim([1000 3*10^5])
    
    subplot(1,3,3)
    semilogx(mu_vector,I_HD(i,:).*clicks(i,:),'--','color', colors(i-5,:),'LineWidth',2)
    set(gca,'FontSize',14)
    hold on
    grid on
    xlabel('Occupation')
    ylabel('Secure bits per second')
    legend('D=2','D=4','D=8','D=16','D=2 model','D=4 model','D=8 model','D=16 model')
    xlim([10^-3 10^0])
    ylim([1000 2*10^5])
end
set(gcf, 'Position',  [0, 0, 500, 800])

%%
figure(17)

yyaxis left
% semilogx(d_vector(1:5),Q_vector(1:5)*100,'.','MarkerSize',20);
% set(gca,'FontSize',12)
errorbar(d_vector(1:5),Q_vector(1:5)*100,1./sqrt(Q_vector(1:5).*[clicks(1,41) clicks(2,41) clicks(3,41) clicks(4,41) clicks(5,41)]),'-','LineWidth',2)
set(gca,'XScale','log');
set(gca,'FontSize',12)
ylim([0 2])
xlabel('Dimension')
ylabel('QBER %')
yyaxis right
% semilogx(d_vector(1:5),V_vector(1:5)*100,'.','MarkerSize',20);
% set(gca,'FontSize',12)
hold on
errorbar(d_vector(1:5),V_vector(1:5)*100,V_vector(1:5)*100 - vis_d_max,'-','LineWidth',2)
% set(gca,'YScale','log');
set(gca,'FontSize',12)
ylim([90 100])
ylabel('Visibility %')
xlim([2,32])
grid on

%%
figure(10)
k = 1000;
Q1 = linspace(0.0001,0.22,k);
clear I_HDmax
clear I_HD
for i = 1:4
    for j = 1
d = 2^i;
Q = Q1./(1+d*Q1);
mu = j/10;
V = 0.98;
p = abs(exp(-mu/2).*sqrt(V)-sqrt(1-exp(-mu)).*sqrt(1-V));
a = (Q/d)*((d-2)*exp(-mu)+1);
b = (Q/d)*(1-exp(-mu));
a1 = ((1-(d-1)*Q)/d) * ((d-1)*p.^2+1);
b1 = ((1-(d-1)*Q)/d) * (1-p.^2);
holevo_hd = d*(S(a)+(d-2)*S(b)-((d-1)/d)*S(Q))+S(a1)+(d-1)*S(b1)-S(1-(d-1)*Q);
I_HD(log2(d),:,j) = (log2(d)+(d-1)*Q.*log2(Q)+(1-Q*(d-1)).*log2(1-Q*(d-1))-holevo_hd);
    end
end
I_HDmax(:,:) = max(real(I_HD),[],3);
plot(100*Q1,I_HDmax','linewidth',2)
set(gca,'FontSize',12)
ylim([0 3])
grid minor
legend('D=2','D=4','D=8','D=16')
ylabel('Secure bits per photon')
xlabel('QBER (%)')

%%
function ans = S(x)
    ans = -x.*log2(x);
end
