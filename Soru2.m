%Hande Yildirim, 141201047,Ele514,Odev2, Soru2

%Adım1
%a.
fs = 1024; %örnekleme frekansı
f = (0:1023);

f1 = fir1(60,[20 45]/(fs/2),'bandpass');
f2 = fir1(90,[20 45]/(fs/2),'bandpass');
f3 = fir1(120,[20 45]/(fs/2),'bandpass');
f4 = fir1(150,[20 45]/(fs/2),'bandpass');
f5 = fir1(180,[20 45]/(fs/2),'bandpass');

%Adım2
%b.
subplot(3,2,1)
u = [zeros(1,499) , 1 , zeros(1,500)]; %impulse oluşturuldu
t = 0:0.001:1-0.001;
plot(t,u,'r','linewidth',2)
grid minor
title('impulse input')
xlabel('t(s)')
ylabel('u(t)')

%c.
subplot(3,2,2)
U = fft(u,fs);
plot(f,abs(U),'r','linewidth',2)
grid minor
ylim([0 1.2])
title('impulse input in Freq. Domain')
ylabel('|U(f)|')
xlabel('f(Hz)')

%d.
subplot(3,2,3)
F1 = filter(f1,1,u);
F2 = filter(f2,1,u);
F3 = filter(f3,1,u);
F4 = filter(f4,1,u);
F5 = filter(f5,1,u);
plot(t,F1,'b','linewidth',2)
hold on
plot(t,F2,'g','linewidth',2)
plot(t,F3,'r','linewidth',2)
plot(t,F4,'c','linewidth',2)
plot(t,F5,'m','linewidth',2)
hold off
legend('60order','90order','120order','150order','180order')
grid minor
ylim([-0.06 0.06])
title('impulse output')
ylabel('y(t)')
xlabel('t(s)')

%e.
subplot(3,2,4)
F1_fft = fft(F1,fs);
F2_fft = fft(F2,fs);
F3_fft = fft(F3,fs);
F4_fft = fft(F4,fs);
F5_fft = fft(F5,fs);
plot(f,abs(F1_fft),'b','linewidth',2)
hold on
plot(f,abs(F2_fft),'g','linewidth',2)
plot(f,abs(F3_fft),'r','linewidth',2)
plot(f,abs(F4_fft),'c','linewidth',2)
plot(f,abs(F5_fft),'m','linewidth',2)
hold off
legend('60order','90order','120order','150order','180order')
grid minor
ylim([-0.2 1.2])
title('impulse output in Freq. Domain')
ylabel('|Y(f)|')
xlabel('f(Hz)')

%f.
subplot(3,2,5)
plot(f,abs(F1_fft)./abs(U),'b','linewidth',2)
hold on
plot(f,abs(F2_fft)./abs(U),'g','linewidth',2)
plot(f,abs(F3_fft)./abs(U),'r','linewidth',2)
plot(f,abs(F4_fft)./abs(U),'c','linewidth',2)
plot(f,abs(F5_fft)./abs(U),'m','linewidth',2)
hold off
legend('60order','90order','120order','150order','180order')
grid minor
xlim([0 100])
ylim([-0.2 1.2])
title('Amplitude Gain')
ylabel('|Y(f)|/|U(f)|')
xlabel('f(Hz)')

%g.
subplot(3,2,6)
db1 = 20*log10(abs(F1_fft)./abs(U));
db2 = 20*log10(abs(F2_fft)./abs(U));
db3 = 20*log10(abs(F3_fft)./abs(U));
db4 = 20*log10(abs(F4_fft)./abs(U));
db5 = 20*log10(abs(F5_fft)./abs(U));
semilogx(f,db1,'b','linewidth',2)
hold on
semilogx(f,db2,'g','linewidth',2)
semilogx(f,db3,'r','linewidth',2)
semilogx(f,db4,'c','linewidth',2)
semilogx(f,db5,'m','linewidth',2)
hold off
legend('60order','90order','120order','150order','180order')
grid minor
ylim([-100 0])
title('Amplitude Gain dB')
ylabel('|Y(f)|/|U(f)| dB')
xlabel('f(Hz)')

%Adım3
%a.%filtreler oluşturuldu
[y1, f1] = butter(2,[20 45]/(fs/2),'bandpass');
[y2, f2] = butter(3,[20 45]/(fs/2),'bandpass');
[y3, f3] = butter(4,[20 45]/(fs/2),'bandpass');

%b.dürtü sinyali üretildi
figure();
subplot(3,2,1)
u = [zeros(1,499) , 1 , zeros(1,500)]; %impulse oluşturuldu
t = 0:0.001:1-0.001;
plot(t,u,'r','linewidth',2)
grid minor
title('impulse input')
xlabel('t(s)')
ylabel('u(t)')

%c.
subplot(3,2,2)
U = fft(u,fs); %dürtü sinyalinin frekans bileşenleri fft komutu ile elde edildi
plot(f,abs(U),'r','linewidth',2)
grid minor
ylim([0 1.4])
title('impulse input in Freq. Domain')
ylabel('|U(f)|')
xlabel('f(Hz)')

%d.
subplot(3,2,3)
F1 = filter(y1,f1,u);
F2 = filter(y2,f2,u);
F3 = filter(y3,f3,u);
plot(t,F1,'b','linewidth',2)
hold on
plot(t,F2,'g','linewidth',2)
plot(t,F3,'r','linewidth',2)
hold off
legend('2order','3order','4order')
grid minor
ylim([-0.065 0.055])
title('impulse output')
ylabel('y(t)')
xlabel('t(s)')

%e.
subplot(3,2,4)
F1_fft = fft(F1,fs);
F2_fft = fft(F2,fs);
F3_fft = fft(F3,fs);
plot(f,abs(F1_fft),'b','linewidth',2)
hold on
plot(f,abs(F2_fft),'g','linewidth',2)
plot(f,abs(F3_fft),'r','linewidth',2)
hold off
legend('2order','3order','4order')
grid minor
ylim([0 1.2])
title('impulse output in Freq. Domain')
ylabel('|Y(f)|')
xlabel('f(Hz)')

%f.
subplot(3,2,5)
plot(f,abs(F1_fft)./abs(U),'b','linewidth',2)
hold on
plot(f,abs(F2_fft)./abs(U),'g','linewidth',2)
plot(f,abs(F3_fft)./abs(U),'r','linewidth',2)
hold off
legend('2order','3order','4order')
grid minor
xlim([0 1200])
ylim([0 1.2])
title('Amplitude Gain')
ylabel('|Y(f)|/|U(f)|')
xlabel('f(Hz)')

%g.
subplot(3,2,6)
db1 = 20*log10(abs(F1_fft)./abs(U));
db2 = 20*log10(abs(F2_fft)./abs(U));
db3 = 20*log10(abs(F3_fft)./abs(U));
semilogx(f,db1,'b','linewidth',2)
hold on
semilogx(f,db2,'g','linewidth',2)
semilogx(f,db3,'r','linewidth',2)
hold off
legend('2order','3order','4order')
grid minor
ylim([-200 0])
xlim([1 500])
title('Amplitude Gain dB')
ylabel('|Y(f)|/|U(f)| dB')
xlabel('f(Hz)')