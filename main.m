%Hande Yildirim, 141201047,Ele514,Odev2, Soru1

%SORU1

%Adim 1
soru("sound01.wav"); %3 adet dosya icin kodu calistirdik 

% Asagida yorum halinde bulunan fonksiyonlar scriptin el altındadir. Kodun
% calisabilmesi icin bu durum yapildi.

% function soru(dos)
% [oku, fs] = audioread(dos);%Ornekleme frekansi
% t = linspace(0, length(oku)/fs, length(oku));
% %Sinyal zamani sampleSayisi/fs
% %t ekseninde data sinyali kadar nokta icin linspace
% 
% figure(1);
% %Grafigi cizdir
% plot(t, oku);
% %baslik
% title("(THE AUDIO DATA) Sampling Frequency: " + fs + " Hz");
% %x ekseni tanimla
% xlabel('t(s)');
% %y ekseni tanimla
% ylabel('Amplitude');
% 
% %Adım2
% 
% %Negatif frekanslarin alinmamasi icin frekans eksenini ikiye bol
% 
% f = (0:fs/2-1);
% 
% %Dongu ilk ornekten son ornege kadardir.
% 
% tus_sayisi = 0; %Kayitli olan tus sayisi, baslangicta 0 
% sessiz = 1; %Dongu sessiz bolgede ise 1 degerindedir, degilse 0 degerindedir.
% sifir = 5; %Sintalin bittigine karar verilmesi icin okunmasi gereken 0 degeri.
% basla = 0; %Bir onceki tus sinyalinin baslama indeksi.
% figure(2)
% 
% for i=1:length(oku)
%     
%     if (sessiz)
%         
%        if (oku(i) ~= 0) %Dongu sessizde iken 0 olmayan bir deger okundugunda
%            sessiz = 0; %Artik dongu sessizde degil.
%            basla = i; %Tus sinyali baslangici su anki indekse ayarlandi.
%        end
%        
%     else
%         
%         if (oku(i) == 0) %Dongu tus sinyalinde 0 okunursa 
%             sifir = sifir - 1; %Sinyali bitirmek icin 0 degerini azalt
%             
%             if (sifir == 0) %5 tane 0 sinyali alindi mi diye bak
%                 %Bir sinyal tamamen alinmistir sinyali kaydederiz.
%                 
%                 tus_sayisi = tus_sayisi + 1; %Tus sayisi artir
%                 sifir = 5; %0 sayisi 
%                 sessiz = 1; %sessize al yine
%                 
%                 s = oku(basla:i-5); %Tus sinyali baslangictan baslanarak 
%                 %sondaki 5 sifiri atmak icin i-5 e kadardir
%                 
%                 fourier = abs(fft(s, fs)); %Fourier al
%                 fourier = fourier(1:fs/2); %Sinyal grafigini ortadan boleriz
%                 
%                 subplot(3, 3, tus_sayisi); %Figure de grafigin konumunu belirle
%                 plot(f, fourier); %Grafikleri cizdir
%                 
%                 %1000den fazla ve az maksimum genlikli frekanslari bulma
%                 fourier1 = find(fourier == max(fourier(1000:2000))); %1000den fazla icin
%                 fourier2 = find(fourier == max(fourier(1:1000))); %1000den az icin
%                 
%                 %bastaki 2 frekans degeri ile tiklanan tusu veren fonk()
%                 %fonksiyonu kullanilarak basilan tus ve frekans degerleri
%                 %bulundu
%           
%                 %grafiklere baslik yazdir
%                 title("Number: " + fonk(fourier1, fourier2) + ", freq1: " + fourier1 + " Hz, freq2: " + fourier2 + " Hz");
%                 %grafiklere x eksen adi yazdir
%                 xlabel("f (Hz)");
%                 %grafiklere y eksen adi yazdir
%                 ylabel("Magnitude");
%             end
%             
%         else
%             %Gerekli sayida 0 gelmeden sinyal geldi sessiz alana girmiyoruz
%             %zerosize eski degerine restore edilir.
%             sifir = 5;
%         end
%         
%     end
% end
% end

% function num = fonk(frekans1, frekans2)
% 
% % Bu fonk yazilirken dusuk ve yuksek freq matris haline getirilir
% % Girilen frekanslarin hangisinin dusuk frekansa ve hangisinin yuksek frekansa
% % yakin olduklari bulunur. 
% 
% dusuk = [697 770 852 941];% Dusuk frekans icin az olan dusuk
% %vektorundeki her elemandan tek tek cikarilarak farki en dusuk olanin
% %indeksini alindi
% 
% yuksek = [1209 1336 1477];% Yuksek frekans icin fazla olan yuksek
% %vektorundeki her elemandan tek tek cikarilarak farki en dusuk olanin
% %indeksini alindi
% 
% tus = [1 2 3;4 5 6; 7 8 9; 0 0 0];
% 
% a = abs(dusuk - min(frekans1,frekans2)); %negatif cikmamasi icin mutlak deger
% b = min(abs(dusuk - min(frekans1,frekans2)));
% dusukfreq = find((a==b), 1); % elde edilen indexe gore dusuk frekansi bul
% 
% c = abs(yuksek - max(frekans1,frekans2));
% d = min(abs(yuksek - max(frekans1,frekans2)));
% yuksekfreq = find((c==d), 1); % elde edilen indexe gore yuksek frekansi bul
% 
% num = tus(dusukfreq, yuksekfreq);%bulunan frekanslara gore hangi tusa basildigini sec
% 
% end

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
figure(1)
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
figure(2);
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

%Hande Yildirim, 141201047,Ele514,Odev2, Soru3

fs = 1000;%ornekleme zamani
% 2 saniye boyunca saniyede 1000 ornek alinacak sekilde zaman vektorunu
% belirle
zaman = 0:1/fs:2-1/fs;%zamani belirle
sinyal = chirp(zaman,0,2,100);%chirp komutuyla sinyali elde et
figure(3)
plot(sinyal)%sinyali cizdir
%baslik koy
title('the sound')
%eksenleri adlandir
xlabel('t(s)')
ylabel('Amplitude')
%grid ac
grid on
%pencere(window) isimli 500 örnekli sinyalleri oluşturma

save = zeros(250,61);%Fourierlerin kaydedilecegi matris
adim = 25; %adim sayisi
pencere = sinyal(1:500);
U = abs(fft(pencere,fs));

for i = 1:60
    % 61 adet donusum sinyali kaydedilecek hepsi 250 elemana sahip
    pencere = sinyal((1+(adim*i)):(500+(adim*i)));
    U = abs(fft(pencere)/length(pencere));
    save(:,i+1) = U(1:250)';  
    
end

y_axs = (0:length(pencere)/2-1) * (fs/length(pencere)); % y eksenini hesapla
x_axs = linspace(0,60*25/fs,61); % x eksenini esit parcalara bol

figure(4) %figure ac
contourf(x_axs,y_axs,save,'edgecolor','none');
%renklendirme cubugunu grafigin sagina cizdir
colorbar
%grid ac
grid on
%eksenleri adlandir
xlabel('Time(s)')
ylabel('Frequency(Hz)')


%Hande Yildirim, 141201047,Ele514,Odev2, Soru4

%Adim 1

%bilinmeyenler yazildi
simulas = 100;%simülasyon suresi
fs = 1000;%ornekleme freq
sig = 10; 
ro = 28; 
b = 8/3;
% baslangic koordinatlari
x0 = 0; 
y0 = 1; 
z0 = 20;

% bu asamada ilk olarak ode45 kullanarak denklem cozdurulur ve hatalar
% ayarlanir. Sonrasinda 1000 Hz icin zaman ayarlanip 1/1000'er 1/1000'er
% artis yapilir.

figure(5)
f = @(t, x) [sig*(x(2)-x(1)); x(1)*(ro-x(3))-x(2); x(1)*x(2)-b*x(3)];
sec = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
[t, xt] = ode45(f, 0:1/fs:simulas-1/fs, [x0 y0 z0], sec);


% Adim 2

% Daha sonra cozumler icin xt degiskenine 100000 elemanli olan 3 adet vektor
% verilip grafik cizdirildi.

plot3(xt(:, 1), xt(:, 2), xt(:, 3)); % baslangic noktalari belirlenir
hold on; % ikinci grafikte ayni figure ayni grafige cizilsin diye
plot3(xt(1,1),xt(1,2),xt(1,3),'r*');

%Iki cemberin orta noktalari sec

yaricap = 11;
% 11 yaricapli saydam daireler cizdir

[x,y,z] = sphere(50);
%x,y,z baslangic koordinatlar
x0 = 8.2; 
y0 = 8.2; 
z0 = 26.8;

%yaricapa gore x,y,z koordinatlar
x = x * yaricap + x0; 
y = y * yaricap + y0; 
z = z * yaricap + z0;
surf(x,y,z, 'FaceColor', 'black', 'LineStyle', 'none', 'FaceAlpha', 0.2)

[x,y,z] = sphere(50);
%x,y,z baslangic koordinatlar
x0 = -8.2; 
y0 = -8.2; 
z0 = 26.8;

%yaricapa gore x,y,z koordinatlar
x = x * yaricap + x0; 
y = y * yaricap + y0; 
z = z * yaricap + z0;
surf(x,y,z, 'FaceColor', 'red', 'LineStyle', 'none', 'FaceAlpha', 0.2)

%Grid ac
grid on;
%eksenleri adlandir
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');

hold off;


% Adim 3

f = (0: length(t) - 1) * (fs / length(t)); %Frekans ekseni

%grafikleri yeni figure cizdirme
figure(6);

subplot(3,2,1); %3x2'lik figure 1. grafigi
plot(t, xt(:,1)); % 1. vektor x1 
grid on; %grid ac
%eksenleri adlandir
xlabel('t (s)');
ylabel('x_{1}'); 
a = abs(fft(xt(:,1))); % sinyali normalize etme
fourierx = a / length(xt(:,1)); %fouriersini al

subplot(3,2,2); %3x2'lik figure 2. grafigi
plot((f(1:2000)), (fourierx(1:2000))); %Ilk 2000 ornegi cizdir
grid on;%grid ac
%eksenleri adlandir
xlabel('Freq(Hz)');
ylabel('|X_1(jw)|'); 

subplot(3,2,3); %3x2'lik figure 3. grafigi
plot(t, xt(:,2)); % 2. vektor x2 
grid on;%grid ac
%eksenleri adlandir
xlabel('t(s)');
ylabel('x_2');
b = abs(fft(xt(:,2)));
fouriery = b / length(xt(:,2)); % fouriersini al

subplot(3,2,4); % 3x2'lik figure 4. grafigi
plot((f(1:2000)) , (fouriery(1:2000))); % Ilk 2000 ornegi cizdir
grid on;% grid ac
%eksenleri adlandir
xlabel('Freq(Hz)');
ylabel('|X_2(jw)|');

subplot(3,2,5); % 3x2'lik figure 5. grafigi
plot(t,xt(:,3)); % 3. vektor x3
grid on;% grid ac
%eksenleri adlandir
xlabel('t (s)');
ylabel('x_3');
c = abs(fft(xt(:,3))); % normalize et
fourierz = c / length(xt(:,3)); % fouriersini al

subplot(3,2,6); % 3x2'lik figure 6. grafigi
plot( (f(1:2000)) , (fourierz(1:2000)) ); % Ilk 2000 ornegi cizdir
grid on;% grid ac
% eksenleri adlandir
xlabel('Freq (Hz)');
ylabel('|X_3(jw)|');

siyah1 = [8.2 8.2 26.8]; % siyah daire icin koordinat eksenleri
kirmizi1 = [-8.2 -8.2 26.8]; % kirmizi daire icin cizdir
yaricap = 11; % dairelerin yaricaplari

sin = xt'; % sinyalleri donustur
elemansayisi = size(xt); % sinyallerin eleman sayilarini bul
% 3 adet degiskenler kadar elamani olan  NaN vektoru 3 kez olusturulur
matris = NaN(3, elemansayisi(1), 3); % 3x10000 matris, her bir farkli renk icin.
renk = 3; %cizilecek renk 

% renk kodlari icin: 3 = maviyi, 2 = kirmiziyi, 1 = yesili temsil eder.

% Once 3 degisken kullanilarak vektorler tarandi ve noktalarin hangi cember icinde
% olduklari bulundu. 

% NaN'lar birlesmedigi icin renkler ayri ayri cizdirilince bir renk
% cizginin bittigi yerde digeri cizdirilmeye baslanir

for i = 1:elemansayisi(1) %Her elemani itere edecek.
    
    % renk degisimi olan nokta, kesinti olmamasi icin iki vektorede basilmalidir.
    
    matris(:, i, renk) = sin(:, i); % baslangicta bir onceki renk matrisine noktayi basildi
    
    if ((xt(i,1)-siyah1(1))^2+(xt(i,2)-siyah1(2))^2+(xt(i,3)-siyah1(3))^2 <= yaricap^2) %Siyah cember icinde mi?
        
        if (renk ~= 1) % bir onceki renk yesil degil
            
            renk = 1; % yesil yap
            matris(:, i,renk) = sin(:, i); % yesil vektore nokta basildi
            
        end
        
    elseif ((xt(i,1)-kirmizi1(1))^2+(xt(i,2)-kirmizi1(2))^2+(xt(i,3)-kirmizi1(3))^2 <= yaricap^2) % kirmizi cemberin icindeyse
        
        if (renk ~= 2) % bir onceki renk kirmizi degil
            
            renk = 2; % kirmizi yap
            matris(:,i,renk) = sin(:,i); % kirmizi vektore nokta basildi
            
        end
        
    else % iki cemberde de degil.
        
        if (renk ~= 3) % bir onceki renk mavi degil
            
            renk = 3; % mavi yap
            matris(:,i,renk) = sin(:,i); % mavi vektore nokta basildi
        
        end
        
    end
    
end

figure(7);

hold on;
plot3(xt(1,1),xt(1,2),xt(1,3),'r*'); %kirmizi ile baslama noktasıni * ile isaretle

hold on;
plot3(matris(1,:,1), matris(2,:,1), matris(3,:,1), 'Color', [0 1 0], 'LineWidth', 3); % yesil icin cizdir

hold on;
plot3(matris(1,:,2), matris(2,:,2), matris(3,:,2), 'Color', [1 0 0], 'LineWidth', 3); % kirmizi icin cizdir

hold on;
plot3(matris(1,:,3), matris(2,:,3), matris(3,:,3), 'Color', [0 0 1], 'LineWidth', 3); % mavi icin cizdir

%grid ac
grid on;

%eksenleri adlandir
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');

hold off;

figure(8);

for i=1:3 %Her degisken icin bir iterasyon
    
    subplot(3,1,i);
    
    hold on;
    plot(t', matris(i,:,1), 'Color', [0 1 0], 'LineWidth', 3); % yesil icin cizdir
    
    hold on
    plot(t', matris(i,:,2), 'Color', [1 0 0], 'LineWidth', 3); % kirmizi icin cizdir
    
    hold on
    plot(t', matris(i,:,3), 'Color', [0 0 1], 'LineWidth', 3); % mavi icin cizdir
    
    %grid ac
    grid on;
    
    %eksenleri adlandir
    xlabel('t(s)');
    ylabel("x_{" + i + "_}");
   
    hold off;
    hold off;
    hold off;
    
end

%Adim 4

%Bu kisim yorum olarak yazilmistir.



% 1. Sorunun fonksiyonu

function soru(dos)
[oku, fs] = audioread(dos);%Ornekleme frekansi
t = linspace(0, length(oku)/fs, length(oku));
%Sinyal zamani sampleSayisi/fs
%t ekseninde data sinyali kadar nokta icin linspace

figure(9);
%Grafigi cizdir
plot(t, oku);
%baslik
title("(THE AUDIO DATA) Sampling Frequency: " + fs + " Hz");
%x ekseni tanimla
xlabel('t(s)');
%y ekseni tanimla
ylabel('Amplitude');

%Adım2

%Negatif frekanslarin alinmamasi icin frekans eksenini ikiye bol

f = (0:fs/2-1);

%Dongu ilk ornekten son ornege kadardir.

tus_sayisi = 0; %Kayitli olan tus sayisi, baslangicta 0 
sessiz = 1; %Dongu sessiz bolgede ise 1 degerindedir, degilse 0 degerindedir.
sifir = 5; %Sintalin bittigine karar verilmesi icin okunmasi gereken 0 degeri.
basla = 0; %Bir onceki tus sinyalinin baslama indeksi.
figure(10)

for i=1:length(oku)
    
    if (sessiz)
        
       if (oku(i) ~= 0) %Dongu sessizde iken 0 olmayan bir deger okundugunda
           sessiz = 0; %Artik dongu sessizde degil.
           basla = i; %Tus sinyali baslangici su anki indekse ayarlandi.
       end
       
    else
        
        if (oku(i) == 0) %Dongu tus sinyalinde 0 okunursa 
            sifir = sifir - 1; %Sinyali bitirmek icin 0 degerini azalt
            
            if (sifir == 0) %5 tane 0 sinyali alindi mi diye bak
                %Bir sinyal tamamen alinmistir sinyali kaydederiz.
                
                tus_sayisi = tus_sayisi + 1; %Tus sayisi artir
                sifir = 5; %0 sayisi 
                sessiz = 1; %sessize al yine
                
                s = oku(basla:i-5); %Tus sinyali baslangictan baslanarak 
                %sondaki 5 sifiri atmak icin i-5 e kadardir
                
                fourier = abs(fft(s, fs)); %Fourier al
                fourier = fourier(1:fs/2); %Sinyal grafigini ortadan boleriz
                
                subplot(3, 3, tus_sayisi); %Figure de grafigin konumunu belirle
                plot(f, fourier); %Grafikleri cizdir
                
                %1000den fazla ve az maksimum genlikli frekanslari bulma
                fourier1 = find(fourier == max(fourier(1000:2000))); %1000den fazla icin
                fourier2 = find(fourier == max(fourier(1:1000))); %1000den az icin
                
                %bastaki 2 frekans degeri ile tiklanan tusu veren fonk()
                %fonksiyonu kullanilarak basilan tus ve frekans degerleri
                %bulundu
          
                %grafiklere baslik yazdir
                title("Number: " + fonk(fourier1, fourier2) + ", freq1: " + fourier1 + " Hz, freq2: " + fourier2 + " Hz");
                %grafiklere x eksen adi yazdir
                xlabel("f (Hz)");
                %grafiklere y eksen adi yazdir
                ylabel("Magnitude");
            end
            
        else
            %Gerekli sayida 0 gelmeden sinyal geldi sessiz alana girmiyoruz
            %zerosize eski degerine restore edilir.
            sifir = 5;
        end
        
    end
end

end

function num = fonk(frekans1, frekans2)

% Bu fonk yazilirken dusuk ve yuksek freq matris haline getirilir
% Girilen frekanslarin hangisinin dusuk frekansa ve hangisinin yuksek frekansa
% yakin olduklari bulunur. 

dusuk = [697 770 852 941];% Dusuk frekans icin az olan dusuk
%vektorundeki her elemandan tek tek cikarilarak farki en dusuk olanin
%indeksini alindi

yuksek = [1209 1336 1477];% Yuksek frekans icin fazla olan yuksek
%vektorundeki her elemandan tek tek cikarilarak farki en dusuk olanin
%indeksini alindi

tus = [1 2 3;4 5 6; 7 8 9; 0 0 0];

a = abs(dusuk - min(frekans1,frekans2)); %negatif cikmamasi icin mutlak deger
b = min(abs(dusuk - min(frekans1,frekans2)));
dusukfreq = find((a==b), 1); % elde edilen indexe gore dusuk frekansi bul

c = abs(yuksek - max(frekans1,frekans2));
d = min(abs(yuksek - max(frekans1,frekans2)));
yuksekfreq = find((c==d), 1); % elde edilen indexe gore yuksek frekansi bul

num = tus(dusukfreq, yuksekfreq);%bulunan frekanslara gore hangi tusa basildigini sec

end