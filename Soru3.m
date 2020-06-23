%Hande Yildirim, 141201047,Ele514,Odev2, Soru3

fs = 1000;%ornekleme zamani
% 2 saniye boyunca saniyede 1000 ornek alinacak sekilde zaman vektorunu
% belirle
zaman = 0:1/fs:2-1/fs;%zamani belirle
sinyal = chirp(zaman,0,2,100);%chirp komutuyla sinyali elde et
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

figure() %figure ac
contourf(x_axs,y_axs,save,'edgecolor','none');
%renklendirme cubugunu grafigin sagina cizdir
colorbar
%grid ac
grid on
%eksenleri adlandir
xlabel('Time(s)')
ylabel('Frequency(Hz)')