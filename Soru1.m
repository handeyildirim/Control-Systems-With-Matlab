%Hande Yildirim, 141201047,Ele514,Odev2, Soru1

%SORU1

%Adim 1
soru("sound01.wav"); %3 adet dosya icin kodu calistirdik 

function soru(dos)
[oku, fs] = audioread(dos);%Ornekleme frekansi
t = linspace(0, length(oku)/fs, length(oku));
%Sinyal zamani sampleSayisi/fs
%t ekseninde data sinyali kadar nokta icin linspace

figure(1);
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
figure(2)

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


