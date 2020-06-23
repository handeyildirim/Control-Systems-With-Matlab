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

figure(1)
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
figure(2);

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

figure(3);

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

figure(4);

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