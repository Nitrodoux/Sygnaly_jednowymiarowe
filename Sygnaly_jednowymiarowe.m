%gr. 1., projekt 1.,PSIO
%%

%I. Wybrane transformacje opis�w dynamicznych
L=[1 0 0 2 0]; %wsp�czynniki licznika
M=[1 1 2 4 12 16]; %wsp�czynniki mianownika
H=tf(L,M) %opis transmitancyjny
%a)przej�cie do opisu stanowego
[A,B,C,D]=ssdata(H)
%b)transformacja opisu transmisyjnego w zero-biegunowy
[z1,p1,k1]=tf2zp(L,M)
H_zpk=zpk(z1,p1,k1);
%c)rozk�ad na u�amki proste
[r2, p2, k2] = residue(L,M)
%d)posta� bikwadratowa
[Hbk,g]=tf2sos(L,M)
%%
% II Wyznaczy� odpowied� impulsow� oraz skokow� dla dw�ch filtr�w dyskretnych
%Wyznaczenie transmitancji
LSOI=[1 1 -2 1]; %licznik transmitancji filtru SOI
MSOI=[1 0 0 0 0]; %mianownik transmitancji filtru SOI
LNOI=[0.12 0.8 0.8 0.12]; %licznik transmitancji filtru NOI
MNOI=[1 -0.4 0.8 0.2]; %mianownik transmitancji filtru NOI
%Odpowiedzi impulsowe
iSOI=dimpulse(LSOI,MSOI,20); %odpowied� impulsowa filtru SOI
iNOI=dimpulse(LNOI,MNOI,30); %odpowied� impulsowa filtru NOI
figure(1)
subplot(2,1,1); stem(iSOI); grid;
title ('Odpowied� impulsowa filtru dyskretnego SOI');
xlabel('Pr�bki');ylabel('Amplituda');
subplot(2,1,2); stem(iNOI); grid;
title ('Odpowied� impulsowa filtru dyskretnego NOI');
xlabel('Pr�bki');ylabel('Amplituda');
%Odpowiedzi skokowe
sSOI=dstep(LSOI,MSOI,20); %odpowied� skokowa filtru NOI
sNOI=dstep(LNOI,MNOI,30); %odpowied� skokowa filtru SOI
figure(2)
subplot(2,1,1); stem(sSOI); grid;
title ('Odpowied� skokowa filtru dyskretnego SOI');
xlabel('Pr�bki');ylabel('Amplituda');
subplot(2,1,2); stem(sNOI); grid;
title ('Odpowied� skokowa filtru dyskretnego NOI');
xlabel('Pr�bki');ylabel('Amplituda');
%%
% III Filtracja sygna�u
%1)Wyznaczenie parametr�w
N=3; %litera B
I=13; %litera J
%Warto�ci amplitud
A1=1/N; %1/3
A2=1/I; %1/13
A3=1; 
A4=0.75;
%Warto�ci cz�stotliwo�ci [Hz]
F1=50+N; %53 Hz
F2=80+I; %93 Hz
F3=105+2*N; %111 Hz
F4=160+2*(N+I); %192 Hz
%Czas
Fs = 750; % czestotliwosc probkowania
T = 1/Fs; % okres probkowania
L = 500; % dlugosc sygnalu(liczba probek)
t = (0:L-1)*T; % wektor czasu
%Posta� sygna�u
x1=A1*sin(2*pi*F1*t)+A2*sin(2*pi*F2*t)+A3*sin(2*pi*F3*t)+A4*sin(2*pi*F4*t);
figure(3)
plot(t,x1)
title('Przebieg sygna�u oryginalnego');
xlabel('Czas [s]');ylabel('Amplituda');
%2)zak��cenie
xz=2/(N+I)*randn(size(t));
x11=x1+xz;
figure(4)
subplot(2,1,1)
plot(t,x11)
title('Przebieg sygna�u zak��conego');
xlabel('Czas [s]');ylabel('Amplituda');
subplot(2,1,2)
plot(t,xz)
title('Zak��cenie');
xlabel('Czas [s]');ylabel('Amplituda');
%3)por�wnanie
figure(5)
plot(t,x1,t,x11,'r')
legend('sygna� oryginalny','sygna� zak��cony')
xlabel('Czas [s]');ylabel('Amplituda');
%% 4) analiza widmowa
%sygna� oryginlany
NFFT = 2^nextpow2(L); %liczba punkt�w transformaty b�d�ca pot�g� liczby 2
Y = fft(x1, NFFT)/L; % generujemy NFFT-punktowa dyskretna transformate Fouriera
f = (Fs/2)*linspace(0, 1, NFFT/2+1); % tworzymy wektor czestotliwosci, skladajacy sie z NFFT/2+1 wartosci
% pomiedzy wartosciami 0 a Fs/2
figure(6)
plot(f, 2*abs(Y(1:NFFT/2+1))); % kreslimy wykres widma, pomijajac druga jego polowe
% kreslony jest modul; mnozenie przez 2 ma za zadanie uwzglednic fakt,
% ze pominelismy jego polowe
title('Charakterystyka amplitudowa sygna�u oryginalnego')
xlabel('Czestotliwo�� [Hz]');ylabel('|Y(f)|');
%sygna� zak��cony
Y1 = fft(x11, NFFT)/L;
figure(7)
plot(f, 2*abs(Y1(1:NFFT/2+1)));
title('Charakterystyka amplitudowa sygna�u zak��conego')
xlabel('Czestotliwo�� [Hz]');ylabel('|Y(f)|');
%% 5) okna czasowe
%% tr�jk�tne
%wvtool(triang(L))
%% okno Hanninga
%wvtool(hanning(L))
%% okno Blackmana
%wvtool(blackman(L))
%% 6) Projektowanie filtru
% a)filtr eliptyczny 10 rz�du, przepuszcza F2, blokuje F1,F3,F4
[b,a] = ellip(5,0.1,40,[60 97]*2/Fs); %filtr eliptyczny
[H,W] = freqz(b,a,512);
figure(8)
plot(W*Fs/(2*pi),unwrap(angle(H)));
title('Charakterystyka fazowa filtru eliptycznego');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Faza[Rad]'); grid;
figure(9)
plot(W*Fs/(2*pi),abs(H));
title('Charakterystyka amplitudowa filtru eliptycznego');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Amplituda'); 
axis([0 375 0 1.2])
grid;
spfe = filter(b,a,x1);
figure(10)
plot(t,x1,t,spfe);
xlabel('Czas (sekundy)');
ylabel('Amplituda')
title('Przebieg czasowy sygna�u przed i po filtracji'); 
legend({'Przed filtracj�','Po filtracji'})
S = fft(x1,512)/L;
SF = fft(spfe,512)/L;
w = (0:255)/256*(Fs/2);
figure(11)
plot(w,abs([2*S(1:256)' 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna�u przed i po filtracji');
xlabel('Cz�stotliwo�� (Hz)');
ylabel('Modu� transformaty Fouriera');grid;
legend({'Przed filtracj�','Po filtracji'})
%b) filtry Czebyszewa
%% czebyszew 1
[b1,a1] = cheby1(6,0.1,80*2/Fs,'high');
[H1,W1]= freqz(b1,a1,512);
figure(12)
plot(W1*Fs/(2*pi),unwrap(angle(H1)));
title('Charakterystyka fazowa filtru Czebyszewa I rodzaju');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Faza[Rad]'); 
grid;
figure(13)
plot(W1*Fs/(2*pi),abs(H1));
title('Charakterystyka amplitudowa filtru Czebyszewa I rodzaju');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Amplituda'); 
grid;
figure(14)
spfc1 = filter(b1,a1,x11);
plot(t,x11,t,spfc1);
xlabel('Czas (sekundy)'); 
ylabel('Amplituda') 
title('Przebieg czasowy sygna�u przed i po filtracji'); 
legend({'Przed filtracj�','Po filtracji'})
S = fft(x11,512)/L;
SF = fft(spfc1,512)/L;
w = (0:255)/256*(Fs/2);
figure(15)
plot(w,abs([2*S(1:256)' 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna�u przed i po filtracji');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Modu� transformaty Fouriera'); grid; 
legend({'Przed filtracj�','Po filtracji'})
%% czebyszew II
[b2,a2] = cheby2(4,40,40*2/Fs,'high');
[H2,W2]= freqz(b2,a2,512);
figure(16)
plot(W2*Fs/(2*pi),unwrap(angle(H2)));
title('Charakterystyka fazowa filtru Czebyszewa II rodzaju');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Faza[Rad]'); 
grid;
figure(17)
plot(W2*Fs/(2*pi),abs(H2));
title('Charakterystyka amplitudowa filtru Czebyszewa II rodzaju');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Amplituda'); 
grid;
spfc2 = filter(b2,a2,x11);
figure(18)
plot(t,x11,t,spfc2);
xlabel('Czas (sekundy)'); 
ylabel('Amplituda') 
title('Przebieg czasowy sygna�u przed i po filtracji'); 
legend({'Przed filtracj�','Po filtracji'})
S = fft(x11,512)/L;
SF = fft(spfc2,512)/L;
w = (0:255)/256*(Fs/2);
figure(19)
plot(w,abs([2*S(1:256)', 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna�u przed i po filtracji');
xlabel('Cz�stotliwo�� (Hz)');  
ylabel('Modu� transformaty Fouriera'); grid; 
legend({'Przed filtracj�','Po filtracji'});
