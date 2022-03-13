%gr. 1., projekt 1.,PSIO
%%

%I. Wybrane transformacje opisów dynamicznych
L=[1 0 0 2 0]; %wspó³czynniki licznika
M=[1 1 2 4 12 16]; %wspó³czynniki mianownika
H=tf(L,M) %opis transmitancyjny
%a)przejœcie do opisu stanowego
[A,B,C,D]=ssdata(H)
%b)transformacja opisu transmisyjnego w zero-biegunowy
[z1,p1,k1]=tf2zp(L,M)
H_zpk=zpk(z1,p1,k1);
%c)rozk³ad na u³amki proste
[r2, p2, k2] = residue(L,M)
%d)postaæ bikwadratowa
[Hbk,g]=tf2sos(L,M)
%%
% II Wyznaczyæ odpowiedŸ impulsow¹ oraz skokow¹ dla dwóch filtrów dyskretnych
%Wyznaczenie transmitancji
LSOI=[1 1 -2 1]; %licznik transmitancji filtru SOI
MSOI=[1 0 0 0 0]; %mianownik transmitancji filtru SOI
LNOI=[0.12 0.8 0.8 0.12]; %licznik transmitancji filtru NOI
MNOI=[1 -0.4 0.8 0.2]; %mianownik transmitancji filtru NOI
%Odpowiedzi impulsowe
iSOI=dimpulse(LSOI,MSOI,20); %odpowiedŸ impulsowa filtru SOI
iNOI=dimpulse(LNOI,MNOI,30); %odpowiedŸ impulsowa filtru NOI
figure(1)
subplot(2,1,1); stem(iSOI); grid;
title ('OdpowiedŸ impulsowa filtru dyskretnego SOI');
xlabel('Próbki');ylabel('Amplituda');
subplot(2,1,2); stem(iNOI); grid;
title ('OdpowiedŸ impulsowa filtru dyskretnego NOI');
xlabel('Próbki');ylabel('Amplituda');
%Odpowiedzi skokowe
sSOI=dstep(LSOI,MSOI,20); %odpowiedŸ skokowa filtru NOI
sNOI=dstep(LNOI,MNOI,30); %odpowiedŸ skokowa filtru SOI
figure(2)
subplot(2,1,1); stem(sSOI); grid;
title ('OdpowiedŸ skokowa filtru dyskretnego SOI');
xlabel('Próbki');ylabel('Amplituda');
subplot(2,1,2); stem(sNOI); grid;
title ('OdpowiedŸ skokowa filtru dyskretnego NOI');
xlabel('Próbki');ylabel('Amplituda');
%%
% III Filtracja sygna³u
%1)Wyznaczenie parametrów
N=3; %litera B
I=13; %litera J
%Wartoœci amplitud
A1=1/N; %1/3
A2=1/I; %1/13
A3=1; 
A4=0.75;
%Wartoœci czêstotliwoœci [Hz]
F1=50+N; %53 Hz
F2=80+I; %93 Hz
F3=105+2*N; %111 Hz
F4=160+2*(N+I); %192 Hz
%Czas
Fs = 750; % czestotliwosc probkowania
T = 1/Fs; % okres probkowania
L = 500; % dlugosc sygnalu(liczba probek)
t = (0:L-1)*T; % wektor czasu
%Postaæ sygna³u
x1=A1*sin(2*pi*F1*t)+A2*sin(2*pi*F2*t)+A3*sin(2*pi*F3*t)+A4*sin(2*pi*F4*t);
figure(3)
plot(t,x1)
title('Przebieg sygna³u oryginalnego');
xlabel('Czas [s]');ylabel('Amplituda');
%2)zak³ócenie
xz=2/(N+I)*randn(size(t));
x11=x1+xz;
figure(4)
subplot(2,1,1)
plot(t,x11)
title('Przebieg sygna³u zak³óconego');
xlabel('Czas [s]');ylabel('Amplituda');
subplot(2,1,2)
plot(t,xz)
title('Zak³ócenie');
xlabel('Czas [s]');ylabel('Amplituda');
%3)porównanie
figure(5)
plot(t,x1,t,x11,'r')
legend('sygna³ oryginalny','sygna³ zak³ócony')
xlabel('Czas [s]');ylabel('Amplituda');
%% 4) analiza widmowa
%sygna³ oryginlany
NFFT = 2^nextpow2(L); %liczba punktów transformaty bêd¹ca potêg¹ liczby 2
Y = fft(x1, NFFT)/L; % generujemy NFFT-punktowa dyskretna transformate Fouriera
f = (Fs/2)*linspace(0, 1, NFFT/2+1); % tworzymy wektor czestotliwosci, skladajacy sie z NFFT/2+1 wartosci
% pomiedzy wartosciami 0 a Fs/2
figure(6)
plot(f, 2*abs(Y(1:NFFT/2+1))); % kreslimy wykres widma, pomijajac druga jego polowe
% kreslony jest modul; mnozenie przez 2 ma za zadanie uwzglednic fakt,
% ze pominelismy jego polowe
title('Charakterystyka amplitudowa sygna³u oryginalnego')
xlabel('Czestotliwoœæ [Hz]');ylabel('|Y(f)|');
%sygna³ zak³ócony
Y1 = fft(x11, NFFT)/L;
figure(7)
plot(f, 2*abs(Y1(1:NFFT/2+1)));
title('Charakterystyka amplitudowa sygna³u zak³óconego')
xlabel('Czestotliwoœæ [Hz]');ylabel('|Y(f)|');
%% 5) okna czasowe
%% trójk¹tne
%wvtool(triang(L))
%% okno Hanninga
%wvtool(hanning(L))
%% okno Blackmana
%wvtool(blackman(L))
%% 6) Projektowanie filtru
% a)filtr eliptyczny 10 rzêdu, przepuszcza F2, blokuje F1,F3,F4
[b,a] = ellip(5,0.1,40,[60 97]*2/Fs); %filtr eliptyczny
[H,W] = freqz(b,a,512);
figure(8)
plot(W*Fs/(2*pi),unwrap(angle(H)));
title('Charakterystyka fazowa filtru eliptycznego');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Faza[Rad]'); grid;
figure(9)
plot(W*Fs/(2*pi),abs(H));
title('Charakterystyka amplitudowa filtru eliptycznego');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Amplituda'); 
axis([0 375 0 1.2])
grid;
spfe = filter(b,a,x1);
figure(10)
plot(t,x1,t,spfe);
xlabel('Czas (sekundy)');
ylabel('Amplituda')
title('Przebieg czasowy sygna³u przed i po filtracji'); 
legend({'Przed filtracj¹','Po filtracji'})
S = fft(x1,512)/L;
SF = fft(spfe,512)/L;
w = (0:255)/256*(Fs/2);
figure(11)
plot(w,abs([2*S(1:256)' 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna³u przed i po filtracji');
xlabel('Czêstotliwoœæ (Hz)');
ylabel('Modu³ transformaty Fouriera');grid;
legend({'Przed filtracj¹','Po filtracji'})
%b) filtry Czebyszewa
%% czebyszew 1
[b1,a1] = cheby1(6,0.1,80*2/Fs,'high');
[H1,W1]= freqz(b1,a1,512);
figure(12)
plot(W1*Fs/(2*pi),unwrap(angle(H1)));
title('Charakterystyka fazowa filtru Czebyszewa I rodzaju');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Faza[Rad]'); 
grid;
figure(13)
plot(W1*Fs/(2*pi),abs(H1));
title('Charakterystyka amplitudowa filtru Czebyszewa I rodzaju');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Amplituda'); 
grid;
figure(14)
spfc1 = filter(b1,a1,x11);
plot(t,x11,t,spfc1);
xlabel('Czas (sekundy)'); 
ylabel('Amplituda') 
title('Przebieg czasowy sygna³u przed i po filtracji'); 
legend({'Przed filtracj¹','Po filtracji'})
S = fft(x11,512)/L;
SF = fft(spfc1,512)/L;
w = (0:255)/256*(Fs/2);
figure(15)
plot(w,abs([2*S(1:256)' 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna³u przed i po filtracji');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Modu³ transformaty Fouriera'); grid; 
legend({'Przed filtracj¹','Po filtracji'})
%% czebyszew II
[b2,a2] = cheby2(4,40,40*2/Fs,'high');
[H2,W2]= freqz(b2,a2,512);
figure(16)
plot(W2*Fs/(2*pi),unwrap(angle(H2)));
title('Charakterystyka fazowa filtru Czebyszewa II rodzaju');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Faza[Rad]'); 
grid;
figure(17)
plot(W2*Fs/(2*pi),abs(H2));
title('Charakterystyka amplitudowa filtru Czebyszewa II rodzaju');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Amplituda'); 
grid;
spfc2 = filter(b2,a2,x11);
figure(18)
plot(t,x11,t,spfc2);
xlabel('Czas (sekundy)'); 
ylabel('Amplituda') 
title('Przebieg czasowy sygna³u przed i po filtracji'); 
legend({'Przed filtracj¹','Po filtracji'})
S = fft(x11,512)/L;
SF = fft(spfc2,512)/L;
w = (0:255)/256*(Fs/2);
figure(19)
plot(w,abs([2*S(1:256)', 2*SF(1:256)']));
title('Charakterystyka amplitudowa sygna³u przed i po filtracji');
xlabel('Czêstotliwoœæ (Hz)');  
ylabel('Modu³ transformaty Fouriera'); grid; 
legend({'Przed filtracj¹','Po filtracji'});
