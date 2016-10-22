% Simulación acústica de una Sala.
%% Parámetros de entrada:
% -------------------------------------------------------------------------
clear all;
profundidad = 20;       %[m]
anchura = 20;           %[m]
altura = 4;             %[m]
superficie_puerta = 3;  %[m^2]
superficie_v1 = 4;      %[m^2]
superficie_v2 = 4;      %[m^2]
superficie_but = 1;     %[m^2]
ro = 1.29;
numero_but = 10;
temperatura = 24;
humedad = 10;
potencia_fuente = 10^-5;
P0 = 2*10^-5;
f_directividad = 2;
emisor = [1,2,1];
receptor = [1,9,1];
% Matriz coeficientes de absorcion:
coeficientes_absorcion = zeros(10:6);
% Cada columna corresponde a una banda de frecuencias (125, 250, 500, 1000,
% 2000, 4000)Hz
% Techo
coeficientes_absorcion(1,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Suelo
coeficientes_absorcion(2,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Pared1
coeficientes_absorcion(3,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Pared2
coeficientes_absorcion(4,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Pared3
coeficientes_absorcion(5,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Pared4
coeficientes_absorcion(6,:) = [0.30,0.45,0.30,0.25,0.40,0.25];
% Ventana1
coeficientes_absorcion(7,:) = [0.035,0.04,0.027,0.03,0.02,0.02];
% Ventana2
coeficientes_absorcion(8,:) = [0.035,0.04,0.027,0.03,0.02,0.02];
% Puerta
coeficientes_absorcion(9,:) = [0.05,0.16,0.13,0.01,0.06,0.05];
% Butacas
coeficientes_absorcion(10,:) = [0.02,0.02,0.03,0.035,0.038,0.038];
% -------------------------------------------------------------------------

%%Cálculo de valores:
% Matriz de superficies:
s = zeros(1,10);
s = [profundidad*anchura,profundidad*anchura,altura*anchura-superficie_v1-superficie_v2,altura*profundidad,altura*anchura-superficie_puerta,altura*profundidad,superficie_puerta,superficie_v1,superficie_v2,superficie_but];

% Velocidad de propagacion dependiendo de la temperatura:
c = 331.3*sqrt(1+temperatura/273.15);

% Calculo de la superficie total de la sala:
superficie_total = 0;
for i=1 : 10
    superficie_total = superficie_total + s(i);
end

% Absorcion media para cada banda de frecuencias y constante de sala R:
absorcion_media = zeros(1,6);
numerador = zeros(1,6);
denominador = zeros(1,6);
R = zeros(1,6);
for j=1 : 6
    for i=1 : 10
        numerador(j) = numerador(j) + s(i)*coeficientes_absorcion(i,j);
        denominador(j) = denominador(j) + s(i);
    end
    absorcion_media(j) = numerador(j)/denominador(j);
end
for j=1 : 6
    R(j) = superficie_total*absorcion_media(j)/(1-absorcion_media(j));
end

% Calculamos el RT60 para cada banda de frecuencias:
volumen = altura*anchura*profundidad;
RT60 = zeros(1,6);
for i=1 : 6
    RT60(i) = 60*volumen / (1.086*superficie_total*c*absorcion_media(i));
end


%% Ejemplo práctico con un audio:
% Vector de tiempo:
[y,fs] = audioread('Emisor2.wav');
y = y(:,1)';
t = 0 : 1/fs : 2;

% Creamos un señal de ruido:
noise = wgn(1,fs,1);
% Calculamos la recta de decaimiento:
decaimiento = 1 - ( 1/mean(RT60) ) * t;
posicion = 1;
min = abs(decaimiento(1));
for i=2 : length(decaimiento)
    if abs(decaimiento(i)) < min
        min = abs(decaimiento(i));
        posicion = i;
    end
end
decaimiento = decaimiento(1:fs).*noise;
decaimiento = decaimiento(1:posicion);
y_reverb = conv(y,decaimiento);

% Normalizamos:
y_reverb = y_reverb / max(y_reverb);
% Audio que simula lo que escucha el receptor:
audiowrite('Receptor.wav',y_reverb,fs);

Y = fft(y);
Y_reverb = fft(y_reverb);

% Nivel de presion sonora campo reverberante y directo:
SPL_rev = zeros(1,6);
SPL_dir = 0;

% Primero calculamos la distancia crítica y la distancia entre emisor/receptor.
distancia = sqrt((emisor(1)-receptor(1))^2 + (emisor(2)-receptor(2))^2 + (emisor(3)-receptor(3))^2);
distancia_critica = sqrt(f_directividad*mean(R)/(16*pi));
SPL_max = 10*log10(potencia_fuente/(10^-12));
for i=1 : 6
    SPL_rev(i) = 10*log10( (potencia_fuente*ro*c*(4/R(i)))/(P0^2) );
end
SPL_dir = 10*log10( (potencia_fuente*ro*c*(f_directividad/(4*pi*distancia*distancia)))/(P0^2) );

if distancia >= distancia_critica
    y_reverb = y_reverb*mean(SPL_rev)/SPL_max;
else
    y_reverb = y_reverb*SPL_dir/SPL_max;
end

sound(y_reverb,fs);
subplot(2,2,1);
plot(y);
subplot(2,2,2);
plot(y_reverb);
subplot(2,2,3);
plot(abs(Y));
subplot(2,2,4);
plot(abs(Y_reverb));