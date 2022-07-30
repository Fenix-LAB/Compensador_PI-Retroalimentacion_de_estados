%Matriz de un PI
% A                                 B              
%     [ (A - BK)     BKe ]              [ 0 ]
%     [    -C         0  ]              [ 1 ]

% C
%     [ C  0 ]

%Ejempllo de un PI para 3er orden
%API=[x      x      x   y;
%     x      x      x   y;
%     x-K1  x-K2  x-K3  Ke;
%     z      z      z   0];

%BPI=[0;
%     0;
%     0;
%     1];

%CPI=[z z z 0];

clc;
%MAtrices del sistema original
'Compensador PI'
A=[-5 1 0;0 -2 1;0 0 -1;];
B=[0;0;1;];
C=[-1 1 0];

%Revisaoms el comportamiento del sistema y si es un sistema controlable
'Sistema'
[NUM,DEN]=ss2tf(A,B,C,0);
T=tf(NUM,DEN)
step(T)

Co=ctrb(A,B)
Rango=rank(Co)

%Variables simbolicas que usaremo en el programa
syms s K1 K2 K3 Ke P 

%Resultados necesarios para formar la matriz del PI
K=[K1 K2 K3];
Ap=A-B*K
Bp=B*Ke
Cp=-C

%Con los valores obtenidos formamos la matriz PI
'Matriz A del sistema con PI'
API=[-5 1 0 0;0 -2 1 0;-K1 -K2 -1-K3 Ke;1 -1 0 0]
BPI=[0;0;0;1];
CPI=[-1 1 0 0];
sI=[s 0 0 0;0 s 0 0;0 0 s 0;0 0 0 s];
r=sI-API;
'Ecuacion caracteristica del sistema'
De=det(r)
%Este denominador es la ecuaicion caracteristica del sistema y la vamos a
%comparar con la ecuacion deseada

%Para la ecuacion deseada debemos realizar los siguiente:

%Se obtienen los polos de un sitema de 2do orden que cumplan con los
%requisitos solicitados (OS y Ts)
OS=20.8;
Ts=4; 
z=(-log(OS/100))/(sqrt(pi^2+log(OS/100)^2)); %Factor de amortiguamiento
wn=4/(z*Ts);                                 %frecuancia natural
%VEmos el sistema y anotamos el denominador
[num,den]=ord2(wn,z);
'Polonimio de los polos dominantes'
den

%Polos dominantes
'Polos dominantes'
r=roots(den) 
%Para est ejemplo tenemos -1 +- 2i

%Para un Pi de un ssistema de 3er orden necesitamos 4 polos, 2 son de los
%polos dominantes y 2 los debemos proponer.

%2 Polos dominantes, 1 que proponemos y un ultimo polo el cual vamos a
%poder modificar, lllamado P

%Para este ejemplo Usaremos los 2 polos, un polo en -4 para eliminar al
%cero y el polo P

% Ec = Ecuacion caracteristica (Deseada)
% Pc = Polo ue cancela al cero

Pc=4;
Ec=(s+P)*(s+Pc)*(s^2+2*s+5)
'Ecuacion caracteristica deseada'
Ec2=expand(Ec)

%VEmos los coeficiented de la ecuacion caracteristica el sistema y la
%ecuacion caracteristica deaseada

CofS=coeffs(De,s)  %Coeficientes del sistema
CofD=coeffs(Ec2,s) %Coeficientes deseados

%Aqui se plantean las ecuaciones para encontrar los valores de K1, K2, K3 y
%Ke

%K3 + 8 = P + 6
%K2 + 7*K3 + 17 = 6*P + 13
%K1 + 5*K2 + 10*K3 + Ke + 10 = 13*P + 20
%4*Ke = 20*P

%Resolvemos el sistema de acuaciones
'Solucion de las ecuaciones caracteristicas'
Sol=CofS==CofD
'Vallores de Ke, K1, K2y K3 en funcion de un tercer polo en P'
[Ken, K1n, K2n, K3n] = solve(Sol,Ke, K1, K2, K3)

%Es neceasrio un valor de P, se recomienda 3 veces la parte real de los
%polos dominantes¨
Kenw=subs(Ken, P, 100);
K1nw=subs(K1n, P, 100);
K2nw=subs(K2n, P, 100);
K3nw=subs(K3n, P, 100);
%P=100;

%Ahora usamos estos valores para terminar la matris del PI
'Sistema compenado con PI'
Acom=[-5 1 0 0;0 -2 1 0;-K1nw -K2nw -1-K3nw Kenw;1 -1 0 0]
BPI
CPI

%Analizamos el resultado
%[numG,denG]=ss2tf(Acom,BPI,CPI,0);
%G=tf(numG,denG)
%step(G)

