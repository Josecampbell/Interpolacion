<<<<<<< HEAD
function[]=Diferencias_divididas_interpolante_newton()
clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Entrada de Datos--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=[1.0 1.3 1.6 1.9 2.2];%nodos del eje x
yn=[0.7651977 0.6200860 0.4554022 0.2818186 0.1103623];%imagenes de los nodos
F=Diferencias_divididas_newton(xn,yn);%llamada al metodo interpolante de Newton
P=Polinomio(F,xn);%la llamada a esta funci�n es para construir el polinomio
P=matlabFunction(P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Parametros de salidad-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F es una matrix digonal inferior cuyas entradas son los polinomios de las iteraciones
% F=diferencias_progresivas(x,y);
x0=1.5;%punto de prueba
y0=P(x0);%Imagen de punto de Prueba
close all
alfa=xn(length(xn))-xn(length(xn)-1);%extremo izquierdo del intervalo donde se quiere visualizar el estudio
alfa1=xn(length(xn))+xn(length(xn)-1);%extremo derecho del intervalo donde se quiere visualizar el estudio
hold on;fplot(P,[alfa,alfa1]);scatter(xn,yn);scatter(x0,y0);legend('Polinomio de Newton','nodos interpolantes','Nodo evaludo en el polinomio');hold off;
xlabel('eje x') %coamando xlabel para colocar nombre a el eje x
ylabel('eje y') %coamando label para colocar nombre a el eje y
title('Polinomio de Newton')%el comando title para colocar nombre a las graficas
end

function[F]=Diferencias_divididas_newton(x,y)
F=y';%primeras iteraciones del metodo
n=length(x);%n�mero de iteraciones del metodo
for i=2:n
    for j=2:i
        F(i,j)=(F(i,j-1)-F(i-1,j-1))/(x(i)-x(i-j+1));%diferencias divididas que dependen de los nodos x(i) y x(i-j+1)
    end
end
end

function[P]=Polinomio(F,xn)
syms x;
coeficientes=diag(F);%coeficientes del Polinomio
n=length(coeficientes);%grado del Polinomio
P=coeficientes(1);%termino independiente del polinomio

for i=2:n
    variable=kron(x,ones(1,i-1));% vector [x,x...,x] de dimension i-1
    variable=variable-xn(1:i-1);% vector [x-xn(1),x-xn(2),...,x-xn(i-1)]
    variable=prod(variable);% producto de los terminos del vector variable (x-xn(1))*(x-xn(2))*...(x-xn(i-1))
    P=P+coeficientes(i)*variable;%polinomio de las diferencias divididas
end

end

=======
function[]=Diferencias_divididas_interpolante_newton()
clc
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Entrada de Datos--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=[1.0 1.3 1.6 1.9 2.2];%nodos del eje x
yn=[0.7651977 0.6200860 0.4554022 0.2818186 0.1103623];%imagenes de los nodos
F=Diferencias_divididas_newton(xn,yn);%llamada al metodo interpolante de Newton
P=Polinomio(F,xn);%la llamada a esta funci�n es para construir el polinomio
P=matlabFunction(P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Parametros de salidad-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% F es una matrix digonal inferior cuyas entradas son los polinomios de las iteraciones
% F=diferencias_progresivas(x,y);
x0=1.5;%punto de prueba
y0=P(x0);%Imagen de punto de Prueba
close all
alfa=xn(length(xn))-xn(length(xn)-1);%extremo izquierdo del intervalo donde se quiere visualizar el estudio
alfa1=xn(length(xn))+xn(length(xn)-1);%extremo derecho del intervalo donde se quiere visualizar el estudio
hold on;fplot(P,[alfa,alfa1]);scatter(xn,yn);scatter(x0,y0);legend('Polinomio de Newton','nodos interpolantes','Nodo evaludo en el polinomio');hold off;
xlabel('eje x') %coamando xlabel para colocar nombre a el eje x
ylabel('eje y') %coamando label para colocar nombre a el eje y
title('Polinomio de Newton')%el comando title para colocar nombre a las graficas
end

function[F]=Diferencias_divididas_newton(x,y)
F=y';%primeras iteraciones del metodo
n=length(x);%n�mero de iteraciones del metodo
for i=2:n
    for j=2:i
        F(i,j)=(F(i,j-1)-F(i-1,j-1))/(x(i)-x(i-j+1));%diferencias divididas que dependen de los nodos x(i) y x(i-j+1)
    end
end
end

function[P]=Polinomio(F,xn)
syms x;
coeficientes=diag(F);%coeficientes del Polinomio
n=length(coeficientes);%grado del Polinomio
P=coeficientes(1);%termino independiente del polinomio

for i=2:n
    variable=kron(x,ones(1,i-1));% vector [x,x...,x] de dimension i-1
    variable=variable-xn(1:i-1);% vector [x-xn(1),x-xn(2),...,x-xn(i-1)]
    variable=prod(variable);% producto de los terminos del vector variable (x-xn(1))*(x-xn(2))*...(x-xn(i-1))
    P=P+coeficientes(i)*variable;%polinomio de las diferencias divididas
end

end

>>>>>>> f695588f6ac19f2dcfd3289f0e7eb5e7a1480266
