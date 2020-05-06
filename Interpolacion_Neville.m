function[]=Interpolacion_Neville()
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Entrada de Datos--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=[1.0 1.3 1.6 1.9 2.2];%nodos de el eje x
yn=[0.7651977 0.6200860 0.4554022 0.2818186 0.1103623];%nodos de el eje y
x0=1.5;%punto a evaluar en la interpolacion de Neville
Q=Neville(xn,yn);%la matriz Q es la matrix de la interpolaciones con los nodos
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Parametros de salidad-------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Q es una matrix digonal inferior cuyas entradas son los polinomios de las iteraciones
n=length(xn);%grado del polinomio
P=Q(n,n);%polinomio interpolante de Neville
P=matlabFunction(P);
Px0=P(x0);%imagen de x0 en el polinomio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Grafica de Resultados--------------------------------
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
alfa=xn(length(xn))-xn(length(xn)-1);%extremo izquierdo del intervalo donde se quiere visualizar el estudio
alfa1=xn(length(xn))+xn(length(xn)-1);%extremo derecho del intervalo donde se quiere visualizar el estudio
hold on;fplot(P,[alfa,alfa1]);scatter(xn,yn);scatter(x0,Px0);legend('Polinomio de Neville','nodos interpolantes','Nodo evaludo');hold off;
xlabel('eje x') %coamando xlabel para colocar nombre a el eje x
ylabel('eje y') %coamando label para colocar nombre a el eje y
title('Polinomio de Neville')%el comando title para colocar nombre a las graficas
save('Neville.mat','P') %guardando los resultados de Neville
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------- Metodo ----------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q]=Neville(xn,yn)
clc
syms x
format long
n=length(xn);%número de iteraciones de el algoritmo de Neville
%Q = sym('Q'); %este comando es pardefinir un arreglo de dimensiones n con parmetros x
Q=sym(zeros(n));
 Q(:,1)=sym(yn');%la primera columna de Q son las imagenes de los nodos x
  for i=2:n
      for j=2:i
         Q(i,j)=((x-xn(i-j+1))*Q(i,j-1)-(x-xn(i))*Q(i-1,j-1))/(xn(i)-xn(i-j+1));%construyendo la interpolacion con el nodos x(i) y x(i-j+1)
      end
  
  end

end

