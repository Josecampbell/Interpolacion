function[]=Interpolacion_puntos()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %-----------------------          ------------------------------------%
 %--------Metodo de Polinomio de Interpolante--------------------------%
 %--------------De Lagrange Para Puntos--------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Entrada de Datos--------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=2007:2017;
xn(10)=[];
yn=[5.3 1.5 -3.2 -1.5 4.2 5.6 4.2 -3.9 -5.7 -12];
x0=2016;%puntos de prueba
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Esquema iterativo-------------------------------%
 %---------------Para construir el polinomio----------------------------%
 %----------------------De Lagrange-------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=0;
    for i=1:length(xn)%esquema iterativo para de i-esimo termino de el polinomio de Lagrange
        L=Polinomio(xn,i);%polinomio de Lagrange para el nodo xn(i) 
        P=P+L*yn(i);%construyendo el polinomio a diferencia de P_Lagrange este no requiere una funciÃ³n f solo las imagenes de los xn
    end 
F=simplify(P);
P=matlabFunction(P); %polinomio de Lagrange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------parametros de salidad--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la Interpolación„ retorna un polinomio P de grado length(xn)-1 generico
% para todo punto x0
save('Interpolacion_puntos.mat','P') %comando save para guardar los resultados en un archivo .mat
Px0=P(x0);%imagen de x0 en el polinomio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Grafica de Resultados--------------------------%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
%alfa=min(yn);%extremo izquierdo del intervalo donde se quiere visualizar el estudio
%alfa1=max(yn);%extremo derecho del intervalo donde se quiere visualizar el estudio
hold on;fplot(P,[xn(1),xn(length(xn))]);scatter(xn,yn,'MarkerEdgeColor','r');scatter(x0,Px0,'MarkerEdgeColor','g');legend('Polinomio de Lagrange','nodos interpolantes','punto evaluado');hold off;
xlabel('Año') %coamando xlabel para colocar nombre a el eje x
ylabel('Producto Interno Bruto') %coamando label para colocar nombre a el eje y
title('Interpolación Puntos')%el comando title para colocar nombre a las graficas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------Informe de los resultados ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID=fopen('Lagrange puntos.txt','w'); %nombre de el archivo txt llamado Lagrange puntos
Resultados(1,:)=xn;
Resultados(2,:)=yn;
Resultados(3,:)=P(xn);
Resultados(4,:)=abs(P(xn)-yn);
fprintf(fileID,'%5s\t \t \t %15s\t \t %10s\t \t %10s\t','%10s\n','xn','yn','P(xn)','|p(xn)-yn|','P(x)');
fprintf(fileID,'%20.15f\t %20.15f\t %20.15f\t %20.15f\t %20s\n',Resultados(:,1),char(F));
fprintf(fileID,'%20.15f\t %20.15f\t %20.15f\t %20.15f\n',Resultados(:,2:size(Resultados,2)));
fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Polinomio Lnk -----------------------------------%
%--------- formula obtenida del teroema 3.2 cap 3-----------------------%
%---------------------Libro Richard L. Burden --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[L]=Polinomio(xn,i)
syms x;
n=length(xn)-1; %n indica el grado de los polinomios
numerador=kron(x,ones(1,n));% vector de variables [x,...,x] vector n entradas
%el comando kron es usando para construir vector de k elementos repetidos nveces
xk=xn(i);%el termino xk utilizado para calcular Lk
denominador=kron(xk,ones(1,n)); % vector de terminos xk [xk,...,xk]
xn(i)=[];%vector de [x1,x2,...,xk-1,xk+1,...,xn]
numerador=prod(numerador-xn);%producto de (x-xn(i)) para todo i=1...k-1,k+1,...n
%el comando prod es usado para calcular el producto
denominador=prod(denominador-xn);%producto de (x-xk) producto realizado n veces
L=numerador/denominador; %define explicitamente el polinomio de lagrange
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Polinomio de Lagrange ---------------------------%
%--------- formula obtenida del teroema 3.5 cap 3-----------------------%
%---------------------Libro Richard L. Burden --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%