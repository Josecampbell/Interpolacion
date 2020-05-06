function[]=Interpolacion_curvas()
syms x;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------Metodo---------------------------------------%
%---------------------Polinomio De Lagrange---------------------------%
%----------------------Para Curva--------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Entrada de Datos--------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xn=[0 pi/4 pi/2]; %conjunto de puntos a interpolar
f=sin(x)^2;
g=f;
f=matlabFunction(f); %funciï¿½n a interpolar 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Esquema iterativo------------------------------%
 %----------------------Para construir---------------------------------%
 %----------------------Polinomio de Lagrange--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=0;
    for i=1:length(xn) %esquema iterativo para de i-esimo termino de el polinomio de Lagrange
        L=Lagrange(xn,i);%polinomio de Lagrange para el nodo xn(i) 
        P=P+L*f(xn(i));%construyendo el polinomio
    end
F=simplify(P);    
P=matlabFunction(P); %polinomio de Lagrange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Representacion Grafica-------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
alfa=xn(1)-1;%%intervalo izquierdo donde se quiere visualizar el estudio
alfa1=xn(length(xn))+1;%intervalo derecho donde se quiere visualizar el estudio
fxn=P(xn);%imagenes de de los nodos a interpolar
hold on;fplot(P,[alfa,alfa1]);fplot(f,[alfa,alfa1],'MarkerEdgeColor','k');scatter(xn,fxn,'MarkerEdgeColor','r');legend('Polinomio de Lagrange','funcion original','Puntos Interpolante');hold off;
xlabel('eje x') %coamando xlabel para colocar nombre a el eje x
ylabel('eje y') %coamando label para colocar nombre a el eje y
title('Interpolación de curvas')%el comando title para colocar nombre a las graficas
%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------parametros de salidad--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% la InterpolaciÃ³Å„ retorna un polinomio P de grado length(xn)-1 generico
save('Interpolacion_curvas.mat','P') %comando save para guardar los resultados en un archivo .mat
%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Error estimado---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=Error(g,xn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------Informe de los resultados ----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID=fopen('Lagrange funcion.txt','w'); %nombre de el archivo txt llamado Lagrange funcion
Resultados(1,:)=xn;
Resultados(2,:)=f(xn);
Resultados(3,:)=P(xn);
Resultados(4,:)=abs(Resultados(3,:)-Resultados(2,:));
fprintf(fileID,'%5s\t \t \t %15s\t \t \t \t %10s\t \t \t \t %15s\t %15s\t %15s\n','xn','f(xn)','P(xn)','|P(xn)-f(xn)|','Error estimado','P(x)');
fprintf(fileID,'%20.15f\t %20.15f\t %20.15f\t %20.15f\t %20.15f\t %20s\n',Resultados(:,1),error,char(F));
fprintf(fileID,'%20.15f\t %20.15f\t %20.15f\t %20.15f\n',Resultados(:,2:size(Resultados,2)));
fclose(fileID);
end

function[L]=Lagrange(xn,i)
syms x;
n=length(xn)-1; %n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%----------------------Polinomio Lnk -----------------------------------%
%--------- formula obtenida del teroema 3.2 cap 3-----------------------%
%---------------------Libro Richard L. Burden --------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% indica el grado de los polinomios
numerador=kron(x,ones(1,n));% vector de variables [x,...,x] vector n entradas
%el comando kron es usando para construir vector de k elementos repetidos n
%veces
xk=xn(i);%el termino xk utilizado para calcular Lk
denominador=kron(xk,ones(1,n)); % vector de terminos xk [xk,...,xk]
xn(i)=[];%vector de [x1,x2,...,xk-1,xk+1,...,xn]
numerador=prod(numerador-xn);%producto de (x-xn(i)) para todo i=1...k-1,k+1,...n
%el comando prod es usado para calcular el producto
denominador=prod(denominador-xn);%producto de (x-xk) producto realizado n veces
L=numerador/denominador; %define explicitamente el polinomio de lagrange
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %----------------------Calculo del Error-------------------------------%
  %----------------------Para el metodo de-------------------------------%
  %----------------------Lagrange Teorema 3.3 Cap 3----------------------%
  %----------------------------Libro Richard L. Burden ------------------%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[error]=Error(f,xn)
syms x;
df=diff(f,length(xn)+1);%calculo de la  derivada 
alfa=max(xn);%
df=matlabFunction(df);
    polinomio=kron(x,ones(1,length(xn))); % vector de terminos xk [xk,...,xk] de tamaño n
    polinomio=polinomio-xn;%vector de(x-xn(i)) para todo i=1...k-1,k+1,...n
    polinomio=prod(polinomio);%producto de (x-xn(i)) para todo i=1...k-1,k,k+1,...n
    polinomio=simplify(polinomio);%simplificando el vector
    coeficientes=coeffs(polinomio);%calculando los coeficientes del vector
    raices=roots(coeficientes);%calculando las raices del vector
    maximo=double(max(raices));%calculando las raiz mas grande
    error=df(alfa)/2*maximo*factorial(length(xn)+1);% calculando el error
end