%% Matheus, aí em baixo tem outras maneiras de calcular o transporte, que eu fiz pra comparar os métodos, mas utilizei pro meu tg o último

%%%%%%%% Calculo da integral da velocidade em todos os pontos e anos p 30S e LON FIXA(46)

%clear v_col int_vz soma t_dep dx  vsverd b%%

%v_col=NaN(24,176,1223);   
%int_vz=NaN(24,176,1223);

%t_dep=diff(dep(1:25));    % diferença da profundidade 
%dx= diff(lon)*111000*cosd(30);    %multiplicando p corrigir a dist (equador) e em km

%%%%%%%%%% Calculo da integral na coluna

%for j=1:1223 %loop do tempo

%for i=1:69  %loop dos pontos
%v_col(:,i,j)= (vv(1:24,i,j)+ vv(2:25,i,j))/2;
%int_vz(:,i,j)= v_col(:,i,j) .* t_dep’;
%end
%end

%soma= squeeze(nansum(int_vz,1)); % ficar com um valor de integral por ponto

%b= (soma(1:end-1,:)+soma(2:end,:))/2;   % para igualar numero de colunas com dx

%for q=1:1223   %ter um fluxo por tempo
%vsverd(:,q)= nansum(b(:,q).*dx)/1e6; %div para ter valor em sverdrup
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Calculo p encontrar lon que integral na coluna >0

%t_dep=diff(dep(1:24));    % diferença da profundidade 
%dx= diff(lon)*111000*(30);    %multiplicando p corrigir a dist (equador) e em km

%% Calculo da integral na coluna até 1500

%for j=1:1223 %loop do tempo

%for i=1:176  %loop dos pontos
%v_col(:,i,j)= (vv(1:23,i,j)+ vv(2:24,i,j))/2;
%int_vz(:,i,j)= v_col(:,i,j) .* t_dep’;
%end
%end

%soma= squeeze(nansum(int_vz,1)); % ficar com um valor de integral por ponto

%for i=1:1223       %%% definir lon pela integral
%a(i)=find((soma(:,i)>0),1);
%if a(i)>76
 %  a(i)=76;
%else
%  if a(i)<51
%     a(i)=76;
%end
%end 
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculo do transporte com longitude variável


%%% encontrar a longitude que v fica positivo
%%% coloquei 48.5 como limite inferior, pq tem velocidades bem prox à costa que são positivas
%%% coloquei 46 como máx, pq corrente do brasil não chega tão longe e há alguns valores que ultrapassam

for i=1:1223         %%% Definir lon pelas velocidades
g(i)=find((vv(1,:,i)>0),1);
if g(i)<45          %lon(45)=48.5, lon(69)= 46.5, lon(76)=46
g(i)=69;
else
   if g(i)>76
      g(i)=69;
end
end
end


clear v_col int_vz soma t_dep dx  vsverd b

v_col=NaN(24,176,1223);   
int_vz=NaN(24,176,1223);

t_dep=diff(dep(1:25));    % diferença da profundidade - 800m

dx= diff(lon(1:176)*111000*cosd(30));    %multiplicando p corrigir a dist (equador) e em km

for j=1:1223 %loop do tempo
for i=1:g(j)  %loop dos pontos ate lon(v>0)
%for i=1:69 %lon fixa

v_col(:,i,j)= (vv(1:24,i,j)+ vv(2:25,i,j))/2;  %%deltav
int_vz(:,i,j)= v_col(:,i,j) .* t_dep’;   %%deltaz


end
end


soma= squeeze(nansum(int_vz,1)); % ficar com um valor de integral por ponto

b= (soma(1:end-1,:)+soma(2:end,:))/2;   % para igualar numero de colunas com dx

for q=1:1223   %ter um fluxo por tempo
vsverd(:,q)= nansum(b(:,q).*dx)/1e6; %div para ter valor em sverdrup
end


%clear vsverd b soma int_vz v_col dx t_dep g — -8.77sv
