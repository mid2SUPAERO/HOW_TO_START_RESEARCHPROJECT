function plotConfig(X,Yc,np,nY,nely)

factor=0.45;
el = zeros((nely*np),3); % Total number of printable layers
for i=1:length(Yc)-np
    if mod(i,np)~=0
        mult = floor(i/np);
    else
        mult = i/np-1;
    end
%     A = linspace(X(2*i-1),X(2*(i+1)-1),4);
%     B = linspace(y(i),y(i+1),4);
%     C = linspace(X(2*i),X(2*(i+1)),4);
    
    A = linspace(X(2*i-1),X(2*(i+1)-1),4);
    B = linspace(Yc(i),Yc(i+1),4);
    C = linspace(X(2*i),X(2*(i+1)),4);
    el(10*mult+i:np:10*mult+i+3*np,:) = [A' B' C'];
end
pos = [];
j=1;
for i=1:length(el)
    if el(i,3) < 0.5
        continue
    end
    pos(j,:) = [el(i,1)-el(i,3)/2 el(i,2) el(i,3) 1];
    j = j+1;
end
pos = [pos; 0 0 0 0];

for i=1:length(pos)
    rectangle('Position', pos(i,:), 'FaceColor', 'g', 'EdgeColor', 'r');
end
end
