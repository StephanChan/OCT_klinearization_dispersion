function[yp]=get_interp(x,xp,y,indice_)
yp=zeros(size(y));
for ii = 1:length(x)
    xt = xp(ii);
    x0=x(indice_(ii,1));
    x1=x(indice_(ii,2));
    y0=y(indice_(ii,1));
    y1=y(indice_(ii,2));
    yp(ii)=y0+(xt-x0)*(y1-y0)/(x1-x0);
end