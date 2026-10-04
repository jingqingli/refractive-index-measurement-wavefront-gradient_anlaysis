function [p2,upwrap_phase]=phase_fit(phase_matrix,delta)
[nx,ny]=size(phase_matrix);
valid_idx = ~isnan(phase_matrix);
% delta=2.4e-6;
% delta=150e-6;
% x0=1:nx;y0=1:ny;x0=x0*delta;y0=y0*delta;
x0=-nx/2:nx/2-1;y0=-ny/2:ny/2-1;x0=x0*delta;y0=y0*delta;
[X0,Y0]=meshgrid(y0,x0);%物面二维网格
% X0 = X0 / max(abs(X0(:)));
% Y0 = Y0 / max(abs(Y0(:)));
x_valid = X0(valid_idx);
y_valid = Y0(valid_idx);
z_valid = phase_matrix(valid_idx);
[p2,gof2]=fit([x_valid(:),y_valid(:)],z_valid(:),'poly22','Robust', 'on'); %进行二次
% if output.iterations >= 50
%     disp('达到稳健迭代上限，程序终止。');
%     return;
% end
% [p2,gof2]=fit([x_valid(:),y_valid(:)],z_valid(:),'pchipinterp','Robust', 'on') %进行二次
% unwrap_phase= fit([x_valid(:),y_valid(:)], z_valid(:), 'lowess', 'Span', 0.2);
upwrap_phase=p2.p10*X0+p2.p01*Y0+p2.p00+p2.p20*X0.^2+p2.p02*Y0.^2+p2.p11*X0.*Y0;
end
