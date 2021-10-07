function WN = compute_winding_number(c,h,M,rho,delta,epsilon,RAY_ANGLE,MESH,intplot)
% WN = compute_winding_number(c,h,M,rho,delta,epsilon,RAY_ANGLE,MESH)
% Computes winding number relative to the ray with angle RAY_ANGLE.
% Inputs:
% 1) c,h,M,rho: parameters from the model.
% 2) delta,epsilon: contour "fattening" parameters. Should be positive, small.
% 3) RAY_ANGLE: 1,2,3,4, for angles exp(i*pi*(RAY_ANGLE-1)/2)
% 4) MESH: mesh width for path-covering on the vertical segments. Lateral
% segments will have mesh scaled relative to ratio between vertical/lateral
% legnths to maintain consistent mesh throughout.
ih = intval(num2str(h)); iM = intval(num2str(M));
irho = intval(num2str(rho)); ic = c;
idelta = intval(num2str(delta)); iepsilon = intval(num2str(epsilon));
hat_omega = (irho*(abs(ih)+abs((1-ih)/iM)))^(1/intval(3));
southeast = idelta - 1i*iepsilon;
northeast = idelta + 1i*(iepsilon + hat_omega);
northwest = -idelta + 1i*(iepsilon + hat_omega);
southwest = -idelta - 1i*iepsilon;
corner = [southeast, northeast, northwest, southwest];
curve_tags = {'g(SouthEast -> NorthEast)' 'g(NorthEast -> NorthWest)' ...
    'g(NorthWest -> SouthWest)' 'g(SouthWest -> SouthEast)'};
g = @(z) z.^3 - ic*z.^2 + irho*(ih*(1-exp(z)) + (1-ih)/iM*(exp(-iM*z)-1));
dg = @(z) 3*z.^2 - 2*ic*z + irho*(ih*(-exp(z)) + (1-ih)*(-exp(-iM*z)));
crossings = [];
% Validate curves 1,....4
hold on
scale = mid(2*idelta/(hat_omega+2*iepsilon));
for m=1:4
    if m==1
        z = corner(m) + 1i*intval([(0:MESH:1),1])*(hat_omega+2*iepsilon);
    elseif m==2
        z = corner(m) - intval([(0:MESH/scale:1),1])*2*idelta;
    elseif m==3
        z = corner(m) - 1i*intval([(0:MESH:1),1])*(hat_omega+2*iepsilon);
    elseif m==4
        z = corner(m) + intval([(0:MESH/scale:1),1])*2*idelta;
    end
    dz = z(2)-z(1);
    u = midrad(mid(z),sup(abs(dz))*0.55).';
    gg = g(u).';
    if intplot==1
        plotintval(gg);
    else
        plot(mid(g(z)),'b','LineWidth',2);
    end
    cdir = [];
    % Check for intersections with zero
    fail = (inf(real(gg))<=0).*(sup(real(gg))>=0).*(inf(imag(gg))<=0).*(sup(imag(gg))>=0);
    if sum(fail)>0
        error([curve_tags{m},' contains zero. Decrease MESH.']);
    end
    % Determine all indices of g(u) that intersect the ray.
    if RAY_ANGLE==1     % positive real axis
        intersect = (sup(imag(gg))>=0).*(inf(imag(gg))<=0).*(inf(real(gg))>0);
    elseif RAY_ANGLE==2 % positive imaginary axis
        intersect = (sup(real(gg))>=0).*(inf(real(gg))<=0).*(inf(imag(gg))>0);
    elseif RAY_ANGLE==3 % negative real axis
        intersect = (sup(imag(gg))>=0).*(inf(imag(gg))<=0).*(sup(real(gg))<0);
    elseif RAY_ANGLE==4 % negative imaginary axis
        intersect = (sup(real(gg))>=0).*(inf(real(gg))<=0).*(sup(imag(gg))<0);
    end
    idx_intersect = intersect.*(1:length(z));
    if sum(idx_intersect)==0
        disp([curve_tags{m},' does not cross the ray.']);
    else
        v = gather_consecutive(idx_intersect.');
        vcol = size(v,2);
        for j=1:vcol
            vcol_j = v(:,j); vcol_j(isnan(vcol_j))=[];
            if RAY_ANGLE==1 || RAY_ANGLE==3
                dgg = imag(dg(u(vcol_j))*dz);
            elseif RAY_ANGLE==2 || RAY_ANGLE==4
                dgg = real(dg(u(vcol_j))*dz);
            end
            if min(dgg>0)==1
                if intplot==1
                    plotintval(gg(vcol_j),'y');
                else
                    plot(mid(gg(vcol_j)),'y','LineWidth',2);
                end
                if RAY_ANGLE==1 || RAY_ANGLE==4
                    crossings = [crossings;+1];
                elseif RAY_ANGLE==3 || RAY_ANGLE==2
                    crossings = [crossings;-1];
                end
                cdir = [cdir,crossings(end)];
            elseif min(dgg<0)==1
                if intplot==1
                    plotintval(gg(vcol_j),'y');
                else
                    plot(mid(gg(vcol_j)),'y','LineWidth',2);
                end
                if RAY_ANGLE==1 || RAY_ANGLE==4
                    crossings = [crossings;-1];
                elseif RAY_ANGLE==3 || RAY_ANGLE==2
                    crossings = [crossings;+1];
                end
                cdir = [cdir,crossings(end)];
            else
                clf
                error(['Cannot resolve a crossing of ',curve_tags{m},' with the ray. Decrease MESH or try a different RAY_ANGLE.']);
            end
        end
        disp([curve_tags{m}, ' crosses the ray ',num2str(vcol),' time(s) with direction(s) [',num2str(cdir),'].']);
    end
end
xlabel('$$Re(z)$$','FontSize',20,'Interpreter','LaTeX')
ylabel('$$Im(z)$$','FontSize',20,'Interpreter','LaTeX')
WN = sum(crossings);
disp(['Success! Winding number = ',num2str(WN)]);
axis square tight
plot(xlim,[0,0],'color',[0,0,0,0.35]);
plot([0,0],ylim,'color',[0,0,0,0.35]);
end