phi = (-180:.5:180).*pi/180;
theta = (-90:.5:90).*pi/180;

ux = sin(phi).'*cos(theta);
uy = cos(phi).'*cos(theta);
uz = repmat(sin(theta),[1 length(phi)]);

u(1,:) = ux(:);
u(2,:) = uy(:);
u(3,:) = uz(:);