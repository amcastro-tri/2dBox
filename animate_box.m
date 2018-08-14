close all;

skip = floor(0.02/h); 
scale = 0.5;
%wsz = [-0.6 -0.2 0 0.4];
wsz = [-1.0 0 0 1.2];
start_step = 1;
end_step = size(tt);
%skip=1;


a = lengths(1);
b = lengths(2);

X       = [-a a  a -a] / 2;
Y       = [-b -b b  b] / 2; 
 
 hSquare = fill(X,Y,'r');

 axis(wsz)
 axis equal;
 %h = gca;
 
 V0 = get(hSquare,'Vertices')';
 
 nsteps = length(xx);
 
 for istep=start_step:skip:end_step  
     x = xx(istep, 1);
     z = xx(istep, 2);
     theta = xx(istep, 3);
     
     % Rotation matrix:
     c = cos(theta);
     s = sin(theta);
     R_WB = [c, -s; s, c];
     
     C       = repmat([x z], 4, 1)';
     
     V = R_WB * V0 + C;             % do the rotation relative to the centre of the square
     
     fill(V(1,:),V(2,:),'r');
     %set(hSquare,'Vertices',V');    % update the vertices        
     hold on     
     
     % Forces
     %fn = YY(istep,10:12) * scale;
     %ft = zeros(size(fn)) * scale;     
     %xx = V(1, 1:3);
     %yy = V(2, 1:3);     
     quiver(V(1,:), V(2,:), ft(istep,:), fn(istep,:), 'AutoScale','off')   

      axis equal;
      axis(wsz)
      hold off;
     
     pause(0.1);
 end