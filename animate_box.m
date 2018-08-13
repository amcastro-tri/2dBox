close all;

skip = 5; 
scale = 0.5;
wsz = [-1 0 0 1.2];

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
 
 for istep=1:skip:nsteps  
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
%      fn = YY(istep,10:12) * scale;
%      ft = zeros(size(fn)) * scale;     
%      xx = V(1, 1:3);
%      yy = V(2, 1:3);     
%      quiver(xx, yy, ft, fn, 'AutoScale','off')   

      axis equal;
      axis(wsz)
      hold off;
     
     pause(0.02);
 end