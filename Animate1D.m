function [MovAnim,zlims] = Animate1D(state,spgrid,tgrid,Tpause,record,zlims)
% function [MovAnim,zlims] = Animate1D(state,spgrid,tgrid,BCtype,Tpause,record,zlims)
%
% Plot the solution a 1D partial differential equation (heat, wave, or beam)
% state = state of the closed-loop system at times 'tgrid'
% spgrid = spatial grid
% tgrid = grid for time
% Tpause = pause time in the animation
% record = 1 if animation is recorded into a movie
% zlims = limits for the z-axis (optional)
%
% Copyright (C) 2020 by Lassi Paunonen (lassi.paunonen@tuni.fi)
% Licensed under GNU GPLv3 (see LICENSE.txt).

if max(max(abs(imag(state)))) > 1e-8
  warning('Solution may contain imaginary parts that are ignored in the plot.')
end

state = real(state);

if nargin <= 6
  zlims(1) = min(min(min(state)));
  zlims(2) = max(max(max(state)));
end
axlims = [spgrid(1) spgrid(end) zlims];
    


if record == 1
  
  MovAnim = struct('cdata',[],'colormap',[]);
  
  for ind = 1:size(state,2)
    
    plot(spgrid,state(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    MovAnim(ind) = getframe(gcf);
    pause(Tpause)
    
  end
  
else
  
  MovAnim = [];
  for ind = 1:size(state,2)
    
    plot(spgrid,state(:,ind).','Linewidth',2)
    axis(axlims)
    xlabel('$\xi$','Interpreter','latex','Fontsize',20)
    title(['Time $=\; ' num2str(tgrid(ind),'%.1f') '$'],'Interpreter','latex','Fontsize',20)
    set(gcf,'color',1/255*[252 247 255])
    drawnow
    pause(Tpause)
    
  end
end
