function [] = gifanim(fname, framerate, fighandl, loopind)

% Capture frame
frame = getframe(fighandl);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);

if loopind == 1

    imwrite(imind,cm,[fname '.gif'],'gif', ...
        'LoopCount',inf, ...
        'DelayTime',1/framerate);

else

    imwrite(imind,cm,[fname '.gif'],'gif', ...
        'WriteMode','append', ...
        'DelayTime',1/framerate);

end

end