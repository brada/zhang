
colours = 'rgbcmykw';
for i = 1:length(Calibration.Images)
    figure;
    imshow( Calibration.Images(i).Image );
    hold on
    iPt = Calibration.Images(i).iPt;
    wPt = Calibration.Images(i).wPt;
    for j = 1:length(iPt(1,:))
        code = mod( wPt(1,j) + wPt(2,j), 8 ) + 1;
        col = colours(code);
        plot( iPt(1,j), iPt(2,j), [col,'O'], 'MarkerSize',10, 'LineWidth',3 );
        s = sprintf( '%d,%d', wPt(1,j), wPt(2,j) );
        text( iPt(1,j)+10, iPt(2,j)+10, s );
    end
    hold off
end