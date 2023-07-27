function pallete = getPallete(palleteSize)
    c = cell(1);
    c{1} = [1 1 0];
    c{2} = [1 0 0];
    c{3} = [0 1 1];
    c{4} = [0 0 1];
    C = reshape(vertcat(c{:}), 2, 2, []);
    pallete = zeros(palleteSize, palleteSize, 3);
    for i=1:3
        pallete(:,:,i) = interp2([0 1], [0 1], C(:,:,i), linspace(0,1,palleteSize), linspace(0,1,palleteSize).');
    end

    HSV = rgb2hsv(pallete);
    HSV(:, :, 2) = HSV(:, :, 2) * 10;
    HSV(HSV > 1) = 1;  % Limit values
    pallete = hsv2rgb(HSV);
end

