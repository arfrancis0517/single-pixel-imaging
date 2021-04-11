function y = scaleimg(x, newdims)
	newh = newdims(1);
	neww = newdims(2);
	
	[h, w] = size(x);
	
	nx = zeros(h+1,w+1);
	nx(1:h, 1:w) = x;
	x = nx;
	
	dh = h / newh;
	dw = w / neww;
	area = dh * dw;
	
	party = zeros(newh, w+1);
	for i = 1:newh
		hh = dh * i;
		lh = dh * (i-1);
		lby = floor(lh);
		hby = floor(hh);
		
		party(i, :) = x(lby+1, :) .* (lby + 1 - lh);
		party(i, :) = party(i,:) + x(1+hby, :) .* (hh - hby);
		for li = lby+1 : hby -1
			party(i,:) = party(i,:) + x(li+1, :);
		end
	end
	
	partx = zeros(newh, neww);
	for i = 1:neww
		hw = dw * i;
		lw = dw * (i-1);
		lbx = floor(lw);
		hbx = floor(hw);
		
		partx(:, i) = party(:, lbx+1) .* (lbx + 1 - lw);
		partx(:, i) = partx(:,i) + party(:, hbx+1) .* (hw - hbx);
		for li = lbx+1 : hbx -1
			partx(:, i) = partx(:,i) + party(:, li+1);
		end
	end
	
	y = partx ./ area;
end

