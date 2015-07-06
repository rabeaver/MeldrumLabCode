function imagehsv(x,y,data)
	datasize = size(data);
	data = reshape(data,1,[]);
	maxdata = max(abs(data));
	colorindex = fix(round(mod(angle(data),2*pi)*255/2/pi))+1;
	amplitude = abs(data)/maxdata;
	hsvmap = hsv(256);
    alldata = hsvmap(colorindex,:).*repmat(amplitude',[1 3]);
	image(x,y,reshape(alldata,[datasize 3]));
end
