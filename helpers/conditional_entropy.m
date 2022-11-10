function h = conditional_entropy(p)

function l = goodlog(s)
	if s == 0
		l = 0;
	else
		l = log2(s);
	end
end

d = length(p);
h = 0;
pB = sum(p,1);

for x = 1:d
	for y= 1:d
		h = h - p(x,y)*(goodlog(p(x,y)) - goodlog(pB(y)));
	end
end


end
