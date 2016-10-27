function p = exact_proj(z)
	d = size(z);
	sup_i = z > 1;
	inf_i = z < 0;
	p = z;
	p(sup_i) = 1;
	p(inf_i) = 0;
