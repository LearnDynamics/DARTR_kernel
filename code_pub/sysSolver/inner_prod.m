function v = inner_prod(f1, f2, range)
dx = 0.00001; xmin = range(1); xmax = range(end);
xmesh = xmin:dx:xmax;

v =  sum(f1(xmesh).* f2(xmesh)*dx);

end