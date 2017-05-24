A = randn(10, 20, 30, 40);

B = util.ndslice(A, 3, 15);

A == B
