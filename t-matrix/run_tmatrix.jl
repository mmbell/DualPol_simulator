diam = 0.5
lam = 100.0
mrr = 8.88
mri = 0.63
phi = 180.0
s11 = Ref{ComplexF64}(0.0 + 0.0im)
s12 = Ref{ComplexF64}(0.0 + 0.0im)
s21 = Ref{ComplexF64}(0.0 + 0.0im)
s22 = Ref{ComplexF64}(0.0 + 0.0im)
ccall((:tmatrix_, "./libtmatrix.so"), Cvoid, (Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}, Ref{ComplexF64}), diam, lam, mrr, mri, phi, s11, s12, s21, s22)
println("s11 = $(s11)")
println("s11 = $(s12)")
println("s21 = $(s21)")
println("s22 = $(s22)")

