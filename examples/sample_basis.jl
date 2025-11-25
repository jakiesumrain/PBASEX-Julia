include("../src/pbasex_core.jl")
using .pBASEX
#using CSV, DataFrames, HDF5

# Define the basis function parameters
l_values = collect(0:2:4)  # Legendre polynomial orders
x_values = collect(0:511)  # radial grid points

#@time basis_file = BASISData{Float64}("./", x = x_values, l = l_values)
basis_file = BASISData{Float64}("./basis.h5", x = x_values, l = l_values)


#println(minimum(basis_file.S))
#println(maximum(basis_file.S))

# add a histogram of the singular values
#histogram(basis_file.S)
#savefig("singular_values_histogram_1.png")

# Save singular values to text file
#open("singular_values_1.txt", "w") do io
    #writedlm(io, basis_file.S)
#end

# add a plot of the singular values
#plot(basis_file.S)
#savefig("singular_values_plot_1.png")


