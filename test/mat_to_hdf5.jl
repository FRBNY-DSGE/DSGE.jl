# Convert a mat file to an hdf5 file. No compression is used.
# To use, run julia mat_to_hdf5.jl matFile1 matFile2 ...

using HDF5
#using MATLAB #uncomment to actually use this script; comment out to suppress warnings when not using

for matfile in ARGS
    
    # The absolute filepath and name of the .mat file
    fn_mat = normpath(joinpath(dirname(@__FILE__()),matfile))
    
    # The name of the file without the file extension
    fn_base = rsplit(fn_mat,".",2)[1]                     # for v0.4, change to rsplit(mat_name,".",limit=2)
    
    # HDF5 file name
    fn_h5 = "$fn_base.h5"
  
    # Open mat file and copy variables to HDF5 file
    if(ispath(fn_mat))

        println("Copying $matfile to hdf5...")
        
        mf = MatFile(fn_mat, "r")
        vars = variable_names(mf)

        h5 = h5open(fn_h5, "w") 

        for varname in vars
            var = get_variable(mf, varname)
              try
                write(h5, varname, var)
            catch
                # This can happen when the variable is of type "Any"  
                println("Warning: could not write variable $varname in $fn_mat; you probably want to investigate")
            end
        end
        
        close(h5)
        
        println("Wrote $fn_h5")
        MATLAB.close(mf)

    else
    
    println("$matfile not found")
        
    end
end

println("Done.")
