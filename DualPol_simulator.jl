using ArgParse
using NetCDF
using DataStructures

# Global constants
xam_r = pi*1000/6.
xbm_r = 3.
xmu_r = 0.
xobmr = 1./xbm_r
xcre = 1. + 2.*xbm_r + xmu_r

function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--output","-o"
    help = "path to output file"
    arg_type = String
    default = "dualpol_radar.nc"
    "input"
    help = "path to input file"
    arg_type = String
    required = true
  end
  return parse_args(s)
end

function read_WRF()
  ### Read in cartesian analysis
  centerX = 0; centerY = 0 ### TC center in km (in relation to x,y arrays)  # How was this defined??
  dims = ["XLAT","XLONG","ZNU","Times"]
  vars = ["QRAIN","QNRAIN","QVAPOR","PB","P","T","REFL_10cm"]
  lat,lon,eta,times = read_nc_var(filein,dims)
  qrain,qnrain,qv,pb,pp,tempk,refl = read_nc_var(filein,vars)
  d1,d2,d3,d4 = size(qrain)

  # Initialize new radar variables
  ZDR      = similar(qrain); fill!(ZDR, -999)
  ZV      = similar(qrain); fill!(ZV, -999)
  ZH      = similar(qrain); fill!(ZH, -999)
  DBZ      = similar(qrain); fill!(DBZ, -999)

  gamma1 = gamma(1. + xbm_r + xmu_r)
  gamma2 = gamma(1. + xmu_r)
  for t in [1:d4]
    for k in [1:d3]
      for j in [1:d2]
        for i in [1:d1]
          pressure = pb[i,j,k,t]+pp[i,j,k,t]
          rho = 0.622*pressure/(287.0*tempk[i,j,k,t]*(qv[i,j,k,t]+0.622))
          if (qrain[i,j,k,t] .gt. 1.0e-9) then
            qrv = qrain[i,j,k,t]*rho
            nrv = qnrain[i,j,k,t]*rho
          else
            qrv = 1.0e-12
            nrv = 1.0e-12
         end
         lambda = (xam_r*gamma1*xorg2*nrv/qrv)^xobmr
         N0 = nrv*xorg2*lambda^gamma2
         ZDR[i,j,k,t], ZV[i,j,k,t], ZA[i,j,k,t], DBZ[i,j,k,t] = calc_radar_variables(N0,lambda)
       end
     end
   end
 end
end


function calc_radar_variables()
  # Define radar constants
  eps = 8.88 + 0.63im
  wavelength = 100
  k = (abs2(eps) - 1)/(abs2(eps) +2)

  # Define the scattering amplitudes
  s_amp = Array(Complex64,16,2)
  s_amp = [ [5.93590e-05+(3.07940e-07)im] [-5.94150e-05+(-3.08530e-07)im]
            [4.70300e-04+(2.29610e-06)im] [-4.76460e-04+(-2.35930e-06)im]
            [1.55640e-03+(6.77020e-06)im] [-1.61730e-03+(-7.37320e-06)im]
            [3.58870e-03+(1.28470e-05)im] [-3.86500e-03+(-1.54380e-05)im]
            [6.77460e-03+(1.72310e-05)im] [-7.62290e-03+(-2.45720e-05)im]
            [1.12570e-02+(1.34310e-05)im] [-1.33110e-02+(-2.91580e-05)im]
            [1.71230e-02+(-9.24210e-06)im] [-2.13580e-02+(-1.72640e-05)im]
            [2.44090e-02+(-6.71130e-05)im] [-3.21760e-02+(3.40790e-05)im]
            [3.31110e-02+(-1.84310e-04)im] [-4.61260e-02+(1.67320e-04)im]
            [4.31850e-02+(-3.96140e-04)im] [-6.34650e-02+(4.59880e-04)im]
            [5.45400e-02+(-7.54470e-04)im] [-8.42680e-02+(1.05350e-03)im]
            [6.70130e-02+(-1.33600e-03)im] [-1.08320e-01+(2.21210e-03)im]
            [8.03320e-02+(-2.25460e-03)im] [-1.34920e-01+(4.44170e-03)im]
            [9.40660e-02+(-3.67760e-03)im] [-1.62650e-01+(8.75570e-03)im]
            [1.07580e-01+(-5.84290e-03)im] [-1.88770e-01+(1.73100e-02)im]
            [1.20010e-01+(-9.05960e-03)im] [-2.08160e-01+(3.50790e-02)im] ];

  # Calculate the reflectivity using 0.5 mm bins
  zh = 0.0
  zv = 0.0
  mu = 0
  N0 = 8000
  D0 = 1.5
  lambda = (3.67 + mu)/D0
  for i in [1:16]
      D = i*0.5
      N = N0*(D^mu)*exp(-lambda*D)
      zv += abs2(s_amp[i,1])*N
      zh += abs2(s_amp[i,2])*N
      println(D, "mm : ",N," #/m^3, H: ", zh, ", V:", zv)
  end
  zv *= (4*wavelength^4)/(pi^4*k^2)
  zh *= (4*wavelength^4)/(pi^4*k^2)
  za = N0*gamma(xcre)/(lambda^xcre)

  Zdr = 10*log10(zh/zv)
  Zv = 10*log10(zv)
  Zh = 10*log10(zh)
  Za = 10*log10(za)
  return Zdr, Zv, Zh, Za
end

function read_nc_var(filename::ASCIIString,varnames::Array)
  ###Read in Cartesian nc-file
  ###Only works as intended if you read in more than 1 variable
  #println("Read in nc file ...")
  data = OrderedDict()

  for varname in varnames
    data[string(varname)] = ncread (filename, string(varname))
  end

  return collect(values(data))
end

function write_ncfile_generic(filename,x,y,z,varnames)
  ###Write to nc-file
  ### data has to be in data["varname"]
  println("Write to nc file ...")

  ncvars = NcVar[]
  xatts = {"long_name" => "x (longitude)", "units" => "km", "missing_value" => -999, "_FillValue" => -999}
  yatts = {"long_name" => "y (latitude)",  "units" => "km", "missing_value" => -999, "_FillValue" => -999}
  zatts = {"long_name" => "z (altitude)",  "units" => "km", "missing_value" => -999, "_FillValue" => -999}
  tatts = {"long_name" => "time",  "units" => "s", "missing_value" => -999, "_FillValue" => -999}
  t = [1]
  x_dim = NcDim("longitude",lon,xatts)
  y_dim = NcDim("latitude",lat,yatts)
  z_dim = NcDim("eta",eta,zatts)
  t_dim = NcDim("t",times,tatts)

  for varname in varnames
    atts  = {"long_name" => varname, "units" => "???", "missing_value" => -999, "_FillValue" => -999}
    push!(ncvars,NcVar(varname,[x_dim,y_dim,z_dim,t_dim],atts,Float64))
  end

  nc = NetCDF.create(filename,ncvars)

  for varname in varnames
    NetCDF.putvar(nc,varname,data[varname])
  end

  NetCDF.close(nc)
  return 0
end

# Main program
args = parse_commandline()
filein = args["input"]
fileout = args["output"]
read_WRF()
vars_original = ["qrain,qnrain,qv,pb,pp,tempk,refl"]
vars_calculated = ["ZDR,ZV,ZH,DBZ"]
vars_out = [vars_original, vars_calculated]
write_ncfile(fileout,x,y,z,vars_out)
