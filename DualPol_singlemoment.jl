using ArgParse
using NetCDF
using DataStructures
using Debug

# Global constants
xam_r = pi*997.0/6.0
xbm_r = 3.0
xmu_r = 0.0
xobmr = 1.0/xbm_r
xcre = 1. + 2.*xbm_r + xmu_r

# Angular moments from Jung et al (2010,JAMC) Eqns (4)
# canting angle distribution width (sigma) equal 10 degrees only for rain & oblate crystals
sig = deg2rad(10.0)
r = exp(-2.0*sig*sig)
A = ( 0.375 + (0.5)*r + (0.125)*r^4 )
B = ( 0.375 - (0.5)*r + (0.125)*r^4 )
C = (0.125)*( 1 - r^4 )

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

function calc_radar_variables(N0,lambda)
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
            [1.20010e-01+(-9.05960e-03)im] [-2.08160e-01+(3.50790e-02)im]
            [0.0+(0.0)im] [0.0+(0.0)im] ];

  # Calculate the reflectivity using 0.5 mm bins
  zh = 0.0
  zv = 0.0
  ql = 0.0
  mu = 0.0
  h = 8.0/48.0
  for i in [1:17]
    if (i == 1) || (i == 17)
      intcoeff = 1.0
    elseif (i%2 == 0)
      intcoeff = 4.0
    else
      intcoeff = 2.0
    end
      D = i*0.5
      N = N0*(D^mu)*exp(-lambda*D)
      ql += xam_r*N0*(D^(3.0+mu))*exp(-lambda*D)*intcoeff
      # Jung et al (2010,JAMC) Eqns (4) and (3)
      zv += ( B*abs2(s_amp[i,2]) + A*abs2(s_amp[i,1]) + 2*C*real( conj(s_amp[i,1])*s_amp[i,2] ) )*N*intcoeff
      zh += ( A*abs2(s_amp[i,2]) + B*abs2(s_amp[i,1]) + 2*C*real( conj(s_amp[i,1])*s_amp[i,2] ) )*N*intcoeff
  end
  ql *= h*1.0e-9
  zv *= h*(4*wavelength^4)/(pi^4*k^2)
  zh *= h*(4*wavelength^4)/(pi^4*k^2)
  za = N0*gamma(xcre)/(lambda^xcre)

  if (zh > 10.0^(-3.5)) && (zv > 10.0^(-3.5))
   Zdr = 10*log10(zh/zv)
   Zv = 10*log10(zv)
   Zh = 10*log10(zh)
  else
   Zdr = 0.0
   Zv = -35.0
   Zh = -35.0
  end

  if (za > 10.0^(-3.5))
   Za = 10*log10(za)
  else
   Za = -35.0
  end

  return Zdr, Zv, Zh, Za, ql
end

function read_nc_var(filename,varnames::Array)
  ###Read in Cartesian nc-file
  ###Only works as intended if you read in more than 1 variable
  #println("Read in nc file ...")
  data = OrderedDict()

  for varname in varnames
    data[string(varname)] = ncread (filename, string(varname))
  end

  return collect(values(data))
end

@debug function write_ncfile(filename,lat,lon,eta,times,varnames)
  ###Write to nc-file
  println("Write to nc file ...")

  ncvars = NcVar[]
  xatts = {"long_name" => "x (longitude)", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  yatts = {"long_name" => "y (latitude)",  "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  zatts = {"long_name" => "z (eta)",  "units" => "unitless", "missing_value" => -999, "_FillValue" => -999}
  tatts = {"long_name" => "time (minutes)",  "units" => "min", "missing_value" => -999, "_FillValue" => -999}
  x_dim = NcDim("east_west",[1:length(lon[:,1,1])],xatts)
  y_dim = NcDim("south_north",[1:length(lat[:,1,1])],yatts)
  z_dim = NcDim("bottom_top",[1:length(eta[:,1])],zatts)
  t_dim = NcDim("Time",[1:length(times)],tatts)

  for varname in varnames
    atts  = {"long_name" => varname, "units" => "???", "missing_value" => -999, "_FillValue" => -999}
    push!(ncvars,NcVar(varname,[x_dim,y_dim,z_dim,t_dim],atts,Float64))
  end
  atts  = {"long_name" => "Latitude", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  push!(ncvars,NcVar("XLAT",[x_dim,y_dim,t_dim],atts,Float64))
  atts  = {"long_name" => "Longitude", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  push!(ncvars,NcVar("XLON",[x_dim,y_dim,t_dim],atts,Float64))
  atts  = {"long_name" => "Time", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  push!(ncvars,NcVar("XTIME",[t_dim],atts,Float64))
  atts  = {"long_name" => "Eta Levels", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  push!(ncvars,NcVar("ZNU",[z_dim,t_dim],atts,Float64))
  nc = NetCDF.create(filename,ncvars)

  NetCDF.putvar(nc,"REFL_10CM",refl)
  NetCDF.putvar(nc,"ZDR",ZDR)
  NetCDF.putvar(nc,"ZV",ZV)
  NetCDF.putvar(nc,"ZH",ZH)
  NetCDF.putvar(nc,"Rayleigh",DBZ)
  NetCDF.putvar(nc,"Zdiff",Zdiff)
  NetCDF.putvar(nc,"QRAIN",qrain)
  NetCDF.putvar(nc,"XLAT",lat)
  NetCDF.putvar(nc,"XLON",lon)
  NetCDF.putvar(nc,"XTIME",times)
  NetCDF.putvar(nc,"ZNU",eta)

  NetCDF.close(nc)
  return 0
end

# Main program
args = parse_commandline()
filein = args["input"]
fileout = args["output"]
#read_WRF()
 ### Read in cartesian analysis
  println("Reading ",filein)
  dims = ["XLAT","XLONG","ZNU","XTIME"]
  vars = ["QRAIN","QVAPOR","PB","P","T","REFL_10CM"]
  lat,lon,eta,times = read_nc_var(filein,dims)
  qrain,qv,pb,pp,theta,refl = read_nc_var(filein,vars)
  d1,d2,d3,d4 = size(qrain)

  # Initialize new radar variables
  ZDR      = similar(qrain); fill!(ZDR, -999)
  ZV       = similar(qrain); fill!(ZV, -999)
  ZH       = similar(qrain); fill!(ZH, -999)
  DBZ      = similar(qrain); fill!(DBZ, -999)
  Zdiff    = similar(qrain); fill!(Zdiff, -999)
  ql       = similar(qrain); fill!(ql, -999)
  qdiff    = similar(qrain); fill!(qdiff, -999)

  println("Calculating radar variables...")
  gamma1 = gamma(1. + xbm_r + xmu_r)
  gamma2 = gamma(1. + xmu_r)
  xorg2 = 1.0/gamma2
  for t in [1:d4]
    for k in [1:d3]
      for j in [1:d2]
        for i in [1:d1]
          qvapor = max(1.0e-10,qv[i,j,k,t])
          pressure = pb[i,j,k,t]+pp[i,j,k,t]
          tempk = (theta[i,j,k,t]+300.0)*(pressure/100000.0)^(2.0/7.0)
          rho = 0.622*pressure/(287.15*tempk*(qvapor+0.622))
          if (qrain[i,j,k,t] > 1.0e-9)
            qrv = qrain[i,j,k,t]*rho
	          N0 = 8000
            lambda = ((xam_r*gamma1*N0*1000.0/qrv)^(1./(1. + xbm_r)))/1000.0
            if (lambda < 1./2800.E-3)
              lambda = 1./2800.E-3
            elseif (lambda > 1./20.E-3)
              lambda = 1./20.E-3
            end
            ZDR[i,j,k,t], ZV[i,j,k,t], ZH[i,j,k,t], DBZ[i,j,k,t], ql[i,j,k,t] = calc_radar_variables(N0,lambda)
            Zdiff[i,j,k,t] = DBZ[i,j,k,t] - ZV[i,j,k,t]
            ql[i,j,k,t] = 1.e-9*xam_r*N0*gamma1/(rho*lambda^(1. + xbm_r + xmu_r))
            if (DBZ[i,j,k,t] > 100.0)
              println("Large Z :",DBZ[i,j,k,t])
              println(lat[i,j,t],",",lon[i,j,t])
              println(i,",",j,",",k,":", qrv, " ", rho, " ", refl[i,j,k,t])
              println(N0,",",lambda)
            elseif (DBZ[i,j,k,t] < -35.0)
              println("Small Z :",DBZ[i,j,k,t])
              println(lat[i,j,t],",",lon[i,j,t])
              println(i,",",j,",",k,":", qrv, " ", rho, " ", refl[i,j,k,t])
              println(N0,",",lambda)
            end
          else
            qrv = 1.0e-12
            ql[i,j,k,t] = 1.0e-9
            qrain[i,j,k,t] = 1.0e-9
            ZDR[i,j,k,t] = 0.0
            ZV[i,j,k,t] = ZH[i,j,k,t] = DBZ[i,j,k,t] = -35.0
         end
         qdiff[i,j,k,t] = ql[i,j,k,t] / qrain[i,j,k,t]
       end
     end
   end
 end

vars_out = ["REFL_10CM","ZDR","ZV","ZH","Rayleigh","Zdiff","QRAIN"]
println("Writing ", fileout)
write_ncfile(fileout,lat,lon,eta,times,vars_out)
