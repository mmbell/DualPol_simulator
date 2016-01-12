# DualPol_simulator.jl
#
# A dual-polarization radar simulator
#
# By Michael M. Bell
# Copyright 2015
#
# This program will simulate the backscatter from
# weather targets using a single, double, or fully
# explicit bin dropsize distribution. It is currently
# configured to work with WRF output, and was used in
# Brown, B. R., M. M. Bell, and A. J. Frambach, 2016:
# "Validation of Simulated Hurricane Drop Size Distributions
# using Polarimetric Radar", Geophys. Res. Lett., 42,
# doi:10.1002/2015GL067278.
#
# Inputs: Type of DSD (single, double, or full)
#         Name of WRF output file
# Outputs: NetCDF file with polarimetric radar variables
#
# Current output is radar reflectivity factor at horizontal
# and vertical polarization, and differential reflectivity
# for raindrops only. Additional variables and hydrometeor
# species will be added in due course
#

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
# Canting angle distribution width (sigma) equal 10 degrees only for rain & oblate crystals
sig = deg2rad(10.0)
r = exp(-2.0*sig*sig)
A = ( 0.375 + (0.5)*r + (0.125)*r^4 )
B = ( 0.375 - (0.5)*r + (0.125)*r^4 )
C = (0.125)*( 1 - r^4 )

# Parse the command line arguments
function parse_commandline()
  s = ArgParseSettings()

  @add_arg_table s begin
    "--dsdtype","-d"
    help = "single, double, full"
    arg_type = String
    required = true
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

# Calculate the radar variables using a prescribed gamma distribution
function calc_radar_gamma_dsd(N0,lambda)
  # Define radar constants

  # Fixed complex index of refraction
  eps = 8.88 + 0.63im

  # Radar wavelength (mm)
  wavelength = 100

  # Dielectric constant
  k = (abs2(eps) - 1)/(abs2(eps) +2)

  # Define the scattering amplitudes
  # Pre-calculated using T-matrix scattering code from Mishchenko (2000)
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
  h = 8.0/48.0 # Simpson's composite rule (8.0 - 0.0)/ (16 * 3)

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

# Calculate the radar variables using a full bin DSD
# This DSD follows the HUJI bin model output
function calc_radar_fulldsd(dsd, flag, qrv, N0,lambda)
  # Fixed complex index of refraction
  eps = 8.88 + 0.63im

  # Radar wavelength (mm)
  wavelength = 100

  # Dielectric constant
  k = (abs2(eps) - 1)/(abs2(eps) +2)

  # Define the scattering amplitudes
  s_amp = Array(Complex64,33,2)
  s_amp = [ [3.03150e-11+(1.59470e-13)im] [-3.04840e-11+(-1.61260e-13)im] #0.00400
          [6.06310e-11+(3.18960e-13)im] [-6.09670e-11+(-3.22510e-13)im] #0.00500
          [1.21270e-10+(6.37940e-13)im] [-1.21930e-10+(-6.44990e-13)im] #0.00640
          [2.42540e-10+(1.27600e-12)im] [-2.43860e-10+(-1.28990e-12)im] #0.00800
          [4.85100e-10+(2.55210e-12)im] [-4.87720e-10+(-2.57980e-12)im] #0.01000
          [9.70240e-10+(5.10470e-12)im] [-9.75410e-10+(-5.15930e-12)im] #0.01260
          [1.94060e-09+(1.02110e-11)im] [-1.95080e-09+(-1.03180e-11)im] #0.01600
          [3.88150e-09+(2.04240e-11)im] [-3.90140e-09+(-2.06340e-11)im] #0.02020
          [7.76370e-09+(4.08550e-11)im] [-7.80230e-09+(-4.12630e-11)im] #0.02540
          [1.55290e-08+(8.17280e-11)im] [-1.56040e-08+(-8.25140e-11)im] #0.03200
          [3.10630e-08+(1.63500e-10)im] [-3.12050e-08+(-1.65000e-10)im] #0.04040
          [6.21370e-08+(3.27090e-10)im] [-6.24050e-08+(-3.29910e-10)im] #0.05080
          [1.24300e-07+(6.54380e-10)im] [-1.24800e-07+(-6.59620e-10)im] #0.06400
          [2.48660e-07+(1.30920e-09)im] [-2.49560e-07+(-1.31870e-09)im] #0.08060
          [4.97470e-07+(2.61930e-09)im] [-4.99040e-07+(-2.63580e-09)im] #0.10160
          [9.95270e-07+(5.23990e-09)im] [-9.97890e-07+(-5.26760e-09)im] #0.12800
          [1.99120e-06+(1.04810e-08)im] [-1.99540e-06+(-1.05240e-08)im] #0.16120
          [3.98400e-06+(2.09560e-08)im] [-3.98980e-06+(-2.10170e-08)im] #0.20320
          [7.97050e-06+(4.18720e-08)im] [-7.97760e-06+(-4.19470e-08)im] #0.25600
          [1.59440e-05+(8.35630e-08)im] [-1.59520e-05+(-8.36450e-08)im] #0.32260
          [3.18850e-05+(1.66430e-07)im] [-3.18990e-05+(-1.66570e-07)im] #0.40640
          [6.37300e-05+(3.30320e-07)im] [-6.37970e-05+(-3.31020e-07)im] #0.51200
          [1.27260e-04+(6.51900e-07)im] [-1.27630e-04+(-6.55740e-07)im] #0.64500
          [2.53730e-04+(1.27450e-06)im] [-2.55450e-04+(-1.29240e-06)im] #0.81280
          [5.04610e-04+(2.45330e-06)im] [-5.11700e-04+(-2.52600e-06)im] #1.02400
          [9.99600e-04+(4.59930e-06)im] [-1.02620e-03+(-4.86780e-06)im] #1.29020
          [1.96830e-03+(8.22720e-06)im] [-2.06200e-03+(-9.14360e-06)im] #1.62540
          [3.84180e-03+(1.34180e-05)im] [-4.15370e-03+(-1.63240e-05)im] #2.04800
          [7.40440e-03+(1.73950e-05)im] [-8.39530e-03+(-2.58320e-05)im] #2.58040
          [1.40260e-02+(5.27710e-06)im] [-1.70300e-02+(-2.63510e-05)im] #3.25100
          [2.59700e-02+(-8.41070e-05)im] [-3.46010e-02+(5.15370e-05)im] #4.09600
          [4.67000e-02+(-4.92210e-04)im] [-6.97740e-02+(6.08330e-04)im] #5.16060
          [8.03860e-02+(-2.25910e-03)im] [-1.35030e-01+(4.45380e-03)im] #6.50200
          ];

  # Calculate the reflectivity using equal mass bins
  zh = 0.0
  zv = 0.0
  ql = 0.0
  NT = 0.0
  mu = 0.0
  h = 0.5
  r3 = 8.0e-9
  Dprevious = 0.0
  zv1 = 0.0
  zh1 = 0.0
  for i in [1:33]
    intcoeff = 1.0
    radius = r3^(1.0 / 3.0)
    mass = (4.0/3.0) * pi * r3 * 1.0e-6
    D = radius*2.0
    binbottom = (D + Dprevious)/2.0
    r3 = r3 * 2.0
    bintop = (D + 2.0 * r3^(1.0 / 3.0))/2.0
    intcoeff = 0.5 * (bintop - binbottom)
    q_bin = dsd[i]
    N = dsd[i]/(mass * (bintop - binbottom))
    Nexp = N0*(D^mu)*exp(-lambda*D)
    ql += q_bin
    NT += N*intcoeff
    if (flag)
      println("$D, $N, $Nexp, $q_bin")
    end
    # Jung et al (2010,JAMC) Eqns (4) and (3)
    zv2 = ( B*abs2(s_amp[i,2]) + A*abs2(s_amp[i,1]) + 2*C*real( conj(s_amp[i,1])*s_amp[i,2] ) ) * N
    zh2 = ( A*abs2(s_amp[i,2]) + B*abs2(s_amp[i,1]) + 2*C*real( conj(s_amp[i,1])*s_amp[i,2] ) ) * N
    zv += (zv1 + zv2)*intcoeff
    zh += (zh1 + zh2)*intcoeff
    Dprevious = D
    zv1 = zv2
    zh1 = zh2
   end
  #ql *= 1.0e-9
  #println("\n$NT, $nrv, $ql, $qrv, $rho")
  zv *= (4*wavelength^4)/(pi^4*k^2)
  zh *= (4*wavelength^4)/(pi^4*k^2)
  za = 10.0^(-3.5)

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

# Read in the WRF output
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

# Write out the radar variables
function write_ncfile(filename,lat,lon,lev,times,varnames)
  ###Write to nc-file
  println("Write to nc file ...")

  ncvars = NcVar[]
  xatts = {"long_name" => "x (longitude)", "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  yatts = {"long_name" => "y (latitude)",  "units" => "deg", "missing_value" => -999, "_FillValue" => -999}
  zatts = {"long_name" => "z (eta)",  "units" => "unitless", "missing_value" => -999, "_FillValue" => -999}
  tatts = {"long_name" => "time (minutes)",  "units" => "min", "missing_value" => -999, "_FillValue" => -999}
  x_dim = NcDim("east_west",[1:length(lon[:,1,1])],xatts)
  y_dim = NcDim("south_north",[1:length(lat[:,1,1])],yatts)
  z_dim = NcDim("bottom_top",[1:length(lev[:,1])],zatts)
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
  NetCDF.putvar(nc,"QNRAIN",qnrain)
  NetCDF.putvar(nc,"XLAT",lat)
  NetCDF.putvar(nc,"XLON",lon)
  NetCDF.putvar(nc,"XTIME",times)
  NetCDF.putvar(nc,"ZNU",lev)

  NetCDF.close(nc)
  return 0
end

# Main program
args = parse_commandline()
filein = args["input"]
fileout = args["output"]
dsdtype = args["dsdtype"]

### Read in Cartesian analysis
println("Reading ",filein)
fileinfo = ncinfo(filein)

dims = ["XLAT","XLONG","ZNU","XTIME"]
lat,lon,lev,times = read_nc_var(filein,dims)

if (dsdtype == "single")
  # Single moment
  vars = ["QRAIN","QVAPOR","PB","P","T","REFL_10CM"]
  qrain,qv,pb,pp,theta,refl = read_nc_var(filein,vars)
elseif (dsdtype == "double")
  # Double moment
  vars = ["QRAIN","QNRAIN","QVAPOR","PB","P","T","REFL_10CM"]
  qrain,qnrain,qv,pb,pp,theta,refl = read_nc_var(filein,vars)
elseif (dsdtype == "full")
  # Spectral bin
  vars = ["QRAIN","QNRAIN","QVAPOR","PB","P","T","REFL_10CM",
  "ff1i01","ff1i02","ff1i03","ff1i04","ff1i05","ff1i06","ff1i07","ff1i08",
  "ff1i09","ff1i10","ff1i11","ff1i12","ff1i13","ff1i14","ff1i15","ff1i16",
  "ff1i17","ff1i18","ff1i19","ff1i20","ff1i21","ff1i22","ff1i23","ff1i24",
  "ff1i25","ff1i26","ff1i27","ff1i28","ff1i29","ff1i30","ff1i31","ff1i32",
  "ff1i33"]
  dummyvar =  ncread (filein,"QRAIN")
  ff1 = Array(Array{Float32,4},33)
  for d in [1:33]
    ff1[d] = similar(dummyvar); fill!(ff1[d], -999)
  end
  ( qrain,qnrain,qv,pb,pp,theta,refl,
  ff1[1],ff1[2],ff1[3],ff1[4],ff1[5],ff1[6],ff1[7],ff1[8],
  ff1[9],ff1[10],ff1[11],ff1[12],ff1[13],ff1[14],ff1[15],ff1[16],
  ff1[17],ff1[18],ff1[19],ff1[20],ff1[21],ff1[22],ff1[23],ff1[24],
  ff1[25],ff1[26],ff1[27],ff1[28],ff1[29],ff1[30],ff1[31],ff1[32],ff1[33] ) = read_nc_var(filein,vars)
else
  println(dsdtype, " not recognized")
  exit
end

d1,d2,d3,d4 = size(qrain)

# Initialize new radar variables
ZDR      = similar(qrain); fill!(ZDR, -999)
ZV       = similar(qrain); fill!(ZV, -999)
ZH       = similar(qrain); fill!(ZH, -999)
DBZ      = similar(qrain); fill!(DBZ, -999)
Zdiff    = similar(qrain); fill!(Zdiff, -999)
ql       = similar(qrain); fill!(ql, -999)
qdiff    = similar(qrain); fill!(qdiff, -999)
if (dsdtype == "single")
  qnrain = similar(qrain); fill!(qnrain, -999)
elseif (dsdtype == "full")
  dsd = Array(Float32,33)
end

# Calculate the radar variables
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
          if (dsdtype == "single")
            N0 = 8000
            lambda = ((xam_r*gamma1*N0*1000.0/qrv)^(1./(1. + xbm_r)))/1000.0
            qnrain[i,j,k,t] = N0*gamma2/(lambda^(1. + xmu_r))
          else
            if (qnrain[i,j,k,t] > 0.0)
              nrv = qnrain[i,j,k,t]*rho
              lambda = ((xam_r*gamma1*xorg2*nrv/qrv)^xobmr)/1000.0
            else
              nrv = 0.0
              lambda = 1./20.E-3
            end
          end
          if (lambda < 1./2800.E-3)
            lambda = 1./2800.E-3
          elseif (lambda > 1./20.E-3)
            lambda = 1./20.E-3
          end
          if (dsdtype == "single")
            ZDR[i,j,k,t], ZV[i,j,k,t], ZH[i,j,k,t], DBZ[i,j,k,t], ql[i,j,k,t] = calc_radar_gamma_dsd(N0,lambda)
          elseif (dsdtype == "double")
            N0 = nrv*xorg2*lambda^(1. + xmu_r)
            ZDR[i,j,k,t], ZV[i,j,k,t], ZH[i,j,k,t], DBZ[i,j,k,t], ql[i,j,k,t] = calc_radar_gamma_dsd(N0,lambda)
          else
            if (i == 215 && j == 200 && k == 10 && t == 3)
              flag = true
              println("Found it!")
            else
              flag = false
            end
            for d in [1:33]
              dsd[d] = ff1[d][i,j,k,t]*rho
            end
            N0 = nrv*xorg2*lambda^(1. + xmu_r)
            ZDR[i,j,k,t], ZV[i,j,k,t], ZH[i,j,k,t], DBZ[i,j,k,t], ql[i,j,k,t] = calc_radar_fulldsd(dsd, flag, qrv, N0, lambda)
          end

          Zdiff[i,j,k,t] = refl[i,j,k,t] - ZH[i,j,k,t]
          ql[i,j,k,t] = 1.e-9*xam_r*N0*gamma1/(rho*lambda^(1. + xbm_r + xmu_r))
          if (DBZ[i,j,k,t] > 100.0)
            println("Large Z :",DBZ[i,j,k,t])
            println(lat[i,j,t],",",lon[i,j,t])
            println(i,",",j,",",k,":",nrv, " ", qrv, " ", rho, " ", refl[i,j,k,t])
            println(N0,",",lambda)
          elseif (DBZ[i,j,k,t] < -35.0)
            println("Small Z :",DBZ[i,j,k,t])
            println(lat[i,j,t],",",lon[i,j,t])
            println(i,",",j,",",k,":",nrv, " ", qrv, " ", rho, " ", refl[i,j,k,t])
            println(N0,",",lambda)
          end
        else
          qrv = 1.0e-12
          nrv = 1.0e-12
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

vars_out = ["REFL_10CM","ZDR","ZV","ZH","Rayleigh","Zdiff","QRAIN","QNRAIN"]
println("Writing ", fileout)
write_ncfile(fileout,lat,lon,lev,times,vars_out)
