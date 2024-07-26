begin

; read PC and add filepaths here
; including u, v, tem and the itm of them

; set years
  yrSrt = 1980
  yrEnd = 2019
  year  = ispan(yrSrt, yrEnd, 1)
  nyrs  = dimsizes(year)

; read lat, lon, level, time
  f1 = addfile(filepath1, "r")
  f2 = addfile(filepath2, "r")
  f3 = addfile(filepath3, "r")

  lat = f2->lat
  lon = f2->lon
  level = f2->plev 
  time0 = tointeger(f2->time)
  time0@units = "seconds since 1970-01-01"
  time = cd_calendar(time0, -1)

;read u,v,tem climate
  fu = addfile(filepath4, "r")
  fv = addfile(filepath5, "r")
  ft = addfile(filepath6, "r")
  cu = fu ->u
  cv = fv ->v
  ct = ft ->t

  cuvar = dim_avg_n((/cu(2,:,:,:), cu(3,:,:,:), cu(4,:,:,:)/), 0) ; MAM
  cvvar = dim_avg_n((/cv(2,:,:,:), cv(3,:,:,:), cv(4,:,:,:)/), 0)
  ctvar = dim_avg_n((/ct(2,:,:,:), ct(3,:,:,:), ct(4,:,:,:)/), 0)
  copy_VarCoords(cu(0,:,:,:), cuvar)
  copy_VarCoords(cv(0,:,:,:), cvvar)
  copy_VarCoords(ct(0,:,:,:), ctvar)
 
  u_clim = cuvar 
  v_clim = cvvar
  t_clim = ctvar

; read u,v,tem
  t0  = f1 ->t
  u0  = f2 ->u
  v0  = f3 ->v
  u0&time = time 
  v0&time = time
  t0&time = time

  uAnom  = calcMonAnomTLLL(u0({yrSrt+"01":yrEnd+"12"},:,:,:), cu)
  vAnom  = calcMonAnomTLLL(v0({yrSrt+"01":yrEnd+"12"},:,:,:), cv)  
  tAnom  = calcMonAnomTLLL(t0({yrSrt+"01":yrEnd+"12"},:,:,:), ct)      

  uMar = uAnom({yrSrt+"03":yrEnd+"03":12},:,:,:)
  uApr = uAnom({yrSrt+"04":yrEnd+"04":12},:,:,:)
  uMay = uAnom({yrSrt+"05":yrEnd+"05":12},:,:,:)
  vMar = vAnom({yrSrt+"03":yrEnd+"03":12},:,:,:)
  vApr = vAnom({yrSrt+"04":yrEnd+"04":12},:,:,:)
  vMay = vAnom({yrSrt+"05":yrEnd+"05":12},:,:,:)
  tMar = tAnom({yrSrt+"03":yrEnd+"03":12},:,:,:)
  tApr = tAnom({yrSrt+"04":yrEnd+"04":12},:,:,:)
  tMay = tAnom({yrSrt+"05":yrEnd+"05":12},:,:,:)

  u = dim_avg_n((/uMar,uApr,uMay/), 0)
  v = dim_avg_n((/vMar,vApr,vMay/), 0) 
  t = dim_avg_n((/tMar,tApr,tMay/), 0)

;regression anomalies -> u', v', t'
  rc_u = regCoef_n(pc, u, 0, 0)
  rc_v = regCoef_n(pc, v, 0, 0)
  rc_t = regCoef_n(pc, t, 0, 0)
  copy_VarMeta(u_clim, rc_u)
  copy_VarMeta(u_clim, rc_v)
  copy_VarMeta(u_clim, rc_t)

;---------calculate CP------------
  R_cp = 0.286
  press = conform(t_clim, level, 0)
  copy_VarMeta(t_clim, press)

  dT_clim_dp = center_finite_diff_n(t_clim, level, False, 0, 0)

  sigma = (R_cp*t_clim/press) - dT_clim_dp
  copy_VarMeta(t_clim, sigma)

  pi = atan(1.0)*4.
  ff = 2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat) ;2Ωsin(lat)
  f  = conform(t_clim, ff, 1)

  du_clim_dp = center_finite_diff_n(u_clim, level, False, 0, 0)
  dv_clim_dp = center_finite_diff_n(v_clim, level, False, 0, 0)

  dcp = f/sigma*rc_v*rc_t*du_clim_dp + f/sigma*rc_u*rc_t*dv_clim_dp
  copy_VarMeta(t_clim, dcp)

  ; if vertical integral, need
  ;f7 = addfile(filepath7, "r")
  ;pres = f7->msl                                     ;read surface pressure (Pa)
  ;pMar = pres({yrSrt+"03":yrEnd+"03":12},:,:)
  ;pApr = pres({yrSrt+"04":yrEnd+"04":12},:,:)
  ;pMay = pres({yrSrt+"05":yrEnd+"05":12},:,:)
  ;pres_MAM = dim_avg_n_Wrap((pMar+pApr+pMay)/3., 0)  ;climatological MAM pres
  ;copy_VarCoords(pMar(0,:,:), pres_MAM)

  psfc      = pres_MAM
  ptop      = 100.      ; 1hPa
  dp        = dpres_plevel_Wrap(level, psfc, ptop, 0)
  dcp_new   = -dcp*dp
  cp = dim_sum_n_Wrap(dcp_new, 0)
  copy_VarMeta(dcp(0,:,:), cp)

;---------calculate Ck------------
  re = 6378388
  pi = atan(1.0)*4.

  dudx = new(dimsizes(u_clim), typeof(u_clim))
  dvdx = new(dimsizes(v_clim), typeof(v_clim))

  do i = 0, dimsizes(lat)-1
     dx = doubletofloat(lon*pi/180.*re*cos(lat(i)*pi/180))
     dudx(:,i:i,:) = center_finite_diff_n(u_clim(:,i:i,:), dx, True, 0, 2)
     dvdx(:,i:i,:) = center_finite_diff_n(v_clim(:,i:i,:), dx, True, 0, 2)
  end do

  dy = doubletofloat(lat*pi/180.*re)
  dudy = center_finite_diff_n(u_clim, dy, True, 0, 1)
  dvdy = center_finite_diff_n(v_clim, dy, True, 0, 1)

  ck = (rc_v*rc_v-rc_u*rc_u)/2.0*(dudx-dvdy)-rc_u*rc_v*(dudy+dvdx)
  copy_VarMeta(u_clim, ck)

;------------output---------------
  cp_file        = addfile("/Users/Desktop/cp.nc", "c")
  cp_file->cp    = dcp
  cp_file->cp    = cp

  ck_file        = addfile("/Users/Desktop/ck.nc", "c")
  ck_file->ck    = ck

end
