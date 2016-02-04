
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"  

begin
in = addfile("data_obs.nc","r")
t = in->t
T = t(:1462,:,{-90.:90.},{0.:360.})
delete(t)
u = in->u
U = u(:1462,:,{-90.:90.},{0.:360.})
delete(u)

loni = in->longitude({0:360})
lati = in->latitude({-90.:90.})
T_file = addfile("../DART/hpss_files/Temp.nc","r")
newlat = T_file->lat(50:140)
newlon = T_file->lon

wgt_file = "ECMWF_to_DART.nc"

Opt = True
;Opt@InterpMethod = "patch"
;Opt@WgtFileName = wgt_file

;Opt@SrcGridLat = lati
;Opt@SrcGridLon = loni
;Opt@SrcInputFileName = in

;Opt@DstGridLat = newlat
Opt@DstGridLon = newlon

;Opt@ForceOverwrite = True
;Opt@Debug = True
;Opt@PrintTimings = True

obsT = ESMF_regrid_with_weights(T,wgt_file,1)
printVarSummary(obsT)

obsU = ESMF_regrid_with_weights(U,wgt_file,1)
printVarSummary(obsU)

ensemble_T = T_file->T
printVarSummary(ensemble_T)
ensemble_meanT = ensemble_T(time|:,copy|0,lev|(/13,18,28/),lat|50:140,lon|:)
ensemble_memberT = ensemble_T(copy|2:31,time|:,lev|(/13,18,28/),lat|50:140,lon|:)
delete(ensemble_T)

printVarSummary(ensemble_meanT)
printVarSummary(ensemble_memberT)

U_file = addfile("../DART/hpss_files/US.nc","r")
ensemble_U = U_file->US
printVarSummary(ensemble_U)

ensemble_meanU = ensemble_U(time|:,copy|0,lev|(/13,18,28/),slat|50:140,lon|:)
ensemble_memberU = ensemble_U(copy|2:31,time|:,lev|(/13,18,28/),slat|50:140,lon|:)
delete(ensemble_U)

printVarSummary(ensemble_meanU)
printVarSummary(ensemble_memberU)

n_ensemble = 30;
n_t = 1463
n_lon = 288
n_lat = 191
n_lev = 3


Depar = dim_sum_n((ensemble_meanT - obsT)^2,0)/(n_t-1);
printVarSummary(Depar)

Bias = dim_sum_n(ensemble_meanT-obsT,0)^2 /(n_t^2-n_t)
printVarSummary(Bias)

DeparU = dim_sum_n((ensemble_meanU-obsU)^2,0)/(n_t-1);
printVarSummary(DeparU)

BiasU = dim_sum_n(ensemble_meanT-obsT,0)^2/(n_t^2-n_t);
printVarSummary(BiasU)

EnsVar = dim_sum_n(dim_sum_n((ensemble_memberT - conform(ensemble_memberT,ensemble_meanT,(/1,2,3,4/)))^2,0),0)*(n_ensemble+1)/(n_ensemble*n_t*(n_ensemble-1))
printVarSummary(EnsVar)

EnsVarU = dim_sum_n(dim_sum_n((ensemble_memberU - conform(ensemble_memberU,ensemble_meanU,(/1,2,3,4/)))^2,0),0)*(n_ensemble+1)/(n_ensemble*n_t*(n_ensemble-1))
printVarSummary(EnsVarU)

Error_Resid = Depar - Bias - EnsVar
Error_ResidU = DeparU - BiasU - EnsVarU
printVarSummary(Error_ResidU)
 
plot = new(4,"graphic")
wks = gsn_open_wks("ps","Tropopause_Temp_")
gsn_define_colormap(wks,"amwg")
 res                       = True
  res@gsnMaximize           = True
  res@gsnDraw               = False   ; Don't draw plot
  res@gsnFrame              = False   ; Don't advance frame.

  res@gsnAddCyclic          = False
  res@cnFillOn              = True
  res@cnLinesOn             = False
  res@cnInfoLabelOn         = False
  res@gsnSpreadColors       = True
  res@gsnSpreadColorStart   = 12
  res@mpFillOn              = False 
  res@mpProjection          = "LambertEqualArea"
  res@mpGridAndLimbOn       = False

  res@mpDataBaseVersion     = "MediumRes"
  res@mpOutlineBoundarySets = "Geophysical"   ; "AllBoundaries"
 
 ; res@cnLevelSpacingF       = 5
  res@lbOrientation         = "Vertical"
 
  res@mpGridAndLimbDrawOrder = "Predraw"
  res@tiMainString           = "Lat/lon grid"

;---Plot will just be created, and not drawn yet.
  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks,Depar(0,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks,Bias(0,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks,EnsVar(0,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks,Error_Resid(0,:,:),res)

 resP                     = True                ; modify the panel plot
  resP@gsnMaximize         = True                ; make large
resP@txString = "Tropopause Temperature"
gsn_panel(wks,plot,(/4,1/),resP)

delete(plot)
plot = new(4,"graphic")
wks1 = gsn_open_wks("ps","Mid_Trop_Temp_")
gsn_define_colormap(wks1,"amwg")

  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks1,Depar(1,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks1,Bias(1,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks1,EnsVar(1,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks1,Error_Resid(1,:,:),res)
resP@txString = "Mid Troposphere Temperature"
gsn_panel(wks1,plot,(/4,1/),resP)

delete(plot)
plot = new(4,"graphic")
wks2 = gsn_open_wks("ps","Lower_Trop_Temp_")
gsn_define_colormap(wks2,"amwg")

  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks2,Depar(2,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks2,Bias(2,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks2,EnsVar(2,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks2,Error_Resid(2,:,:),res)
resP@txString = "Lower Troposphere Temperature"
gsn_panel(wks2,plot,(/4,1/),resP)

delete(plot)
plot = new(4,"graphic")
wks3 = gsn_open_wks("ps","Tropopause_U_")
gsn_define_colormap(wks3,"amwg")

  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks3,DeparU(0,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks3,BiasU(0,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks3,EnsVarU(O,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks3,Error_ResidU(0,:,:),res)
resP@txString = "Tropopause Zonal Wind"
gsn_panel(wks3,plot,(/4,1/),resP)

delete(plot)
plot = new(4,"graphic")
wks4 = gsn_open_wks("ps","Mid_Trop_U_")
gsn_define_colormap(wks4,"amwg")

  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks4,DeparU(1,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks4,BiasU(1,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks4,EnsVarU(1,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks4,Error_ResidU(1,:,:),res)
resP@txString = "Mid Troposphere Zonal Wind"
gsn_panel(wks4,plot,(/4,1/),resP)


delete(plot)
plot = new(4,"graphic")
wks5 = gsn_open_wks("ps","Lower_Trop_U_")
gsn_define_colormap(wks5,"amwg")

  res@gsnLeftString         = "Departure"
  plot(0) = gsn_csm_contour_map(wks5,DeparU(2,:,:),res)
  res@gsnLeftString         = "Bias"
  plot(1) = gsn_csm_contour_map(wks5,BiasU(2,:,:),res)
  res@gsnLeftString         = "Ensemble Variance"
  plot(2) = gsn_csm_contour_map(wks5,EnsVarU(2,:,:),res)
  res@gsnLeftString         = "Measurment Error and Residual"
  plot(3) = gsn_csm_contour_map(wks5,Error_ResidU(2,:,:),res)

resP@txString = "Lower Troposphere Zonal Wind"
gsn_panel(wks5,plot,(/4,1/),resP)

ncdf = addfile("Reliability.nc","c")
fAtt = True
fAtt@title = "Reliability Budget"
fAtt@source_file = "../DART/hpss_files/Temp.nc"
fAtt@Conventions = "None"
fAtt@creation_date = systemfunc("date")
fileattdef(ncdf,fAtt)

 filedimdef(ncdf,"time",-1,True) 

    ;===================================================================
    ; output variables directly; NCL will call appropriate functions
    ; to write the meta data associated with each variable
    ;===================================================================
       ncdf->T_Depar  = Depar                          ; 4D               
       ncdf->T_Bias = Bias
       ncdf->T_EnsVar = EnsVar
       ncdf->T_Error_Resid = Error_Resid
       ncdf->U_Depar = DeparU
       ncdf->U_Bias = BiasU
       ncdf->U_EnsVar = EnsVarU
       ncdf->U_Error_Resid = Error_ResidU

end
