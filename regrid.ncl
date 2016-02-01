
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"  

begin
in = addfile("data_obs.nc","r")

T = in->t(:,:,{-90.:90.},{0.:360.})
U = in->u(:,:,{-90.:90.},{0.:360.})

lon = in->longitude({0:360})
lati = in->latitude({-90.:90.})
T_file = addfile("../DART/hpss_files/Temp.nc","r")
newlat = T_file->lat
newlon = T_file->lon

obsT = linint2_Wrap(lon,lati,T,True,newlon,newlat,0)
obsU = linint2_Wrap(lon,lati,U,True,newlon,newlat,0)


ensemble_meanT = T_file->T(:,0,:,:,(/13,18,28/))
ensemble_memberT = T_file->T(:,2:,:,:,(/13,18,28/))
U_file = addfile("../DART/hpss_files/US.nc","r")
ensemble_meanU = U_file->US(:,0,:,:,(/13,18,28/))
ensemble_memberU = U_file->US(:,2:,:,:,(/13,18,28/))
n_ensemble = 30;
n_t = 1463
n_lon = 288
n_lat = 191
n_lev = 3

Depar =new((/n_lev,n_lat,n_lon/),"float",0);
Bias =new((/n_lev,n_lat,n_lon/),"float",0);
EnsVar =new((/n_lev,n_lat,n_lon/),"float",0);

DeparU =new((/n_lev,n_lat,n_lon/),"float",0);
BiasU =new((/n_lev,n_lat,n_lon/),"float",0);
EnsVarU =new((/n_lev,n_lat,n_lon/),"float",0);
Error_ResidU =new((/n_lev,n_lat,n_lon/),"float",0);

do ez = 0,n_lev
	do ey = 0,n_lat
		do ex = 0,n_lon
			do j = 0,n_t
				Depar(ez,ey,ex) = Depar(ez,ey,ex)+ (ensemble_meanT(j,ey,ex,ez)-obsT(j,ez,ey,ex))^2;
				Bias(ez,ey,ex) = Bias(ez,ey,ex) + ensemble_meanT(j,ey,ex,ez)-obsT(j,ez,ey,ex);
				DeparU(ez,ey,ex) = DeparU(ez,ey,ex)+ (ensemble_meanU(j,ey,ex,ez)-obsU(j,ez,ey,ex))^2;
				BiasU(ez,ey,ex) = BiasU(ez,ey,ex) + ensemble_meanU(j,ey,ex,ez)-obsU(j,ez,ey,ex);
				do i = 0,n_ensemble
					EnsVar(ez,ey,ex) = EnsVar(ez,ey,ex) + (ensemble_memberT(j,i,ey,ex,ez) - ensemble_meanT(j,ey,ex,ez))^2;
					EnsVarU(ez,ey,ex) = EnsVarU(ez,ey,ex) + (ensemble_memberU(j,i,ey,ex,ez) - ensemble_meanU(j,ey,ex,ez))^2;
				end do
		
			end do
			Depar(ez,ey,ex) = Depar(ez,ey,ex)/(n_t-1)
			Bias(ez,ey,ex) = Bias(ez,ey,ex)^2/(n_t^2-n_t)
			EnsVar(ez,ey,ex) = EnsVar(ez,ey,ex)*(n_ensemble+1)/(n_ensemble*nt*(n_ensemble-1))
			DeparU(ez,ey,ex) = DeparU(ez,ey,ex)/(n_t-1)
			BiasU(ez,ey,ex) = BiasU(ez,ey,ex)^2/(n_t^2-n_t)
			EnsVarU(ez,ey,ex) = EnsVarU(ez,ey,ex)*(n_ensemble+1)/(n_ensemble*nt*(n_ensemble-1))
		end do
	end do
end do
Error_Resid = Depar - Bias - EnsVar
Error_ResidU = DeparU - BiasU - EnsVarU

 
plot = new(4,"graphic")
wks = gsn_open_wks("ps","Tropopause_Temp")
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
wks1 = gsn_open_wks("ps","Mid_Trop_Temp")
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
wks2 = gsn_open_wks("ps","Lower_Trop_Temp")
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
wks3 = gsn_open_wks("ps","Tropopause_U")
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
wks4 = gsn_open_wks("ps","Mid_Trop_U")
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
wks5 = gsn_open_wks("ps","Lower_Trop_U")
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

end
