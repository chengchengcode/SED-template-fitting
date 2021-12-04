c=299792458d2 	; cm/s
AB_zpt_Jy = 3631.																; in unit of Jy

AB_zpt_for_flambda = AB_zpt_Jy / (1.d23 * 1.d-8 / c)

path_to_template = './Brown_sed/'
path_to_filter = './filters/'


cd, path_to_template
spawn, 'ls *.dat', template_name
cd, '..'

;readcol, '../galaxy_name.txt', template_name, format = 'a'
;readcol, '/Users/chengcheng/Jobs/DATA/sed/605/Brown_temp_210317.class.txt', name, a1,a2,a3,typecode, format = 'a'

lambda_filter_obs_name = ['sdssu.filter', 'sdssg.filter', 'sdssr.filter', 'sdssi.filter', 'sdssz.filter', 'J.pb', 'H.pb', 'Ks.pb', 'WISE_WISE.W1.dat', 'WISE_WISE.W2.dat', 'WISE-WISE.W3.dat', 'WISE-WISE.W4.dat']
lambda_filter_obs = fltarr(n_elements(lambda_filter_obs_name))

for i = 0, n_elements(lambda_filter_obs_name) - 1 do begin
	readcol, path_to_filter+lambda_filter_obs_name[i], l_filter, resp_filter, /silen
	lambda_filter_obs[i] = total(l_filter * resp_filter) / total(resp_filter)
endfor

tab = mrdfits('example-sed.fits',1)

redshift = tab.zspec

chi_sq_all = dblarr(n_elements(redshift))

openw, lun_i_xy, 'SED-fitting-chi2list.txt', /get_lun
for i_obs = 0LL, n_elements(redshift) - 1 do begin

;load the sed
mag_obs_all = [tab[i_obs].sdss_u_mag,tab[i_obs].sdss_g_mag,tab[i_obs].sdss_r_mag,tab[i_obs].sdss_i_mag,tab[i_obs].sdss_z_mag,tab[i_obs].J_mag,tab[i_obs].H_mag,tab[i_obs].K_mag,tab[i_obs].w1_mag,tab[i_obs].w2_mag,tab[i_obs].w3_mag,tab[i_obs].w4_mag]
magerr_obs_all = [tab[i_obs].sdss_u_mag_err,tab[i_obs].sdss_g_mag_err,tab[i_obs].sdss_r_mag_err,tab[i_obs].sdss_i_mag_err,tab[i_obs].sdss_z_mag_err,tab[i_obs].J_mag_err,tab[i_obs].H_mag_err,tab[i_obs].K_mag_err,tab[i_obs].w1_mag_err+0.1,tab[i_obs].w2_mag_err+0.1,tab[i_obs].w3_mag_err,tab[i_obs].w4_mag_err]; I add 0.1 to the w1, w2 error since the catalog value seems too low

;remove the band with NAN or -99 mag, or even the low S/N band
ind_filter = where(910. * (1 + redshift[i_obs]) lt lambda_filter_obs and mag_obs_all eq mag_obs_all and mag_obs_all gt 0 and magerr_obs_all gt 0 and mag_obs_all lt 90)

;continue if the SED is bad
if ind_filter[0] eq -1 then begin
	printf, lun_i_xy, i_obs, -1., -999., ' 	FLUX_-99', format = '(i, f40.4, f50.4, a)'
	continue
endif

if n_elements(ind_filter) lt 3 then begin
	printf, lun_i_xy, i_obs, -1., -999., '	fit_bands_lt_3', format = '(i, f40.4, f50.4, a)'
	continue
endif


mag_obs = mag_obs_ALL[ind_filter]
magerr_obs = magerr_obs_ALL[ind_filter]
lambda_template_normal = lambda_filter_obs[ind_filter]

yrange = [max(mag_obs[where(mag_obs gt 0 and mag_obs lt 90)])+2, min(mag_obs[where(mag_obs gt 0 and mag_obs lt 90)])-2]

chi_sq_list = fltarr(n_elements(template_name))
normalize_factor_list = fltarr(n_elements(template_name))

mag_template_normal = fltarr(n_elements(ind_filter))
mag_template_normal_all = fltarr(n_elements(ind_filter), n_elements(template_name))


for i_temp = 0, n_elements(template_name) - 1 do begin

	readcol, 'Brown_sed/'+template_name[i_temp], lambda, flux_lambda, /sil

;	CCM_UNRED, lambda, flux_lambda, ebv, funred;, [ R_V = ]      
;	flux_lambda = funred

	lambda_z = lambda * (1. + redshift[i_obs])
	flux_lambda_z = flux_lambda / (1. + redshift[i_obs]) 									; in unit of erg/s/cm2/A

	;flux_nu = 1.d23 * flux_lambda_z * lambda_z * lambda_z * 1.d-8 / c				; in unit of Jy
	flux_nu = 3.34d4 * flux_lambda_z * lambda_z * lambda_z							; in unit of Jy, from AB mag wiki

    nu = 3.d18 / lambda_z


	;photometry to the template model with the filters
	;cgplot, lambda_z, - 2.5*alog10(flux_nu / 3631.), yran = [18, 8], xran = [2000, 250000], /xlog, xtitl = textoidl('Wavelength [\AA]'), ytitl = 'AB mag'	
	for i_filter = 0, n_elements(ind_filter) - 1 do begin
		readcol, path_to_filter+lambda_filter_obs_name[ind_filter[i_filter]], l_filter, resp_filter, /silen
	    min_l = min(l_filter)
	    max_l = max(l_filter)    
		resp = fltarr(n_elements(lambda_z))
		filter_range=where( lambda_z gt min_l and lambda_z lt max_l)
		resp[filter_range] = interpol(resp_filter, l_filter, lambda_z[filter_range])
		mag_template_normal[i_filter] = -2.5*alog10 ( int_tabulated(nu[filter_range], flux_nu[filter_range] * resp[filter_range] / nu[filter_range] ,/sort) / int_tabulated(nu[filter_range], resp[filter_range] * AB_zpt_Jy / nu[filter_range] ,/sort) )
		;if one would like to integrate across the lambda to get the photometry of the template model, then:
		;mag_template_normal[i_filter] = -2.5 * alog10 (int_tabulated(lambda_z[filter_range], flux_lambda_z[filter_range] * lambda_z[filter_range] * resp[filter_range], /sort) / int_tabulated(lambda_z[filter_range], resp[filter_range] * AB_zpt_for_flambda / lambda_z[filter_range] ,/sort) 	)

	endfor
	;cgoplot,lambda_filter_obs[ind_filter], mag_template_normal, ps = 9, color = 'red'	;check the template magnitude.
	;stop

	mag_template_normal_all[*, i_temp] = mag_template_normal

	; derive the nomilize factor and chi2 of this template
	normalize_factor_list[i_temp] = total((mag_obs - mag_template_normal)/magerr_obs^2)/total(1./magerr_obs^2)
	chi_sq_list[i_temp] = total((mag_obs - mag_template_normal - normalize_factor_list[i_temp])^2. / magerr_obs^2.) / n_elements(ind_filter)

	;show the current SED and the fitting template:
	plot, lambda_template_normal, mag_obs, /xlog, xrange = [min(lambda_filter_obs)-1000, max(lambda_filter_obs)*1.5], yrange = yrange, linesty = 1, xsty = 1, ysty = 1, charsize = 2, charthick = 1.5, /nod, xtitle = textoidl('Wavelength [\AA]'), ytitle = 'AB mag'
	oploterror, lambda_template_normal, mag_obs , magerr_obs, psy = symcat(16), color = cgcolor('green')
	oplot, lambda, - 2.5*alog10(flux_nu / AB_zpt_Jy) + normalize_factor_list[i_temp], ps = 3
;	cgoplot, lambda_filter_obs[ind_filter], mag_template_normal + normalize_factor_list[i_temp], color = 'red', ps = 16

	endfor

	;find the template with the minimal chi2
	chi_sq_all[i_obs] = min(chi_sq_list)

	normalize_factor = normalize_factor_list[where(chi_sq_list eq min(chi_sq_list))]
	normalize_factor = normalize_factor[0]

	;then load this template, and show the best fitting results
	readcol, path_to_template+template_name[where(chi_sq_list eq min(chi_sq_list))], lambda, flux_lambda, /sil
	lambda_z = lambda * (1. + redshift[i_obs])
	flux_lambda_z = flux_lambda / (1. + redshift[i_obs]) 									; in unit of erg/s/cm2/A
	flux_nu = 3.34d4 * flux_lambda_z * lambda_z * lambda_z									; in unit of Jy, from AB mag wiki
	lambda_z = lambda * (1. + redshift[i_obs])
		
	plot, lambda_template_normal, mag_obs, /xlog, xrange = [min(lambda_filter_obs)-1000, max(lambda_filter_obs)*1.5], yrange = yrange, xsty = 1, ysty = 1, charsize = 2, charthick = 1.5, /nod, xtitle = textoidl('Wavelength [\AA]'), ytitle = 'AB mag'
	oplot, lambda_z, - 2.5*alog10(flux_nu / AB_zpt_Jy) + normalize_factor, color = cgcolor('orange');, ps = 3
	oploterror, lambda_template_normal, mag_obs , magerr_obs, psy = symcat(16), color = cgcolor('green')
	oplot, lambda_template_normal, mag_template_normal_all[*, where(chi_sq_list eq min(chi_sq_list))] + normalize_factor, psy = symcat(9), thick = 1, symsize = 2, color = cgcolor('red')
	cgtext, 300,140, 'spec_z = '+strtrim(redshift[i_obs],2), color = cgcolor('red'),charthick = 1, charsize = 2,/device
	cgtext, 300,120, 'chi2 = '+ strtrim(chi_sq_list[where(chi_sq_list eq min(chi_sq_list))],2), color = cgcolor('cyan'),charthick = 1, charsize = 2,/device

	;save the id, specz, minimal chi2, and the template name for this target.
	printf, lun_i_xy, [i_obs], redshift[i_obs], chi_sq_list[where(chi_sq_list eq min(chi_sq_list))], '  '+template_name[[where(chi_sq_list eq min(chi_sq_list))]], format = '(i, f40.4, f50.4, a)'
	flush, lun_i_xy
wait, 1
;cc_pause
endfor
free_lun, lun_i_xy






end




