#!/bin/csh

# WAVES
#foreach i ( `awk 'BEGIN { for (k=1996; k<=2022; k=k+1) print k; exit }'` )
#    echo $i
#    foreach j ( `awk 'BEGIN { for (k=1; k<=12; k=k+1) print k; exit }'` )
#	if($j<10)then
#	    set name=0$j
#	else
#	    set name=$j
#	endif
#    wget -r -l1 -nd -nc -np -e robots=off -A.nc --no-check-certificate https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waves/waverider-17m/$i/FRF-ocean_waves_waverider-17m_$i$name.nc
#    end
#end

# DEMs SURVEY
#foreach i ( `awk 'BEGIN { for (k=2008; k<=2022; k=k+1) print k; exit }'` )
#    foreach m ( `awk 'BEGIN { for (k=1; k<=12; k=k+1) print k; exit }'` )
#	if($m<10)then
#	    set mois=0$m
#	else
#	    set mois=$m
#	endif
#	foreach j ( `awk 'BEGIN { for (k=1; k<=31; k=k+1) print k; exit }'` )
#	    if($j<10)then
#		set jour=0$j
#	    else
#		set jour=$j
#	    endif
#	wget -r -l1 -nd -nc -np -e robots=off -A.nc --no-check-certificate https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/geomorphology/DEMs/surveyDEM/data/FRF_geomorphology_DEMs_surveyDEM_$i$mois$jour.nc
#	end
#    end
#end

# ELEVATION TRANSECTS
#foreach i ( `awk 'BEGIN { for (k=1974; k<=2022; k=k+1) print k; exit }'` )
#    foreach m ( `awk 'BEGIN { for (k=1; k<=12; k=k+1) print k; exit }'` )
#	if($m<10)then
#	    set mois=0$m
#	else
#	    set mois=$m
#	endif
#	foreach j ( `awk 'BEGIN { for (k=1; k<=31; k=k+1) print k; exit }'` )
#	    if($j<10)then
#		set jour=0$j
#	    else
#		set jour=$j
#	    endif
#	wget -r -l1 -nd -nc -np -e robots=off -A.nc --no-check-certificate https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/geomorphology/elevationTransects/survey/data/FRF_geomorphology_elevationTransects_survey_$i$mois$jour.nc
#	end
#    end
#end

# TIDE
foreach i ( `awk 'BEGIN { for (k=1981; k<=2022; k=k+1) print k; exit }'` )
    echo $i
    foreach j ( `awk 'BEGIN { for (k=1; k<=12; k=k+1) print k; exit }'` )
	if($j<10)then
	    set name=0$j
	else
	    set name=$j
	endif
    wget -r -l1 -nd -nc -np -e robots=off -A.nc --no-check-certificate https://chlthredds.erdc.dren.mil/thredds/fileServer/frf/oceanography/waterlevel/eopNoaaTide/$i/FRF-ocean_waterlevel_eopNoaaTide_$i$name.nc
    end
end	    
