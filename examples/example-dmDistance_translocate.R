### get a distance matrix from a list of dmProfile while minimizing distance by translocating
somedmProfile <-  alldmProfile[1:10] 
dblist <- dmDistance_translocate(dmprofileList = somedmProfile, slide = 500)
