from atomworks.io.utils.ccd import get_available_ccd_codes_in_biotite
metals = ["LA","ND","SM","GD","TB","Y","CE","PR","PM","EU","DY","HO","ER","TM","YB","LU","SC"]
available = get_available_ccd_codes_in_biotite()
for m in metals:
    print(f"{m} available: {m in available}")
