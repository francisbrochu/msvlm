import msvlm.msAlign.msAlign as ms
import csv 

list_of_spectra = []

with open('alignment.csv', newline='') as f:
    reader = csv.reader(f)
    for row in reader:
        list_of_strings = list(row)
        spectrum = []
        for s in list_of_strings:
            spectrum.append(float(s))
        list_of_spectra.append(spectrum)

vlm_pts = ms.find_vlm(list_of_spectra, 5e-6)

align_pts = ms.find_alpt(list_of_spectra, 5e-6)

with open('vlmPoints.txt','w') as fich_vlm:
    for x in vlm_pts:
        print(x,file=fich_vlm)

with open('alignmentPoints.txt','w') as fich_alp:
    for x in align_pts:
        print(x,file=fich_alp)

