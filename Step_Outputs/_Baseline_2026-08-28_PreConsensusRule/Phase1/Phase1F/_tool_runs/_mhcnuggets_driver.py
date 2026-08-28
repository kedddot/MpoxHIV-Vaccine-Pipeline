
import sys, csv
from mhcnuggets.src.predict import predict
peptides_path, allele, output_path = sys.argv[1], sys.argv[2], sys.argv[3]
predict(class_='II', peptides_path=peptides_path, mhc=allele, output=output_path, rank_output=True)
