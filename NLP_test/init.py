import os
import pandas as pd

filepath = "C:/Users/mdnur/OneDrive/Desktop/python/LipidNet/NLP_test"
filenames = []
labels =[]
for folder in os.listdir(filepath):
	folpath = os.path.join(filepath, folder)
	print(folpath)
	for file in os.listdir(folpath):
		if file.endswith('.txt'):
			fp = os.path.join(folpath, file)
			filenames.append(fp)
			labels.append(folder)
df = pd.DataFrame({'filename':filenames, 'label':labels})
df.to_csv('filemap.csv')
		
