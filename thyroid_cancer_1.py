#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas, numpy, seaborn
import sklearn, sklearn.decomposition, sklearn.decomposition, sklearn.pipeline, sklearn.preprocessing
import scipy, scipy.signal
import matplotlib, matplotlib.pyplot
matplotlib.rcParams.update({'font.size':20, 'font.family':'FreeSans', 'xtick.labelsize':30, 'ytick.labelsize':30, 'axes.labelsize':40, 'figure.figsize':(12, 8)})
import cycler
import pandas as pd
import re 
import requests
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler


# In[3]:


'''

# Open GDC manifest file.
mf_name = input('Write manifest file name: ') 
data = pd.read_csv(mf_name,delimiter='\t')

# Directory where the files are saved.  
output_directory = "/Users/ASUS/Documents/bioinf/data"
os.makedirs(output_directory, exist_ok=True)

count = 0
for _, row in data.iterrows():
    # Loop through the file IDs and download each file.
    file_id = str(row[0])  
    data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)
    response = requests.get(data_endpt, headers={"Content-Type": "application/json"})

    # Check if the request was successful before proceeding.
    if response.status_code == 200:
        # The file name can be found in the header within the Content-Disposition key.
        response_head_cd = response.headers.get("Content-Disposition")
        if response_head_cd:
            file_name = re.findall("filename=(.+)", response_head_cd)[0]

            # Join the directory path with the file name to create the complete file path.
            output_file_path = os.path.join(output_directory, file_name)

            with open(output_file_path, "wb") as output_file:
                output_file.write(response.content)
                count = count + 1
                print(f"Downloaded file nº {count}: {file_name} " )
        else:
            print(f"Failed to get file name for file ID: {file_id}")
    else:
        print(f"Failed to download file ID: {file_id} - Status Code: {response.status_code}")'''



# In[4]:


#also download metadata and try to associate patient ID's with case ID's


# In[2]:


# Get the list of all files and directories
# in the root directory
path = "/Users/ASUS/Documents/bioinf/data"
dir_list = os.listdir(path)
  
print("Files and directories in '", path, "' :") 
  
# print the list
#print(dir_list)


# In[4]:


#We use a for loop to create a list containing the names of the files in our directory

#for element in dir_list:
    #print(element)


# In[3]:


concat_df=pandas.DataFrame()
print(type(concat_df))
for element in dir_list[:100]: 
    #We first visualize the name of the file 
    #print(element)

    #specify path
    path_to_data="C:/Users/ASUS/Documents/bioinf/data/"
    full_path=path_to_data+element
    

    #Create a dataframe of that file
    df = pandas.read_csv(full_path, sep="\t",skiprows=1, index_col=0)
    df.drop(['gene_name', 'gene_type'], axis=1, inplace = True)
    df.drop(['N_unmapped', 'N_multimapping', 'N_noFeature', 'N_ambiguous'], axis=0, inplace= True)
    
    #specify row and column indexes
    #file_id = str(row[1]).split('.')[0]
    
    #take the first column of each one
    expression_column = df['tpm_unstranded']
    expression_column = expression_column.rename(element.split('.')[0])   ######################
    
    if concat_df.shape == (0,0):
        concat_df = expression_column 
        print('starting to construct the concatenated dataframe')
    else:
        concat_df = pandas.concat([concat_df, expression_column], axis=1)


# In[4]:


df


# In[5]:


concat_df


# In[6]:


#convert df to csv.
concat_df.to_csv('/Users/ASUS/Documents/bioinf/tpm.csv')


# In[7]:


tpm = pandas.read_csv('/Users/ASUS/Documents/bioinf/tpm.csv', index_col=0)
print(tpm.shape)
tpm.head()


# In[8]:


#filter out genes that (almost) never show up   ##############################################################################
expressed_genes = tpm[tpm.max(axis=1) >= 8]
expressed_genes.shape


# In[9]:


#bring data to log2 TPM
log2_tpm_PO = numpy.log2(expressed_genes + 1)
log2_tpm_PO


# In[10]:


'''HISTOGRAM'''

found_max = 8 # especifico el valor màxim aqui. El poso una mica mes alt
number_of_bins = found_max*10
print(number_of_bins)

absolute_max = 0  #defineixo absolute_max = 0
working_samples = log2_tpm_PO.columns.to_list()  #afegeixo les columnes amb les que vull treballar a una llista working_samples

most_likely_expressions = []
all_hats = []
for sample in working_samples:
    
    log2TPM = log2_tpm_PO.loc[:, sample] #selecciono una columna i la guardo a log2TPM
    if max(log2TPM) > absolute_max:   #si el valor maxim de log2TPM és major que l'absolut_max, absolut_max passa a ser aquest valor
        absolute_max = max(log2TPM)
       
    #print(numpy.min(log2TPM), numpy.max(log2TPM))
    
    # agafa el valor de dalt del mig de cada barra
    # agafa els 2 valors dels extrems de baix de cada barra i resta i divideix entre 2 per obtenir un unic valor de baix 
    hist, bin_edges = numpy.histogram(log2TPM, bins=number_of_bins, range=(0, found_max))
    half_bin = (bin_edges[1] - bin_edges[0])/2
    x = bin_edges + half_bin   
    x = x[:-1]   #elimina l'ultim valor pq no el necessitem
    
    #print(x)   #--> fem un 'print' per veure les dades i agafar el segon pic
    #print(hist)
    
    plotting_x = x#[1:500] #depen de com surti lhistograma agafem a pertir de unes dades o altres 
    plotting_hist = hist
    #print(plotting_x)
    
    #matplotlib.pyplot.plot(plotting_x, plotting_hist,'.', alpha=1/3)
    yhat = scipy.signal.savgol_filter(plotting_hist, 51, 3)
    matplotlib.pyplot.plot(plotting_x, yhat, '-', lw=4, alpha = 1/2)
    
    most_likely_expression = x[numpy.argmax(yhat)]
    most_likely_expressions.append(most_likely_expression)
    
    all_hats.append(yhat)
    
    
matplotlib.pyplot.xlim([numpy.min(plotting_x)+0.75, numpy.max(plotting_x)+0.25])
matplotlib.pyplot.ylim([0, 400])

matplotlib.pyplot.xlabel('log2 (TPM+1)')
matplotlib.pyplot.ylabel('Gene count')
matplotlib.pyplot.grid(ls=':')
#matplotlib.pyplot.legend(ncol=10, fontsize=12, bbox_to_anchor=(1.02, 1.25))

matplotlib.pyplot.tight_layout()

print(absolute_max)
    


# In[11]:


#mirem quin és el vallor màxim que tenim abans de fer l'histograma anterior per saber quins paràmetres hem de posar
maxim = max(x)
print(maxim)


# In[12]:


hist, bin_edges = numpy.histogram(most_likely_expressions, bins=20, range=(2, 6))
half_bin = (bin_edges[1] - bin_edges[0])/2
x = bin_edges + half_bin
x = x[:-1]

matplotlib.pyplot.plot(x, hist, '-', lw=4, alpha=1/2)
matplotlib.pyplot.grid(ls=':')
matplotlib.pyplot.xlabel('Most expected log2 TPM + 1')
matplotlib.pyplot.ylabel('Sample count')
matplotlib.pyplot.tight_layout()


# In[13]:


median = numpy.median(most_likely_expressions)
std = numpy.std(most_likely_expressions)
threshold = median - 1.96*std
print(median, std, threshold)
print()

new_threshold = 3.5
suspicious_samples = []
for i in range(len(most_likely_expressions)):
    if most_likely_expressions[i] <= new_threshold:
        print(working_samples[i], new_threshold, most_likely_expressions[i])
        suspicious_samples.append(working_samples[i])
print(suspicious_samples)


# In[14]:


for i in range(len(working_samples)):
    if working_samples[i] in suspicious_samples:
        matplotlib.pyplot.plot(plotting_x, all_hats[i], '-', lw=4, alpha=1/3, color='red')
    else:
        matplotlib.pyplot.plot(plotting_x, all_hats[i], '-', lw=4, alpha=1/3, color='blue')
    
matplotlib.pyplot.xlim([numpy.min(plotting_x)-0.25, numpy.max(plotting_x)+0.25])
matplotlib.pyplot.ylim([0, 400])

matplotlib.pyplot.xlabel('log2 (TPM+1)')
matplotlib.pyplot.ylabel('Gene count')
matplotlib.pyplot.grid(ls=':')
#matplotlib.pyplot.legend(ncol=10, fontsize=12, bbox_to_anchor=(1.02, 1.25))

matplotlib.pyplot.tight_layout()


# In[15]:


#We use '.transpose()' function to change df format so we obtain values that were in columns now in rows and viceversa
transpose = log2_tpm_PO.transpose()
print(transpose.shape)
transpose.head()


# In[16]:


from scipy.stats import zscore

# We have a DataFrame called 'log2_tpm_PO' with gene expresion data and we want to distinguish the top_500 genes frfom the 
#other ones. Top 500 genes are the ones with higher mean and higher variance. 

# Calcula la media y la varianza de expresión para cada gen
average_expression = log2_tpm_PO.mean(axis=1)
variance_by_rows = log2_tpm_PO.var(axis=1)

# Calcula los puntajes Z para las columnas 'Average' y 'Variance'
average_zscores = zscore(average_expression)
variance_zscores = zscore(variance_by_rows)

# Calcula el puntaje de bondad sumando los puntajes Z de media y varianza
goodness_score = average_zscores + variance_zscores

# Crea un DataFrame con los puntajes de bondad y los nombres de los genes
genes_exp = pd.DataFrame({'Goodness_Score': goodness_score, 'Gene_Name': log2_tpm_PO.index})

# Selecciona los 500 genes con los puntajes de bondad más altos
top_genes = genes_exp.nlargest(500, 'Goodness_Score')

#Select the 50 top genes 
top_50_genes = genes_exp.nlargest(50, 'Goodness_Score')


# Muestra los primeros 5 genes en la lista de los mejores 500 genes
print(top_genes.head(5))


# In[72]:


'''scatterplot with standard deviation'''

# Calcular la media y la desviación estándar para cada gen
media = transpose.mean()
std_dev = transpose.std()

# Crear un DataFrame con las medias y desviaciones estándar
data_for_scatter = pd.DataFrame({'Media': media, 'Desviación Estándar': std_dev})

# Crear un scatter plot
plt.figure(figsize=(10, 8))

mida=10

# Crear un scatter plot con todos los genes en azul y forma circular
sns.scatterplot(x='Media', y='Desviación Estándar', data=data_for_scatter, color='blue', alpha=1/2, label='Otros Genes')

# Crear un scatter plot con los 500 mejores genes en rojo y forma cuadrada
sns.scatterplot(x='Media', y='Desviación Estándar', data=data_for_scatter.loc[top_genes.index], color='red', marker='D', alpha=1/2, label='Top 500 Genes')

# Personalizar el scatter plot
plt.title('Scatter Plot de Media vs. Standard Deviation')
plt.xlabel('Media')
plt.ylabel('SD')

# Mostrar el scatter plot con una leyenda
plt.legend()

# Mostrar el scatter plot
plt.show()


# In[71]:


'''SCATTERPLOT WITH COEFICIENT OF VARIANCE'''

# Calcular la media y la desviación estándar para cada gen
media = transpose.mean()
std_dev = transpose.std()
cv = (std_dev / media) * 100  # Calcula el CV y multiplica por 100 para obtener un porcentaje

# Crear un DataFrame con las medias y desviaciones estándar
data_for_scatter = pd.DataFrame({'Media': media, 'Coeficient of variance': cv})

# Crear un scatter plot
plt.figure(figsize=(10, 8))

mida=10

# Crear un scatter plot con todos los genes en azul y forma circular
sns.scatterplot(x='Media', y='Coeficient of variance', data=data_for_scatter, color='blue', alpha=1/2, label='Otros Genes')

# Crear un scatter plot con los 500 mejores genes en rojo y forma cuadrada
sns.scatterplot(x='Media', y='Coeficient of variance', data=data_for_scatter.loc[top_genes.index], color='red', marker='D', alpha=1/2, label='Top 500 Genes')

# Personalizar el scatter plot
plt.title('Scatter Plot de Media vs. Coeficient of variance')
plt.xlabel('Media')
plt.ylabel('CV')

# Mostrar el scatter plot con una leyenda
plt.legend()

# Mostrar el scatter plot
plt.show()


'''fix a threshold '''


# In[35]:


import pandas as pd

# Supongamos que 'top_genes' contiene los nombres de los genes que deseas mantener
important_genes = top_50_genes['Gene_Name'].tolist()

# Filtra las columnas de 'transpose' para mantener solo los genes presentes en 'top_genes'
important_dataframe = transpose[important_genes]
important_dataframe

############################

important_500_genes = top_genes['Gene_Name'].tolist()

important_500_dataframe = transpose[important_500_genes]
important_500_dataframe


# In[34]:


normalized_data


# In[20]:


'''This boxplot shows gene expression  data that has been normalized through Z-score standarizationusing Seaborn
and Matplotlib'''


scaler = StandardScaler()
normalized_data = pd.DataFrame(scaler.fit_transform(important_dataframe), columns=important_dataframe.columns)

# Crea un boxplot para la gene expression data normalizada
plt.figure(figsize=(15, 8))
sns.boxplot(data=normalized_data, orient="v", palette="Set2")

# Personaliza el gráfico
plt.ylabel('Z score')
plt.xlabel('Genes')
plt.title('Boxplot of Top 50 genes expression')

# Puedes agregar plt.xticks([]) si deseas ocultar las etiquetas del eje x
plt.xticks([])

# Muestra el gráfico
plt.show()


# In[40]:


scaler = StandardScaler()
normalized_500 = pd.DataFrame(scaler.fit_transform(important_500_dataframe), columns=important_500_dataframe.columns)
normalized_500


# In[21]:


from sklearn.decomposition import PCA
n_componentes = 2  # Elige el número de componentes principales que deseas
pca = PCA(n_components=n_componentes)
componentes_principales = pca.fit_transform(transpose)

plt.figure(figsize=(8, 6))
plt.scatter(componentes_principales[:, 0], componentes_principales[:, 1])
plt.xlabel('PC1')  
plt.ylabel('PC2')  
plt.title('PCA de Datos de Expresión Génica')
plt.show()



# In[22]:


#remove low quality samples 
strict = transpose.drop(index=['00f54652-691e-4446-869a-4dd51c236f56', '036f7712-a636-4198-b568-6c85d9db4a3b', '0386037f-728c-47b0-93d0-6bb81b8b2a05', '0cfde228-89bc-4454-a894-d6fc4b032892', '144be06b-54da-4b42-a219-8bbe88cc7f78', '153946fb-8375-4d84-ad2f-2988722b845e', '164fe0cc-c07d-4d09-adfd-cca86b737aa3', '185ecfd0-680c-4874-8667-5d7543ec562c', '193da121-43e7-40ab-8e0e-980c62586568', '1c4e7140-e649-4267-92ac-f30d65d795e5', '1e100c4d-13eb-4e5a-9117-ead385902710', '1e31f670-de56-4122-8966-42767135b420', '1ff78fb6-b27a-4946-b39f-c87d6dfa669d', '20c70125-f4ec-4d52-9efd-a207b559fedc', '2231c28a-bf18-4c0c-885c-42896df3b0e8', '284b2aab-f424-43d3-b160-78f2c28ecf94', '28da87d8-0ac9-4a62-9591-67a713bfe4e7', '28ecaa28-d2dc-4988-a484-544d4ad7355d', '2a51bda3-fb5e-4eb6-bc86-e78dca8c4e2b', '2a90f5ad-a9d0-4f5d-a255-3e8233039868', '2b7a660d-f472-476a-9145-cae4f3ba5e4c', '2df38eb8-5350-4951-9159-a8add6474efe', '2ea7afdd-f947-47e6-9516-b6f38e0a1967', '30b75646-82f8-467f-8f2b-b71b97af3bc6', '3682a3e2-17ed-4f6c-b4de-9057e6d7f4c6', '3762f95e-180b-43c2-967c-0e138fb23b63']) #canviar el nom de la mostra que cal eliminar 
print(strict.shape)
strict.head()


'''PORTAR A DALT'''


# In[23]:


features = strict.columns
x = strict.loc[:, features].values
 
x = sklearn.preprocessing.StandardScaler().fit_transform(x)
pca = sklearn.decomposition.PCA(n_components=2)
principalComponents = pca.fit_transform(x)

principalDf = pandas.DataFrame(data = principalComponents, columns = ['PC1', 'PC2'])
x


# In[24]:


pca.explained_variance_ratio_


# In[30]:


'''HEATMAP OF GENE EXPRESSION   ---> ARREGLAR, ferlo sobre el z-score'''

plt.figure(figsize=(12, 8))  # Set the figure size as needed
sns.heatmap(log2_tpm_PO, cmap='coolwarm', annot=False, cbar=True)

# Customize the plot
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.title('Heatmap of Gene Expression')

plt.xticks([])
plt.yticks([])

# Show the plot
plt.show()


# In[53]:


'''Z-SCORE HEATMAP'''


plt.figure(figsize=(12, 8))  # Set the figure size as needed
sns.heatmap(normalized_500, cmap='bwr', annot=False, cbar=True)

# Customize the plot
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.title('Heatmap of top 500 Gene Expression')

plt.xticks([])
plt.yticks([])

# Show the plot
plt.show()


# In[60]:


linkage_method = 'complete'
distance_metric = 'cosine'

seaborn.set(font_scale=0.9)
seaborn.clustermap(normalized_500, cmap='bwr', col_cluster=False, vmin=-1, vmax=1, method=linkage_method, metric=distance_metric, yticklabels=[], xticklabels = [], cbar_kws={'label':'log2FC'})


matplotlib.pyplot.title('{} {}'.format(linkage_method, distance_metric))
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.show()



# In[73]:


import seaborn as sns
import matplotlib.pyplot as plt

linkage_method = 'complete'
distance_metric = 'cosine'

sns.set(font_scale=0.9)
sns.clustermap(
    normalized_500,  # Tus datos normalizados
    cmap='bwr',   # Mapa de colores
    col_cluster=True,  # Realizar clustering en el eje de las columnas
    row_cluster=True,  # Realizar clustering en el eje de las filas
    vmin=-1,  # Valor mínimo en el mapa de colores
    vmax=1,   # Valor máximo en el mapa de colores
    method=linkage_method,  # Método de enlace para el clustering
    metric=distance_metric,  # Métrica de distancia
    yticklabels=[],  # Mostrar etiquetas en el eje y
    xticklabels=[],  # Etiquetas en el eje x
    cbar_kws= {'label': 'log2FC'}  # Etiqueta para la barra de color

)



plt.title('{} {}'.format(linkage_method, distance_metric))
plt.tight_layout()
plt.show()


# In[29]:


'''HEATMAP of top 50 gene expresion'''

plt.figure(figsize=(12, 8))  # Set the figure size as needed
sns.heatmap(normalized_data, cmap='coolwarm', annot=False, cbar=True)

# Customize the plot
plt.xlabel('Samples')
plt.ylabel('Genes')
plt.title('Heatmap of Gene Expression')

plt.xticks([])
plt.yticks([])

# Show the plot
plt.show()


# In[ ]:




