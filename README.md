# Detección de Metástasis en Cáncer de Mama y Cáncer de Próstata mediante Machine Learning a partir de RNA-seq
El objetivo de este trabajo es crear modelos de machine learning de metástasis a partir de datos de expresión génica, con la finalidad de detectar si una persona tiene metástasis o no,  evaluando su rendimiento mediante diversas métricas y herramientas, y comparando sus resultados para identificar el enfoque más eficaz.

### Obtención de Datos

Los datos con los que se trabaja no están disponibles en el repositorio debido a su tamaño, pero se pueden conseguir ejecutando los archivos `PRAD-data-cleaning.R` y `BRCA-data-cleaning.R`.

#### Pasos para Obtener los Datos:

1. **Archivos de Extracción de Datos**: Los datos no están incluidos directamente en el repositorio debido a su gran tamaño. En su lugar, se proporcionan scripts de extracción de datos que se deben ejecutar para obtener los datos procesados. Estos scripts son:
   - `PRAD-data-cleaning.R`
   - `BRCA-data-cleaning.R`

2. **Ejecutar los Scripts**: Para obtener los datos, simplemente hay que ejecutar los scripts mencionados, adaptando el nombre del directorio base. Estos scripts están diseñados para:
   - Descargar los datos brutos.
   - Guardar los datos en el formato requerido para el análisis posterior.

3. **Resultado**: Al finalizar la ejecución de los scripts, los conjuntos de datos están listos para ser utilizados por los archivos de Jupyter Notebook.
